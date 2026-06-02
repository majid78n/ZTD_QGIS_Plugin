"""Iterative Tropospheric Decomposition (ITD) and its Augmented variant (AITD).

This module implements the spatial-interpolation strategy described in:

    Song, Venuti, Monti-Guarnieri, Manzoni (2025),
    "Augmented iterative tropospheric decomposition strategy for GNSS-based
    zenith tropospheric delay map generation",
    Environmental Modelling and Software 194, 106669.

The ZTD at a point P(x, h) is modelled as the sum of
  * a *stratified* component S(h) = L * exp(-B * h)   (height dependent, Eq. 2)
  * a *turbulent* component  T(x)                      (horizontally correlated)
  * noise.

Two predictors are provided, both producing a ZTD map, a stratified map and a
turbulent map sampled on a DEM grid:

  ITD  (Sec. 2.1) - the iterative decomposition is repeated for every DEM knot,
                    using the classic inverse-distance weight 1/d^2 (Eq. 5).
  AITD (Sec. 2.2) - a two-step procedure: the decomposition is run once per
                    observation point, then the stratified parameters (L, B) and
                    the turbulent residuals are interpolated onto the DEM with the
                    Modified IDW (MIDW) weight 1/(1 + d^2/alpha) (Eq. 8).

All horizontal coordinates passed in MUST be in a metric CRS (e.g. UTM): the
decorrelation radius and all distances are expressed in metres. Heights are
scaled to [0, 1] *before* being passed in (see ``scale_heights``).

The code is intentionally free of any PyQGIS dependency so that it can be unit
tested in a plain NumPy environment.
"""

from __future__ import annotations

import numpy as np

# Optional SciPy acceleration for neighbour search. Falls back to NumPy.
try:  # pragma: no cover - exercised implicitly depending on environment
    from scipy.spatial import cKDTree as _KDTree
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover
    _KDTree = None
    _HAVE_SCIPY = False


# --------------------------------------------------------------------------- #
# Helpers: height scaling and stratified-model fitting
# --------------------------------------------------------------------------- #
def scale_heights(h, h_min=None, h_max=None):
    """Scale heights to the [0, 1] range used by the stratified model (Eq. 2).

    Returns ``(scaled, h_min, h_max)``. Pass the returned ``h_min``/``h_max``
    back in when scaling a second array (e.g. the DEM) so both share one scale.
    """
    h = np.asarray(h, dtype=float)
    if h_min is None:
        h_min = float(np.nanmin(h))
    if h_max is None:
        h_max = float(np.nanmax(h))
    span = h_max - h_min
    if span <= 0:
        # Degenerate (flat) terrain: avoid division by zero.
        return np.zeros_like(h), h_min, h_max
    return (h - h_min) / span, h_min, h_max


def fit_stratified(h, z, max_iter=15, tol=1e-10):
    """Fit ``z ~= L * exp(-B * h)`` and return ``(L, B)``.

    A log-linear least squares provides a robust initial guess
    (``ln z = ln L - B h``); a few Gauss-Newton steps then refine the genuine
    non-linear least-squares solution. No SciPy dependency required.
    """
    h = np.asarray(h, dtype=float)
    z = np.asarray(z, dtype=float)
    n = z.size
    if n == 0:
        return np.nan, np.nan
    if n == 1:
        # Cannot resolve a slope from one point: flat stratification.
        return float(z[0]), 0.0

    # --- log-linear initialisation -------------------------------------- #
    zpos = np.clip(z, 1e-6, None)
    A = np.column_stack([np.ones_like(h), -h])
    try:
        coef, *_ = np.linalg.lstsq(A, np.log(zpos), rcond=None)
        L = float(np.exp(coef[0]))
        B = float(coef[1])
    except Exception:
        L, B = float(np.mean(z)), 0.0
    if not np.isfinite(L) or not np.isfinite(B):
        L, B = float(np.mean(z)), 0.0

    # --- Gauss-Newton refinement of the non-linear LS ------------------- #
    for _ in range(max_iter):
        e = np.exp(-B * h)
        model = L * e
        r = z - model
        # Jacobian columns: d/dL = e ; d/dB = -L*h*e
        jL = e
        jB = -L * h * e
        J = np.column_stack([jL, jB])
        JTJ = J.T @ J
        JTr = J.T @ r
        try:
            delta = np.linalg.solve(JTJ + 1e-9 * np.eye(2), JTr)
        except np.linalg.LinAlgError:
            break
        L += delta[0]
        B += delta[1]
        if not np.isfinite(L) or not np.isfinite(B):
            L, B = float(np.mean(z)), 0.0
            break
        if np.max(np.abs(delta)) < tol:
            break
    return float(L), float(B)


# --------------------------------------------------------------------------- #
# Weighting functions (Eq. 5 and Eq. 8)
# --------------------------------------------------------------------------- #
def midw_alpha(radius):
    """alpha so the MIDW weight equals 0.01 at the decorrelation radius (Eq. 8).

    1 / (1 + r^2/alpha) = 0.01  =>  alpha = r^2 / 99
    """
    return (radius * radius) / 99.0


def _normalise_rows(w):
    """Row-normalise a weight matrix; rows that sum to zero stay zero."""
    s = w.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return w / s


def _pair_weights(d2, scheme, alpha):
    """Weight matrix from a squared-distance matrix.

    ``scheme='idw'``  -> 1/d^2          (Eq. 5, classic ITD)
    ``scheme='midw'`` -> 1/(1+d^2/alpha) (Eq. 8, AITD)
    The diagonal (self distance, d2=inf or 0) is forced to zero weight.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        if scheme == "idw":
            w = 1.0 / d2
        else:
            w = 1.0 / (1.0 + d2 / alpha)
    w[~np.isfinite(w)] = 0.0
    return w


# --------------------------------------------------------------------------- #
# The shared iterative decomposition (Eqs. 3-7)
# --------------------------------------------------------------------------- #
def decompose_window(h, z, x, y, scheme, alpha, max_iter=10, tol=1e-6):
    """Iteratively split the ZTDs of a set of reference points.

    Parameters
    ----------
    h : scaled heights of the reference points (already in [0, 1]).
    z : observed ZTDs of the reference points.
    x, y : metric horizontal coordinates of the reference points.
    scheme : 'idw' (ITD) or 'midw' (AITD) turbulent weighting.
    alpha : MIDW alpha (ignored for 'idw').

    Returns
    -------
    (L, B, residual) where ``residual = z - L*exp(-B*h)`` is the turbulent
    component at each reference point after convergence.
    """
    h = np.asarray(h, float)
    z = np.asarray(z, float)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = z.size
    if n == 0:
        return np.nan, np.nan, np.array([])
    if n == 1:
        return float(z[0]), 0.0, np.array([0.0])

    # Pairwise squared distances among reference points (self -> inf).
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    d2 = dx * dx + dy * dy
    np.fill_diagonal(d2, np.inf)
    W = _normalise_rows(_pair_weights(d2, scheme, alpha))

    # Iteration: start with zero turbulence (Eq. 3), then alternate
    # stratified fit (Eq. 6) and IDW turbulent estimate (Eq. 7).
    L, B = fit_stratified(h, z)
    prevL, prevB = L, B
    for _ in range(max_iter):
        resid = z - L * np.exp(-B * h)
        T = W @ resid                      # Eq. 4 / Eq. 7
        L, B = fit_stratified(h, z - T)    # Eq. 6
        if abs(L - prevL) < tol * (abs(prevL) + tol) and abs(B - prevB) < tol:
            break
        prevL, prevB = L, B

    residual = z - L * np.exp(-B * h)
    return float(L), float(B), residual


# --------------------------------------------------------------------------- #
# Neighbour search
# --------------------------------------------------------------------------- #
def _build_index(x, y):
    pts = np.column_stack([x, y])
    if _HAVE_SCIPY:
        return _KDTree(pts)
    return pts  # plain array; brute-force fallback used below


def _query_radius(index, qx, qy, radius):
    """Return a list of neighbour index arrays, one per query point."""
    if _HAVE_SCIPY:
        q = np.column_stack([qx, qy])
        return index.query_ball_point(q, r=radius)
    # Brute-force fallback (fine for modest obs counts).
    out = []
    r2 = radius * radius
    for i in range(qx.size):
        d2 = (index[:, 0] - qx[i]) ** 2 + (index[:, 1] - qy[i]) ** 2
        out.append(np.nonzero(d2 <= r2)[0].tolist())
    return out


# --------------------------------------------------------------------------- #
# Coverage mask (restrict output to the reliable area around the network)
# --------------------------------------------------------------------------- #
def coverage_mask(sx, sy, gx, gy, buffer_m=0.0):
    """Boolean grid mask of knots inside the station coverage.

    A knot is kept if it lies inside the convex hull of the stations OR within
    ``buffer_m`` of the nearest station. This removes the unreliable
    extrapolation/border zone beyond the network edge.
    """
    pts = np.column_stack([np.asarray(sx, float), np.asarray(sy, float)])
    q = np.column_stack([np.asarray(gx, float), np.asarray(gy, float)])
    keep = np.zeros(q.shape[0], dtype=bool)

    if pts.shape[0] >= 3:
        try:
            from scipy.spatial import Delaunay
            keep = Delaunay(pts).find_simplex(q) >= 0
        except Exception:
            keep = np.zeros(q.shape[0], dtype=bool)

    if buffer_m and buffer_m > 0:
        if _HAVE_SCIPY:
            d, _ = _KDTree(pts).query(q)
        else:
            d = np.array([
                np.sqrt(np.min((pts[:, 0] - q[i, 0]) ** 2
                               + (pts[:, 1] - q[i, 1]) ** 2))
                for i in range(q.shape[0])
            ])
        keep = keep | (d <= buffer_m)

    if not keep.any():
        # Safety: never blank the whole map (e.g. SciPy missing, buffer 0).
        keep = np.ones(q.shape[0], dtype=bool)
    return keep


# --------------------------------------------------------------------------- #
# AITD (Sec. 2.2) - the recommended, fast, artefact-free method
# --------------------------------------------------------------------------- #
def run_aitd(obs, dem, radius, max_iter=10, tol=1e-6,
             chunk=50000, progress=None, fill_nearest=True):
    """Augmented Iterative Tropospheric Decomposition.

    Parameters
    ----------
    obs : dict with metric arrays ``x``, ``y``, ``ztd`` and *scaled* ``h``.
    dem : dict with metric arrays ``x``, ``y`` and *scaled* ``h`` (flattened
          DEM knots, all of equal length).
    radius : decorrelation radius in metres.
    progress : optional callable ``progress(fraction, message)``.

    Returns ``dict`` with flat arrays ``ztd``, ``stratified``, ``turbulent``.
    """
    ox, oy, oh, oz = (np.asarray(obs[k], float) for k in ("x", "y", "h", "ztd"))
    gx, gy, gh = (np.asarray(dem[k], float) for k in ("x", "y", "h"))
    alpha = midw_alpha(radius)
    n_obs = oz.size

    def _report(frac, msg):
        if progress:
            progress(frac, msg)

    # ----- Step 1: decompose once per observation point (Sec. 2.2) ------ #
    _report(0.0, "AITD step 1/2: decomposing observations")
    L_o = np.empty(n_obs)
    B_o = np.empty(n_obs)
    res_o = np.empty(n_obs)  # turbulent value at each observation point
    index = _build_index(ox, oy)
    nbrs = _query_radius(index, ox, oy, radius)
    for i in range(n_obs):
        idx = np.asarray(nbrs[i], dtype=int)
        if idx.size == 0:
            idx = np.array([i])
        L, B, resid = decompose_window(
            oh[idx], oz[idx], ox[idx], oy[idx], "midw", alpha, max_iter, tol
        )
        L_o[i], B_o[i] = L, B
        # residual at the window centre i (its position within idx)
        centre = np.nonzero(idx == i)[0]
        res_o[i] = resid[centre[0]] if centre.size else float(oz[i] - L * np.exp(-B * oh[i]))
        if progress and (i % 25 == 0):
            _report(0.05 + 0.35 * (i / max(n_obs, 1)),
                    "AITD step 1/2: decomposing observations")

    # ----- Step 2: interpolate L, B and residuals onto the DEM ---------- #
    _report(0.4, "AITD step 2/2: predicting on the DEM grid")
    n_grid = gx.size
    ztd = np.full(n_grid, np.nan)
    strat = np.full(n_grid, np.nan)
    turb = np.full(n_grid, np.nan)
    r2 = radius * radius

    for start in range(0, n_grid, chunk):
        end = min(start + chunk, n_grid)
        cx = gx[start:end][:, None]
        cy = gy[start:end][:, None]
        ch = gh[start:end]
        d2 = (cx - ox[None, :]) ** 2 + (cy - oy[None, :]) ** 2  # (m, n_obs)
        within = d2 <= r2
        w = _pair_weights(d2, "midw", alpha) * within
        wsum = w.sum(axis=1)
        good = wsum > 0
        wn = np.zeros_like(w)
        wn[good] = w[good] / wsum[good, None]
        Li = wn @ L_o
        Bi = wn @ B_o
        Ti = wn @ res_o
        s = Li * np.exp(-Bi * ch)
        strat[start:end] = np.where(good, s, np.nan)
        turb[start:end] = np.where(good, Ti, np.nan)
        ztd[start:end] = np.where(good, s + Ti, np.nan)
        _report(0.4 + 0.6 * (end / max(n_grid, 1)),
                "AITD step 2/2: predicting on the DEM grid")

    if fill_nearest:
        _fill_holes_nearest(ztd, strat, turb, gx, gy, ox, oy,
                            L_o, B_o, res_o, gh)
    _report(1.0, "AITD complete")
    return {"ztd": ztd, "stratified": strat, "turbulent": turb}


# --------------------------------------------------------------------------- #
# ITD (Sec. 2.1) - the original method, for comparison
# --------------------------------------------------------------------------- #
def run_itd(obs, dem, radius, max_iter=10, tol=1e-6,
            chunk=20000, progress=None, fill_nearest=True):
    """Original Iterative Tropospheric Decomposition.

    The decomposition is, in principle, repeated for every DEM knot. As proven
    in the paper, the iterative decomposition depends *only* on the set of
    reference points inside the window, not on the knot position. We therefore
    group knots that share the same neighbour set and decompose once per unique
    group - an exact optimisation that keeps ITD tractable on large grids.
    """
    ox, oy, oh, oz = (np.asarray(obs[k], float) for k in ("x", "y", "h", "ztd"))
    gx, gy, gh = (np.asarray(dem[k], float) for k in ("x", "y", "h"))
    n_obs = oz.size
    n_grid = gx.size
    r2 = radius * radius

    def _report(frac, msg):
        if progress:
            progress(frac, msg)

    ztd = np.full(n_grid, np.nan)
    strat = np.full(n_grid, np.nan)
    turb = np.full(n_grid, np.nan)

    # Membership of each knot expressed as an integer bitmask over observations
    # (fast exact grouping when n_obs <= 63; otherwise frozenset fallback).
    use_bitmask = n_obs <= 63
    cache = {}  # window signature -> (L, B, residual_full[n_obs])

    _report(0.0, "ITD: predicting on the DEM grid")
    for start in range(0, n_grid, chunk):
        end = min(start + chunk, n_grid)
        cx = gx[start:end][:, None]
        cy = gy[start:end][:, None]
        ch = gh[start:end]
        d2 = (cx - ox[None, :]) ** 2 + (cy - oy[None, :]) ** 2
        within = d2 <= r2

        if use_bitmask:
            powers = (1 << np.arange(n_obs)).astype(np.int64)
            sig = within.astype(np.int64) @ powers
        else:
            sig = None

        for k in range(end - start):
            members = np.nonzero(within[k])[0]
            if members.size == 0:
                continue
            key = int(sig[k]) if use_bitmask else frozenset(members.tolist())
            if key not in cache:
                L, B, resid = decompose_window(
                    oh[members], oz[members], ox[members], oy[members],
                    "idw", 0.0, max_iter, tol,
                )
                full_res = np.full(n_obs, np.nan)
                full_res[members] = resid
                cache[key] = (L, B, full_res)
            L, B, full_res = cache[key]

            gi = start + k
            s = L * np.exp(-B * ch[k])
            # Turbulent prediction at the knot via classic IDW (Eq. 5).
            dd2 = d2[k, members]
            w = _pair_weights(dd2[None, :], "idw", 0.0)[0]
            wsum = w.sum()
            if wsum > 0:
                t = float((w / wsum) @ full_res[members])
            else:
                t = 0.0
            strat[gi] = s
            turb[gi] = t
            ztd[gi] = s + t
        _report(end / max(n_grid, 1), "ITD: predicting on the DEM grid")

    if fill_nearest:
        # For ITD the per-knot L/B vary; fall back to nearest observation's
        # decomposition only to fill empty knots outside every window.
        _fill_holes_itd(ztd, strat, turb, gx, gy, gh, ox, oy, oh, oz,
                        radius, max_iter, tol)
    _report(1.0, "ITD complete")
    return {"ztd": ztd, "stratified": strat, "turbulent": turb}


# --------------------------------------------------------------------------- #
# Hole filling (knots with no observation inside the decorrelation radius)
# --------------------------------------------------------------------------- #
def _fill_holes_nearest(ztd, strat, turb, gx, gy, ox, oy,
                        L_o, B_o, res_o, gh):
    holes = np.isnan(ztd)
    if not holes.any():
        return
    index = _build_index(ox, oy)
    hx, hy, hh = gx[holes], gy[holes], gh[holes]
    if _HAVE_SCIPY:
        _, nn = index.query(np.column_stack([hx, hy]), k=1)
    else:
        nn = np.array([
            int(np.argmin((ox - hx[i]) ** 2 + (oy - hy[i]) ** 2))
            for i in range(hx.size)
        ])
    s = L_o[nn] * np.exp(-B_o[nn] * hh)
    t = res_o[nn]
    strat[holes] = s
    turb[holes] = t
    ztd[holes] = s + t


def _fill_holes_itd(ztd, strat, turb, gx, gy, gh, ox, oy, oh, oz,
                    radius, max_iter, tol):
    holes = np.isnan(ztd)
    if not holes.any():
        return
    # Decompose the global set once to obtain a sensible L/B for far knots.
    alpha = midw_alpha(radius)
    L, B, resid = decompose_window(oh, oz, ox, oy, "idw", alpha, max_iter, tol)
    hx, hy, hh = gx[holes], gy[holes], gh[holes]
    s = L * np.exp(-B * hh)
    # nearest-observation turbulent value
    if _HAVE_SCIPY:
        _, nn = _KDTree(np.column_stack([ox, oy])).query(
            np.column_stack([hx, hy]), k=1)
    else:
        nn = np.array([
            int(np.argmin((ox - hx[i]) ** 2 + (oy - hy[i]) ** 2))
            for i in range(hx.size)
        ])
    t = resid[nn]
    strat[holes] = s
    turb[holes] = t
    ztd[holes] = s + t
