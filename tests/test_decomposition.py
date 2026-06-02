"""Synthetic-data validation of the ITD / AITD core.

Run from the plugin root with:  python -m tests.test_decomposition
(Requires only NumPy; SciPy is used automatically if present.)

The test builds a known truth field
    ZTD(x, h) = L_true * exp(-B_true * h) + turbulent(x)
samples it at scattered "station" locations, predicts it back onto a grid with
both ITD and AITD, and checks that the leave-one-out prediction error is small.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import decomposition as dc  # noqa: E402


def _truth(x, y, h_scaled, L=2.4, B=1.3):
    """A smooth stratified field plus a low-frequency turbulent field."""
    strat = L * np.exp(-B * h_scaled)
    turb = 0.03 * np.sin(x / 40000.0) + 0.02 * np.cos(y / 60000.0)
    return strat + turb, strat, turb


def main():
    rng = np.random.default_rng(42)

    # 60 scattered stations over a ~200 x 200 km metric area.
    n = 60
    sx = rng.uniform(0, 200000, n)
    sy = rng.uniform(0, 200000, n)
    sh = rng.uniform(0, 2000, n)  # heights 0..2000 m
    h_scaled, hmin, hmax = dc.scale_heights(sh)
    sztd, _, _ = _truth(sx, sy, h_scaled)

    obs = {"x": sx, "y": sy, "h": h_scaled, "ztd": sztd}

    # Regular DEM grid (60 x 60) over the same extent.
    gx, gy = np.meshgrid(np.linspace(0, 200000, 60), np.linspace(0, 200000, 60))
    gh = rng.uniform(0, 2000, gx.size)
    gh_scaled, _, _ = dc.scale_heights(gh, hmin, hmax)
    dem = {"x": gx.ravel(), "y": gy.ravel(), "h": gh_scaled}
    truth_grid, strat_grid, turb_grid = _truth(dem["x"], dem["y"], dem["h"])

    radius = 120000.0  # 120 km decorrelation window

    for name, fn in (("ITD", dc.run_itd), ("AITD", dc.run_aitd)):
        out = fn(obs, dem, radius, max_iter=8)
        err = out["ztd"] - truth_grid
        rms = float(np.sqrt(np.nanmean(err ** 2)))
        # stratified/turbulent split should roughly recover the truth pieces
        srms = float(np.sqrt(np.nanmean((out["stratified"] - strat_grid) ** 2)))
        print(f"{name:5s}  ZTD grid RMS = {rms*1000:6.2f} mm   "
              f"stratified RMS = {srms*1000:6.2f} mm   "
              f"finite knots = {np.isfinite(out['ztd']).sum()}/{dem['x'].size}")
        assert np.isfinite(out["ztd"]).all(), f"{name}: NaNs left in ZTD map"
        assert rms < 0.010, f"{name}: ZTD RMS too high ({rms*1000:.1f} mm)"

    # Leave-one-out prediction error on the stations themselves (AITD).
    loo = np.empty(n)
    for i in range(n):
        keep = np.arange(n) != i
        sub = {k: obs[k][keep] for k in obs}
        one = {"x": sx[i:i+1], "y": sy[i:i+1], "h": h_scaled[i:i+1]}
        pred = dc.run_aitd(sub, one, radius, max_iter=8, fill_nearest=True)
        loo[i] = pred["ztd"][0] - sztd[i]
    loo_rms = float(np.sqrt(np.nanmean(loo ** 2)))
    print(f"AITD   leave-one-out station RMS = {loo_rms*1000:6.2f} mm")
    assert loo_rms < 0.015, f"LOO RMS too high ({loo_rms*1000:.1f} mm)"

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
