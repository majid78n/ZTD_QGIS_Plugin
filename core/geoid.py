"""Geoid height conversion: ellipsoidal (WGS-84) -> orthometric (mean sea level).

GNSS station heights are usually *ellipsoidal* (height above the WGS-84
ellipsoid), whereas DEM elevations are *orthometric* (height above the geoid /
mean sea level). The two differ by the geoid undulation ``N``::

    H_msl = h_ellipsoidal - N

The stratified delay model in this plugin assumes station heights and DEM
heights share the same vertical datum, so ellipsoidal station heights must be
converted to mean sea level before the decomposition.

We use GDAL/OSR (always present in QGIS) with a global geoid model, so this
module stays free of any PyQGIS dependency and unit-testable.
"""

from __future__ import annotations

import numpy as np

try:
    from osgeo import osr
    osr.UseExceptions()
except Exception:  # pragma: no cover
    osr = None


# Compound CRS (WGS-84 horizontal + orthometric height) for each geoid model.
# Source is always EPSG:4979 (WGS-84 geographic 3D, ellipsoidal height).
_GEOID_CRS = {
    "EGM2008": 9518,   # WGS 84 + EGM2008 height
    "EGM96": 9707,     # WGS 84 + EGM96 height
}

_SRC_EPSG = 4979


def available_models():
    """Return the geoid model names this module can convert to."""
    return list(_GEOID_CRS.keys())


def _srs(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(int(epsg))
    # Force lon/lat (x, y) axis order regardless of GDAL version.
    try:
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except Exception:
        pass
    return s


def ellipsoidal_to_orthometric(lon, lat, height, model="EGM2008"):
    """Convert ellipsoidal heights to orthometric (mean-sea-level) heights.

    Parameters
    ----------
    lon, lat : array-likes of longitude / latitude in degrees (WGS-84).
    height : array-like of ellipsoidal heights in metres (same shape).
    model : geoid model name (see :func:`available_models`).

    Returns a float array the same shape as ``height``. Raises ``RuntimeError``
    when GDAL/OSR or the required geoid grid is unavailable.
    """
    if osr is None:
        raise RuntimeError("GDAL/OSR is not available; cannot convert heights.")
    if model not in _GEOID_CRS:
        raise ValueError(
            f"Unknown geoid model '{model}'. Choose from {available_models()}.")
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    height = np.asarray(height, float)

    ct = osr.CoordinateTransformation(_srs(_SRC_EPSG), _srs(_GEOID_CRS[model]))
    pts = np.column_stack([lon.ravel(), lat.ravel(), height.ravel()]).tolist()
    res = np.asarray(ct.TransformPoints(pts), float)
    ortho = res[:, 2].reshape(height.shape)

    if not np.isfinite(ortho).any():
        raise RuntimeError(
            f"The {model} geoid grid is not available to PROJ; heights could "
            "not be converted to mean sea level.")
    return ortho


def geoid_undulation(lon, lat, model="EGM2008"):
    """Return the geoid undulation ``N = h_ellipsoidal - H_msl`` in metres."""
    zeros = np.zeros(np.shape(lon), float)
    return zeros - ellipsoidal_to_orthometric(lon, lat, zeros, model)
