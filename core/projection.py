"""Projection helpers - turn geographic coordinates into metres.

The decomposition works with Euclidean distances, so all horizontal
coordinates must be in a metric CRS. We auto-pick the UTM zone covering the
data and reproject with GDAL/OSR (always present in QGIS), keeping this module
free of any PyQGIS dependency.
"""

from __future__ import annotations

import numpy as np

try:
    from osgeo import osr
    osr.UseExceptions()
except Exception:  # pragma: no cover
    osr = None


def utm_epsg_for(lon, lat):
    """Return the EPSG code of the UTM zone containing (lon, lat)."""
    zone = int((lon + 180.0) / 6.0) + 1
    zone = min(max(zone, 1), 60)
    return (32600 if lat >= 0 else 32700) + zone


def _srs(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(int(epsg))
    # Force lon/lat (x, y) axis order regardless of GDAL version.
    try:
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except Exception:
        pass
    return s


def make_transformer(src_epsg, dst_epsg):
    """Return a callable ``f(x, y) -> (x2, y2)`` reprojecting arrays."""
    if osr is None:
        raise RuntimeError("GDAL/OSR is not available.")
    ct = osr.CoordinateTransformation(_srs(src_epsg), _srs(dst_epsg))

    def _f(x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        pts = np.column_stack([x, y, np.zeros_like(x)])
        res = np.array(ct.TransformPoints(pts.tolist()))
        return res[:, 0], res[:, 1]

    return _f


def to_metric(lon, lat, src_epsg=4326, dst_epsg=None):
    """Reproject lon/lat to a metric CRS.

    Returns ``(x_m, y_m, dst_epsg)``. When ``dst_epsg`` is None an appropriate
    UTM zone is chosen from the centroid of the points.
    """
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    if dst_epsg is None:
        dst_epsg = utm_epsg_for(float(np.mean(lon)), float(np.mean(lat)))
    f = make_transformer(src_epsg, dst_epsg)
    x, y = f(lon, lat)
    return x, y, dst_epsg
