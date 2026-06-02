"""Keyless DEM download from the AWS *Terrain Tiles* open dataset.

https://registry.opendata.aws/terrain-tiles/

Two endpoints are used:
  * ``/geotiff/{z}/{x}/{y}.tif``      - single-band 16-bit elevation (preferred)
  * ``/terrarium/{z}/{x}/{y}.png``    - RGB-encoded elevation (fallback)

Tiles covering the requested bounding box are fetched, mosaicked into a VRT and
warped to EPSG:4326 (lon/lat) clipped to the box, matching the GNSS station
coordinates. No API key or registration is required.

Only the standard library, NumPy and GDAL are used (all bundled with QGIS).
"""

from __future__ import annotations

import math
import os
import tempfile
import urllib.request

import numpy as np

try:
    from osgeo import gdal
    gdal.UseExceptions()
except Exception:  # pragma: no cover
    gdal = None

GEOTIFF_URL = "https://s3.amazonaws.com/elevation-tiles-prod/geotiff/{z}/{x}/{y}.tif"
TERRARIUM_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"

# GeoTIFF tiles are 512x512 in Web Mercator.
_TILE_PX = 512
_EARTH_CIRC = 40075016.686


def zoom_for_resolution(target_m, lat=0.0):
    """Pick the smallest zoom whose pixel size is <= ``target_m`` (approx).

    Web-Mercator ground resolution shrinks with cos(lat); we size for the box
    centre so the result is at least as fine as requested.
    """
    res_equator = _EARTH_CIRC * math.cos(math.radians(lat))
    z = math.log2(res_equator / (target_m * _TILE_PX))
    return int(max(0, min(13, math.ceil(z))))


def _lon2tilex(lon, z):
    return int((lon + 180.0) / 360.0 * (1 << z))


def _lat2tiley(lat, z):
    lat = max(min(lat, 85.05112878), -85.05112878)
    r = math.radians(lat)
    return int((1.0 - math.log(math.tan(r) + 1.0 / math.cos(r)) / math.pi)
               / 2.0 * (1 << z))


def tiles_for_bbox(bbox, z):
    """Yield ``(x, y)`` tile indices covering ``bbox=(minlon,minlat,maxlon,maxlat)``."""
    minlon, minlat, maxlon, maxlat = bbox
    x0 = _lon2tilex(minlon, z)
    x1 = _lon2tilex(maxlon, z)
    y0 = _lat2tiley(maxlat, z)   # note: y grows southward
    y1 = _lat2tiley(minlat, z)
    n = 1 << z
    for x in range(min(x0, x1), max(x0, x1) + 1):
        for y in range(min(y0, y1), max(y0, y1) + 1):
            yield x % n, max(0, min(y, n - 1))


def estimate_tile_count(bbox, z):
    return sum(1 for _ in tiles_for_bbox(bbox, z))


def _download(url, dest, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "QGIS-ZTD-AITD/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    with open(dest, "wb") as fh:
        fh.write(data)
    return dest


def _terrarium_to_geotiff(png_path, tif_path):
    """Decode a terrarium PNG (elevation = R*256 + G + B/256 - 32768)."""
    ds = gdal.Open(png_path)
    r = ds.GetRasterBand(1).ReadAsArray().astype(float)
    g = ds.GetRasterBand(2).ReadAsArray().astype(float)
    b = ds.GetRasterBand(3).ReadAsArray().astype(float)
    elev = (r * 256.0 + g + b / 256.0) - 32768.0
    drv = gdal.GetDriverByName("GTiff")
    out = drv.Create(tif_path, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    out.SetGeoTransform(ds.GetGeoTransform())
    out.SetProjection(ds.GetProjection())
    out.GetRasterBand(1).WriteArray(elev.astype(np.float32))
    out.FlushCache()
    out = None
    ds = None
    return tif_path


def download_dem(bbox, out_path, target_resolution_m=90.0, zoom=None,
                 progress=None, max_tiles=600):
    """Download and assemble a DEM for ``bbox`` into ``out_path`` (GeoTIFF).

    Parameters
    ----------
    bbox : (min_lon, min_lat, max_lon, max_lat) in degrees.
    out_path : destination GeoTIFF (EPSG:4326).
    target_resolution_m : desired ground sampling; used to pick the zoom level
        when ``zoom`` is not given.
    zoom : explicit XYZ zoom level (overrides ``target_resolution_m``).
    progress : optional ``progress(fraction, message)`` callback.
    max_tiles : safety cap; raises if the box would need more tiles.

    Returns ``out_path``.
    """
    if gdal is None:
        raise RuntimeError("GDAL is not available - cannot assemble the DEM.")
    minlon, minlat, maxlon, maxlat = bbox
    lat_c = 0.5 * (minlat + maxlat)
    if zoom is None:
        zoom = zoom_for_resolution(target_resolution_m, lat_c)

    tiles = list(tiles_for_bbox(bbox, zoom))
    if len(tiles) > max_tiles:
        raise RuntimeError(
            f"Requested area needs {len(tiles)} tiles at zoom {zoom} "
            f"(limit {max_tiles}). Use a coarser resolution or smaller area."
        )

    tmpdir = tempfile.mkdtemp(prefix="ztd_dem_")
    local = []
    n = len(tiles)
    for i, (x, y) in enumerate(tiles):
        tif = os.path.join(tmpdir, f"{zoom}_{x}_{y}.tif")
        try:
            _download(GEOTIFF_URL.format(z=zoom, x=x, y=y), tif)
        except Exception:
            # fall back to terrarium PNG and decode
            try:
                png = os.path.join(tmpdir, f"{zoom}_{x}_{y}.png")
                _download(TERRARIUM_URL.format(z=zoom, x=x, y=y), png)
                _terrarium_to_geotiff(png, tif)
            except Exception:
                continue  # missing tile (e.g. ocean); skip
        local.append(tif)
        if progress:
            progress(0.7 * (i + 1) / n, f"Downloading DEM tiles {i+1}/{n}")

    if not local:
        raise RuntimeError("No DEM tiles could be downloaded for this area.")

    if progress:
        progress(0.8, "Mosaicking DEM tiles")
    vrt = os.path.join(tmpdir, "mosaic.vrt")
    gdal.BuildVRT(vrt, local)

    if progress:
        progress(0.9, "Warping DEM to lon/lat and clipping")
    gdal.Warp(
        out_path, vrt,
        dstSRS="EPSG:4326",
        outputBounds=[minlon, minlat, maxlon, maxlat],
        resampleAlg="bilinear",
        format="GTiff",
        creationOptions=["COMPRESS=DEFLATE", "TILED=YES"],
    )
    if progress:
        progress(1.0, "DEM ready")
    return out_path
