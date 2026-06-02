"""Input/output helpers: ZTD table reading, bounding boxes and raster writing.

Kept free of PyQGIS imports (uses only the standard library, NumPy and GDAL,
all of which ship with QGIS) so the logic stays unit-testable.
"""

from __future__ import annotations

import csv
import math
import os

import numpy as np

try:
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()
except Exception:  # pragma: no cover
    gdal = None
    ogr = None
    osr = None


# Column names we will try to auto-detect, in priority order.
_LON_KEYS = ["lon", "longitude", "long", "x", "east", "easting"]
_LAT_KEYS = ["lat", "latitude", "y", "north", "northing"]
_HEIGHT_KEYS = ["h", "height", "altitude", "alt", "elevation", "elev",
                "z", "ortho", "orthometric"]
_ZTD_KEYS = ["ztd", "ztd_m", "ztd[m]", "delay", "zenith", "value"]


def _sniff_delimiter(sample):
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t ").delimiter
    except Exception:
        # Fall back to the most common separators by frequency.
        counts = {d: sample.count(d) for d in [",", ";", "\t"]}
        return max(counts, key=counts.get) or ","


def _match_column(headers, keys):
    """Return the index of the first header matching any key (case-insensitive)."""
    low = [h.strip().lower() for h in headers]
    for key in keys:
        for i, h in enumerate(low):
            if h == key:
                return i
    # relaxed: substring match
    for key in keys:
        for i, h in enumerate(low):
            if key in h:
                return i
    return None


def detect_columns(headers):
    """Guess (lon, lat, height, ztd) column indices from a header row.

    Returns a dict ``{'lon': i, 'lat': i, 'height': i, 'ztd': i}`` with values
    possibly ``None`` when no confident match is found (the GUI lets the user
    fix these manually).
    """
    return {
        "lon": _match_column(headers, _LON_KEYS),
        "lat": _match_column(headers, _LAT_KEYS),
        "height": _match_column(headers, _HEIGHT_KEYS),
        "ztd": _match_column(headers, _ZTD_KEYS),
    }


def read_ztd_table(path, mapping=None, delimiter=None):
    """Read a tabular ZTD file into arrays.

    Parameters
    ----------
    path : path to a CSV/TSV/whitespace-delimited text file.
    mapping : optional dict with integer column indices for
              ``lon``, ``lat``, ``height``, ``ztd``. If omitted the columns are
              auto-detected from the header.
    delimiter : optional explicit delimiter; auto-sniffed when omitted.

    Returns ``(data, headers, mapping)`` where ``data`` is a dict of float
    arrays ``lon``, ``lat``, ``height``, ``ztd``.
    """
    with open(path, "r", newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = delimiter or _sniff_delimiter(sample)
        reader = csv.reader(fh, delimiter=delim)
        rows = [r for r in reader if any(c.strip() for c in r)]

    if not rows:
        raise ValueError("The ZTD file is empty.")

    # Detect whether the first row is a header (non-numeric first cell).
    def _is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    has_header = not all(_is_number(c) for c in rows[0])
    headers = rows[0] if has_header else [f"col{i}" for i in range(len(rows[0]))]
    body = rows[1:] if has_header else rows

    if mapping is None:
        mapping = detect_columns(headers)
    missing = [k for k, v in mapping.items() if v is None]
    if missing:
        raise ValueError(
            "Could not identify column(s): " + ", ".join(missing)
            + f". Detected headers: {headers}"
        )

    cols = {k: [] for k in ("lon", "lat", "height", "ztd")}
    for r in body:
        try:
            for k in cols:
                cols[k].append(float(r[mapping[k]]))
        except (ValueError, IndexError):
            continue  # skip malformed rows silently

    data = {k: np.asarray(v, dtype=float) for k, v in cols.items()}
    if data["ztd"].size == 0:
        raise ValueError("No valid numeric rows were read from the ZTD file.")
    return data, headers, mapping


def detect_epoch_columns(headers):
    """Return indices of columns that look like ZTD epochs.

    Matches headers such as ``ZTD``, ``ZTD1``, ``ZTD_2``, ``ztd 03`` etc., i.e.
    any header containing 'ztd' (case-insensitive). Used to pre-select epochs in
    the time-series UI.
    """
    out = []
    for i, h in enumerate(headers):
        if "ztd" in h.strip().lower():
            out.append(i)
    return out


def _read_rows(path, delimiter=None):
    """Return ``(headers, body_rows)`` from a delimited text file."""
    with open(path, "r", newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = delimiter or _sniff_delimiter(sample)
        reader = csv.reader(fh, delimiter=delim)
        rows = [r for r in reader if any(c.strip() for c in r)]
    if not rows:
        raise ValueError("The ZTD file is empty.")

    def _is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    has_header = not all(_is_number(c) for c in rows[0])
    headers = rows[0] if has_header else [f"col{i}" for i in range(len(rows[0]))]
    body = rows[1:] if has_header else rows
    return headers, body


def read_ztd_epochs(path, mapping, ztd_indices, delimiter=None):
    """Read station coordinates plus one or more ZTD epoch columns.

    Parameters
    ----------
    mapping : dict with integer indices for ``lon``, ``lat``, ``height``.
    ztd_indices : list of column indices, each treated as one epoch.

    Returns a dict::

        {
          "lon": (n,), "lat": (n,), "height": (n,),   # stations with valid coords
          "epoch_names": [header, ...],
          "epochs": {header: (n,) array, NaN where the value is missing},
        }

    Stations whose lon/lat/height cannot be parsed are dropped. Per-epoch ZTD
    values that are blank or non-numeric become NaN (the caller filters them out
    for that epoch).
    """
    headers, body = _read_rows(path, delimiter)
    for key in ("lon", "lat", "height"):
        if mapping.get(key) is None:
            raise ValueError(f"Missing column mapping for '{key}'.")
    if not ztd_indices:
        raise ValueError("No ZTD epoch column selected.")
    names = [headers[i] if 0 <= i < len(headers) else f"ZTD{i}" for i in ztd_indices]

    lon, lat, hgt = [], [], []
    epoch_vals = {n: [] for n in names}
    for r in body:
        try:
            lo = float(r[mapping["lon"]])
            la = float(r[mapping["lat"]])
            he = float(r[mapping["height"]])
        except (ValueError, IndexError):
            continue  # invalid coordinates -> drop station
        lon.append(lo)
        lat.append(la)
        hgt.append(he)
        for n, ci in zip(names, ztd_indices):
            try:
                epoch_vals[n].append(float(r[ci]))
            except (ValueError, IndexError):
                epoch_vals[n].append(np.nan)  # missing ZTD for this epoch

    if not lon:
        raise ValueError("No valid stations (lon/lat/height) read from the file.")
    return {
        "lon": np.asarray(lon, float),
        "lat": np.asarray(lat, float),
        "height": np.asarray(hgt, float),
        "epoch_names": names,
        "epochs": {n: np.asarray(v, float) for n, v in epoch_vals.items()},
    }


def bounding_box(lon, lat, margin_deg=0.1):
    """Return ``(min_lon, min_lat, max_lon, max_lat)`` padded by ``margin_deg``.

    The margin prevents stations at the edge from sitting exactly on the DEM
    border (the decorrelation windows need some surrounding terrain).
    """
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    return (
        float(lon.min()) - margin_deg,
        float(lat.min()) - margin_deg,
        float(lon.max()) + margin_deg,
        float(lat.max()) + margin_deg,
    )


# --------------------------------------------------------------------------- #
# Raster reading / writing (GDAL)
# --------------------------------------------------------------------------- #
def read_dem(path):
    """Read a single-band DEM into ``(z, gt, wkt, nodata)``.

    ``z`` is a 2D float array (north-up), ``gt`` the GDAL geotransform and
    ``wkt`` the CRS as WKT.
    """
    if gdal is None:
        raise RuntimeError("GDAL is not available.")
    ds = gdal.Open(path)
    if ds is None:
        raise RuntimeError(f"Could not open DEM: {path}")
    band = ds.GetRasterBand(1)
    z = band.ReadAsArray().astype(float)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        z[z == nodata] = np.nan
    gt = ds.GetGeoTransform()
    wkt = ds.GetProjection()
    ds = None
    return z, gt, wkt, nodata


def dem_knot_coordinates(z, gt):
    """Return flattened (lon, lat, value) arrays for every DEM pixel centre."""
    rows, cols = z.shape
    col_idx = np.arange(cols)
    row_idx = np.arange(rows)
    cc, rr = np.meshgrid(col_idx, row_idx)
    # GDAL geotransform: x = gt0 + col*gt1 + row*gt2 ; y = gt3 + col*gt4 + row*gt5
    x = gt[0] + (cc + 0.5) * gt[1] + (rr + 0.5) * gt[2]
    y = gt[3] + (cc + 0.5) * gt[4] + (rr + 0.5) * gt[5]
    return x.ravel(), y.ravel(), z.ravel()


def write_raster(path, array2d, gt, wkt, nodata=-9999.0):
    """Write a 2D float array to a GeoTIFF."""
    if gdal is None:
        raise RuntimeError("GDAL is not available.")
    rows, cols = array2d.shape
    # Remove any stale file first so GDAL does not fail trying to delete it.
    if os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, cols, rows, 1, gdal.GDT_Float32,
                    options=["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    out = np.where(np.isfinite(array2d), array2d, nodata).astype(np.float32)
    band.WriteArray(out)
    band.SetNoDataValue(nodata)
    band.FlushCache()
    ds = None
    return path


def _safe_field(name):
    """Make a GeoPackage-safe field name (alnum/underscore, not digit-leading)."""
    s = "".join(c if c.isalnum() else "_" for c in str(name)).strip("_")
    if not s:
        s = "f"
    if s[0].isdigit():
        s = "f_" + s
    return s[:60]


def write_points(path, lon, lat, fields=None, layer_name="stations", epsg=4326):
    """Write station points to a GeoPackage.

    ``fields`` is an optional dict ``{name: array}`` of per-point attributes
    (floats); non-finite values are stored as NULL. ``lon``/``lat`` are always
    added as attributes too.
    """
    if ogr is None:
        raise RuntimeError("GDAL/OGR is not available.")
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    fields = fields or {}

    drv = ogr.GetDriverByName("GPKG")
    if os.path.exists(path):
        try:
            drv.DeleteDataSource(path)
        except Exception:
            try:
                os.remove(path)
            except OSError:
                pass
    ds = drv.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    try:
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except Exception:
        pass
    layer = ds.CreateLayer(layer_name, srs, ogr.wkbPoint)

    field_map = {}  # original name -> safe name
    for base in (["lon", "lat"] + list(fields.keys())):
        safe = _safe_field(base)
        field_map[base] = safe
        layer.CreateField(ogr.FieldDefn(safe, ogr.OFTReal))

    defn = layer.GetLayerDefn()
    for i in range(lon.size):
        feat = ogr.Feature(defn)
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint_2D(float(lon[i]), float(lat[i]))
        feat.SetGeometry(pt)
        feat.SetField(field_map["lon"], float(lon[i]))
        feat.SetField(field_map["lat"], float(lat[i]))
        for base, arr in fields.items():
            v = arr[i]
            if np.isfinite(v):
                feat.SetField(field_map[base], float(v))
        layer.CreateFeature(feat)
        feat = None
    ds = None
    return path


def haversine_spacing_m(lat_deg):
    """Approximate metres-per-degree at a given latitude (for grid sizing)."""
    lat = math.radians(lat_deg)
    m_per_deg_lat = 111132.92 - 559.82 * math.cos(2 * lat) + 1.175 * math.cos(4 * lat)
    m_per_deg_lon = 111412.84 * math.cos(lat) - 93.5 * math.cos(3 * lat)
    return m_per_deg_lat, m_per_deg_lon
