"""Quick integration check (run with QGIS's python). Not part of the unit tests."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from core import dem_downloader as dd, io_utils as io, projection as pj, decomposition as dc

bbox = (9.0, 45.3, 9.4, 45.6)
out = dd.download_dem(bbox, "_test_dem.tif", zoom=9, progress=lambda f, m: None)
z, gt, wkt, nd = io.read_dem(out)
print("DEM shape", z.shape, "elev min/max %.0f/%.0f" % (np.nanmin(z), np.nanmax(z)))
lon, lat, elev = io.dem_knot_coordinates(z, gt)
xm, ym, epsg = pj.to_metric(lon, lat)
print("utm epsg", epsg)

rng = np.random.default_rng(0)
slon = rng.uniform(9.05, 9.35, 8); slat = rng.uniform(45.35, 45.55, 8); sh = rng.uniform(100, 800, 8)
sztd = 2.3 * np.exp(-1.2 * ((sh - 100) / 700)) + 0.01 * np.sin(slon * 50)
sxm, sym, _ = pj.to_metric(slon, slat, dst_epsg=epsg)
hs, hmin, hmax = dc.scale_heights(np.concatenate([sh, elev[np.isfinite(elev)]]))
sh_s = hs[:8]
demh_s, _, _ = dc.scale_heights(elev, hmin, hmax)
obs = {"x": sxm, "y": sym, "h": sh_s, "ztd": sztd}
dem = {"x": xm, "y": ym, "h": demh_s}
r = dc.run_aitd(obs, dem, 80000.0)
print("AITD ztd finite %d  min/max %.4f/%.4f" %
      (np.isfinite(r["ztd"]).sum(), np.nanmin(r["ztd"]), np.nanmax(r["ztd"])))
io.write_raster("_test_ztd.tif", r["ztd"].reshape(z.shape), gt, wkt)
print("wrote _test_ztd.tif OK")
