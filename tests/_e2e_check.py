"""End-to-end run mirroring AitdTask, using the bundled sample CSV."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from core import decomposition as dc, dem_downloader as dd, io_utils as io, projection as pj

here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv = os.path.join(here, "sample_data", "ztd_sample.csv")
out_dir = os.path.join(here, "sample_data", "_out")
os.makedirs(out_dir, exist_ok=True)

data, headers, mapping = io.read_ztd_table(csv)
print("auto-mapping:", mapping, "stations:", data["ztd"].size)
lon, lat, h, ztd = data["lon"], data["lat"], data["height"], data["ztd"]

bbox = io.bounding_box(lon, lat, 0.1)
print("bbox", tuple(round(b, 3) for b in bbox))
dem_path = os.path.join(out_dir, "ztd_dem.tif")
dd.download_dem(bbox, dem_path, zoom=8, progress=lambda f, m: None)  # coarse=fast
z, gt, wkt, nd = io.read_dem(dem_path)
dlon, dlat, delev = io.dem_knot_coordinates(z, gt)
print("DEM", z.shape, "knots", dlon.size)

sx, sy, epsg = pj.to_metric(lon, lat)
dx, dy, _ = pj.to_metric(dlon, dlat, dst_epsg=epsg)
finite = np.isfinite(delev)
_, hmin, hmax = dc.scale_heights(np.concatenate([h, delev[finite]]))
sh, _, _ = dc.scale_heights(h, hmin, hmax)
dh, _, _ = dc.scale_heights(delev, hmin, hmax)
dh = np.where(np.isfinite(dh), dh, 0.0)
obs = {"x": sx, "y": sy, "h": sh, "ztd": ztd}
dem = {"x": dx, "y": dy, "h": dh}

for method, fn in (("aitd", dc.run_aitd), ("itd", dc.run_itd)):
    r = fn(obs, dem, 120000.0, max_iter=8)
    for k in r:
        r[k][~finite] = np.nan
    for prod in ("ztd", "stratified", "turbulent"):
        a = r[prod].reshape(z.shape)
        p = os.path.join(out_dir, f"ztd_{method}_{prod}.tif")
        io.write_raster(p, a, gt, wkt)
        v = a[np.isfinite(a)]
        print(f"{method} {prod:10s} mean={np.mean(v):.4f} std={np.std(v)*1000:6.1f}mm -> {os.path.basename(p)}")

print("\nE2E OK. Files in", out_dir)
print(sorted(os.listdir(out_dir)))
