"""QgsTask that runs the full ZTD-AITD pipeline on a background thread.

Pipeline
--------
1. Read the ZTD table (lon, lat, height, ZTD).
2. Build a padded bounding box from the station coordinates.
3. Download a DEM (keyless AWS Terrain Tiles) or read a local one.
4. Reproject stations and DEM knots to a metric (UTM) CRS.
5. Scale heights to [0, 1] using a shared scale.
6. Run ITD and/or AITD.
7. Write the requested rasters (ZTD / stratified / turbulent / DEM) as GeoTIFFs.

Only :meth:`run` executes on the worker thread; it must not touch the GUI.
Layer loading happens back on the main thread in the plugin's callback.
"""

from __future__ import annotations

import os
import traceback

import numpy as np
from qgis.core import QgsTask, QgsMessageLog, Qgis

from ..core import decomposition as dc
from ..core import dem_downloader as dd
from ..core import io_utils as io
from ..core import projection as pj

MESSAGE_CATEGORY = "ZTD-AITD"


def _san(name):
    """Sanitise an epoch label into a filename-safe token."""
    s = "".join(c if c.isalnum() else "_" for c in str(name)).strip("_")
    return s.lower() or "epoch"


def map_filename(prefix, method, prod, epoch=None, multi=False):
    ep = f"_{_san(epoch)}" if (multi and epoch is not None) else ""
    return f"{prefix}_{method.lower()}_{prod}{ep}.tif"


def aps_filename(prefix, method, epoch):
    return f"{prefix}_{method.lower()}_ztd_aps_{_san(epoch)}.tif"


def planned_output_files(params):
    """Full paths the task will (over)write - used by the GUI to release file
    locks before overwriting layers already loaded in QGIS."""
    out_dir = params["output_dir"]
    prefix = params.get("prefix", "ztd")
    epochs = params.get("epoch_names") or [None]
    multi = len(epochs) > 1
    files = []
    if params.get("dem_mode") == "download":
        files.append(os.path.join(out_dir, f"{prefix}_dem.tif"))
    ref = params.get("reference_epoch")
    aps = (params.get("make_aps") and "ztd" in params["products"]
           and ref in epochs)
    for method in params["methods"]:
        for ep in epochs:
            for prod in params["products"]:
                files.append(os.path.join(
                    out_dir, map_filename(prefix, method, prod, ep, multi)))
            if aps and ep != ref:
                files.append(os.path.join(out_dir, aps_filename(prefix, method, ep)))
    if params.get("flag_extrapolation"):
        files.append(os.path.join(out_dir, f"{prefix}_extrapolation_mask.tif"))
    if params.get("add_stations_layer", True):
        files.append(os.path.join(out_dir, f"{prefix}_stations.gpkg"))
    return files


class AitdTask(QgsTask):
    """Run the decomposition pipeline. ``params`` is a plain dict (see GUI)."""

    def __init__(self, params):
        super().__init__("ZTD-AITD map generation", QgsTask.CanCancel)
        self.params = params
        self.outputs = {}      # {layer_title: file_path}
        self.exception = None
        self.message = ""
        self.summary = {}

    # ------------------------------------------------------------------ #
    def _progress(self, frac, msg):
        self.setProgress(max(0.0, min(100.0, frac * 100.0)))
        if msg:
            QgsMessageLog.logMessage(msg, MESSAGE_CATEGORY, Qgis.Info)

    def _log(self, msg, level=Qgis.Info):
        QgsMessageLog.logMessage(msg, MESSAGE_CATEGORY, level)

    # ------------------------------------------------------------------ #
    def run(self):
        try:
            return self._run()
        except Exception as exc:  # capture for the main thread
            self.exception = exc
            self.message = "".join(traceback.format_exception(exc))
            self._log(self.message, Qgis.Critical)
            return False

    def _run(self):
        p = self.params
        out_dir = p["output_dir"]
        os.makedirs(out_dir, exist_ok=True)
        prefix = p.get("prefix", "ztd")

        # --- 1. read ZTD epoch(s) -------------------------------------- #
        self._progress(0.02, "Reading ZTD table")
        mapping = p.get("mapping")
        ztd_columns = p.get("ztd_columns") or [mapping["ztd"]]
        table = io.read_ztd_epochs(
            p["ztd_path"], mapping, ztd_columns, delimiter=p.get("delimiter"))
        lon, lat, height = table["lon"], table["lat"], table["height"]
        epoch_names = p.get("epoch_names") or table["epoch_names"]
        epochs = table["epochs"]
        multi = len(epoch_names) > 1
        self._log("Read %d station(s), %d epoch(s): %s"
                  % (lon.size, len(epoch_names), ", ".join(epoch_names)))
        if self.isCanceled():
            return False

        # --- 2. bounding box ------------------------------------------- #
        bbox = io.bounding_box(lon, lat, margin_deg=p.get("margin_deg", 0.1))
        self._log("Bounding box (deg): %.3f, %.3f, %.3f, %.3f" % bbox)

        # --- 3. DEM ---------------------------------------------------- #
        if p["dem_mode"] == "download":
            dem_path = os.path.join(out_dir, f"{prefix}_dem.tif")
            self._progress(0.05, "Downloading DEM")
            dd.download_dem(
                bbox, dem_path,
                target_resolution_m=p.get("dem_resolution_m", 90.0),
                zoom=p.get("dem_zoom"),
                progress=lambda f, m: self._progress(0.05 + 0.25 * f, m),
                max_tiles=p.get("max_tiles", 600),
            )
        else:
            dem_path = p["dem_path"]
            self._log(f"Using local DEM: {dem_path}")
        if self.isCanceled():
            return False

        dem_z, gt, wkt, nodata = io.read_dem(dem_path)
        dlon, dlat, delev = io.dem_knot_coordinates(dem_z, gt)
        self._log("DEM grid: %d x %d knots." % (dem_z.shape[1], dem_z.shape[0]))

        finite = np.isfinite(delev)
        dem_lo = float(np.nanmin(delev)) if finite.any() else 0.0
        dem_hi = float(np.nanmax(delev)) if finite.any() else 0.0
        self._log("DEM elevation range: %.0f .. %.0f m" % (dem_lo, dem_hi))

        # Terrain tiles include sea/lake bathymetry (negative heights). The
        # stratified model S(h)=L*exp(-B*h) is only meaningful above ground, so
        # by default we clamp non-positive elevations to sea level (0 m); this
        # removes the unrealistic ZTD inflation seen below sea level.
        sea_floor = p.get("sea_floor", 0.0)
        if sea_floor is not None:
            n_below = int(np.nansum(delev < sea_floor))
            if n_below:
                self._log("Clamping %d DEM knot(s) below %.0f m to sea level."
                          % (n_below, sea_floor))
            delev = np.where(finite, np.clip(delev, sea_floor, None), delev)

        # Warn when the DEM spans elevations beyond the GNSS station range:
        # predictions there are extrapolations of the stratified model.
        s_hi = float(np.nanmax(height))
        if dem_hi > s_hi + 50.0:
            self._log(
                "NOTE: DEM reaches %.0f m but the highest station is %.0f m. "
                "ZTD above ~%.0f m is extrapolated and may read low at peaks "
                "(expected behaviour; see README)." % (dem_hi, s_hi, s_hi),
                Qgis.Warning)

        # --- 4. project to metric CRS ---------------------------------- #
        self._progress(0.32, "Reprojecting to metric CRS")
        sx, sy, epsg = pj.to_metric(lon, lat)
        dx, dy, _ = pj.to_metric(dlon, dlat, dst_epsg=epsg)
        self._log(f"Working metric CRS: EPSG:{epsg}")

        # --- 5. scale heights (shared scale) --------------------------- #
        all_h = np.concatenate([height, delev[finite]])
        _, hmin, hmax = dc.scale_heights(all_h)
        sh, _, _ = dc.scale_heights(height, hmin, hmax)
        dh, _, _ = dc.scale_heights(delev, hmin, hmax)
        # DEM knots with no elevation get height 0 -> masked out at the end
        dh = np.where(np.isfinite(dh), dh, 0.0)

        dem = {"x": dx, "y": dy, "h": dh}

        radius_m = p["radius_km"] * 1000.0
        max_iter = p.get("max_iter", 10)
        fill = p.get("fill_nearest", True)
        methods = p["methods"]      # subset of ["AITD", "ITD"]
        products = p["products"]    # subset of ztd/stratified/turbulent
        mask_nodata = ~finite
        clip = bool(p.get("clip_coverage", False))
        coverage_buffer_m = float(p.get("coverage_buffer_km", 5.0)) * 1000.0

        # APS: difference each epoch's ZTD map against a reference epoch's.
        reference = p.get("reference_epoch")
        make_aps = (bool(p.get("make_aps")) and "ztd" in products
                    and reference in epoch_names)
        order = epoch_names
        if make_aps:
            # process the reference epoch first so it is available to subtract
            order = [reference] + [e for e in epoch_names if e != reference]
            self._log(f"APS maps enabled; reference epoch: {reference}")
        ref_ztd = {}  # method -> reference-epoch ZTD map (flat)

        # --- 6. run each epoch x method -------------------------------- #
        total = max(len(order) * len(methods), 1)
        step = 0
        for name in order:
            zt = epochs[name]
            valid = np.isfinite(zt)
            nval = int(valid.sum())
            if nval < 3:
                self._log("Epoch '%s': only %d valid station(s); skipped."
                          % (name, nval), Qgis.Warning)
                step += len(methods)
                continue
            obs = {"x": sx[valid], "y": sy[valid],
                   "h": sh[valid], "ztd": zt[valid]}
            tag = f" [{name}]" if multi else ""
            self._log("Epoch '%s': %d valid station(s)." % (name, nval))

            # coverage clip (convex hull of this epoch's stations + buffer)
            outside = None
            if clip:
                keep = dc.coverage_mask(sx[valid], sy[valid], dx, dy,
                                        coverage_buffer_m)
                outside = ~keep
                self._log("Coverage clip: keeping %d of %d knot(s) "
                          "(hull + %.0f km buffer)."
                          % (int(keep.sum()), keep.size, coverage_buffer_m / 1000.0))

            for method in methods:
                if self.isCanceled():
                    return False
                base = 0.35 + 0.6 * (step / total)
                span = 0.6 / total

                def prog(f, m, _b=base, _s=span, _meth=method, _t=tag):
                    self._progress(_b + _s * f, f"{_meth}{_t}: {m}")

                runner = dc.run_aitd if method == "AITD" else dc.run_itd
                result = runner(obs, dem, radius_m, max_iter=max_iter,
                                progress=prog, fill_nearest=fill)
                for arr in result.values():
                    arr[mask_nodata] = np.nan
                    if outside is not None:
                        arr[outside] = np.nan

                for prod in products:
                    arr2d = result[prod].reshape(dem_z.shape)
                    fpath = os.path.join(
                        out_dir, map_filename(prefix, method, prod, name, multi))
                    io.write_raster(fpath, arr2d, gt, wkt)
                    title = f"{method} {prod}{tag} ({prefix})"
                    self.outputs[title] = fpath
                    vals = arr2d[np.isfinite(arr2d)]
                    if vals.size:
                        self.summary[title] = (float(np.mean(vals)),
                                                float(np.std(vals)))

                if make_aps:
                    if name == reference:
                        ref_ztd[method] = result["ztd"].copy()
                    elif method in ref_ztd:
                        aps2d = (result["ztd"] - ref_ztd[method]).reshape(dem_z.shape)
                        fpath = os.path.join(
                            out_dir, aps_filename(prefix, method, name))
                        io.write_raster(fpath, aps2d, gt, wkt)
                        title = f"{method} ZTD APS [{name}-{reference}] ({prefix})"
                        self.outputs[title] = fpath
                        vals = aps2d[np.isfinite(aps2d)]
                        if vals.size:
                            self.summary[title] = (float(np.mean(vals)),
                                                    float(np.std(vals)))
                step += 1

        # --- optional extrapolation mask ------------------------------- #
        # Pixels whose elevation exceeds the highest GNSS station fall outside
        # the calibration range of the stratified model: predictions there are
        # extrapolations and should be treated with caution. The mask carries
        # data (value 1) only over those knots; everywhere else it is NoData,
        # so it overlays cleanly on any ZTD map.
        if p.get("flag_extrapolation", False):
            s_hi = float(np.nanmax(height))
            extr = np.where(finite & (delev > s_hi), 1.0, np.nan)
            extr2d = extr.reshape(dem_z.shape)
            fpath = os.path.join(out_dir, f"{prefix}_extrapolation_mask.tif")
            io.write_raster(fpath, extr2d, gt, wkt)
            n_extr = int(np.nansum(extr))
            self.outputs[f"Extrapolation mask >{s_hi:.0f} m ({prefix})"] = fpath
            self._log("Extrapolation mask: %d knot(s) above the highest "
                      "station (%.0f m) flagged." % (n_extr, s_hi))

        # --- GNSS station points layer -------------------------------- #
        if p.get("add_stations_layer", True):
            spath = os.path.join(out_dir, f"{prefix}_stations.gpkg")
            fields = {"height": height}
            for name in epoch_names:
                fields[f"ztd_{_san(name)}"] = epochs[name]
            io.write_points(spath, lon, lat, fields, layer_name="stations")
            self.outputs[f"GNSS stations ({prefix})"] = spath
            self._log("Wrote %d GNSS station point(s)." % lon.size)

        # always expose the DEM
        if p.get("add_dem_layer", True):
            self.outputs[f"DEM ({prefix})"] = dem_path

        self._progress(1.0, "Done")
        return True

    # ------------------------------------------------------------------ #
    def finished(self, ok):
        # The plugin connects to taskCompleted/taskTerminated; this hook just
        # logs so the worker thread itself never touches the GUI.
        if ok:
            self._log(f"Task finished, {len(self.outputs)} layer(s) ready.")
        elif self.exception is not None:
            self._log("Task failed.", Qgis.Critical)
        else:
            self._log("Task cancelled.", Qgis.Warning)
