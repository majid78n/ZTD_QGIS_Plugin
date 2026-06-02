"""Standalone control window for the ZTD-AITD plugin.

Implemented as an independent, resizable, minimizable top-level window (QDialog
with the ``Qt.Window`` flag) rather than a docked panel, so the user can move,
resize, maximize and minimize it freely. The whole UI is built in code (no
compiled .ui / .qrc).
"""

from __future__ import annotations

import os

from qgis.PyQt.QtCore import Qt, pyqtSignal
from qgis.PyQt.QtGui import QColor
from qgis.PyQt.QtWidgets import (
    QWidget, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QComboBox, QSpinBox,
    QDoubleSpinBox, QCheckBox, QRadioButton, QButtonGroup, QProgressBar,
    QPlainTextEdit, QFileDialog, QScrollArea, QListWidget, QListWidgetItem,
)
from qgis.core import (
    QgsApplication, QgsProject, QgsRasterLayer, QgsVectorLayer, Qgis,
    QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer,
    QgsStyle, QgsMarkerSymbol,
)

from ..core import io_utils as io
from ..workers.aitd_task import AitdTask, planned_output_files


class ZtdMainWindow(QDialog):
    """Independent window with the full ZTD-AITD workflow."""

    closed = pyqtSignal()

    def __init__(self, iface, parent=None):
        super().__init__(parent)
        self.iface = iface
        self.task = None
        self._mapping_combos = {}
        self.setWindowTitle("ZTD-AITD Map Generator")
        # Independent top-level window with min / max / close buttons.
        self.setWindowFlags(Qt.Window)
        self.resize(580, 820)
        self._build_ui()

    # ------------------------------------------------------------------ #
    def closeEvent(self, event):
        self.closed.emit()
        super().closeEvent(event)

    # ------------------------------------------------------------------ #
    # UI construction
    # ------------------------------------------------------------------ #
    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        container = QWidget()
        scroll.setWidget(container)
        outer.addWidget(scroll)

        root = QVBoxLayout(container)
        root.addWidget(self._build_ztd_group())
        root.addWidget(self._build_dem_group())
        root.addWidget(self._build_params_group())
        root.addWidget(self._build_output_group())
        root.addWidget(self._build_run_group())
        root.addStretch(1)

    # --- ZTD input ------------------------------------------------------ #
    def _build_ztd_group(self):
        g = QGroupBox("1. GNSS ZTD observations")
        lay = QVBoxLayout(g)

        row = QHBoxLayout()
        self.ztd_path = QLineEdit()
        self.ztd_path.setPlaceholderText("Path to a CSV / tabular ZTD file")
        btn = QPushButton("Browse…")
        btn.clicked.connect(self._pick_ztd)
        row.addWidget(self.ztd_path)
        row.addWidget(btn)
        lay.addLayout(row)

        drow = QHBoxLayout()
        drow.addWidget(QLabel("Delimiter:"))
        self.delim = QComboBox()
        self.delim.addItems(["auto", ",", ";", "tab", "space"])
        drow.addWidget(self.delim)
        reload_btn = QPushButton("Read columns")
        reload_btn.clicked.connect(self._load_columns)
        drow.addWidget(reload_btn)
        drow.addStretch(1)
        lay.addLayout(drow)

        grid = QGridLayout()
        labels = [("lon", "Longitude"), ("lat", "Latitude"),
                  ("height", "Height"), ("ztd", "ZTD")]
        for i, (key, text) in enumerate(labels):
            grid.addWidget(QLabel(text + ":"), i // 2, (i % 2) * 2)
            cb = QComboBox()
            cb.setEnabled(False)
            self._mapping_combos[key] = cb
            grid.addWidget(cb, i // 2, (i % 2) * 2 + 1)
        lay.addLayout(grid)

        # --- time-series (multi-epoch) controls ------------------------ #
        self.cb_timeseries = QCheckBox(
            "Time series: process multiple ZTD epoch columns")
        self.cb_timeseries.toggled.connect(self._toggle_timeseries)
        lay.addWidget(self.cb_timeseries)

        self.epoch_box = QWidget()
        ebl = QVBoxLayout(self.epoch_box)
        ebl.setContentsMargins(16, 0, 0, 0)
        ebl.addWidget(QLabel("ZTD epoch columns to process:"))
        self.epoch_list = QListWidget()
        self.epoch_list.setMaximumHeight(110)
        ebl.addWidget(self.epoch_list)
        arow = QHBoxLayout()
        self.cb_aps = QCheckBox("Also output APS-like difference maps")
        self.cb_aps.toggled.connect(self._toggle_aps)
        arow.addWidget(self.cb_aps)
        arow.addWidget(QLabel("reference epoch:"))
        self.ref_combo = QComboBox()
        self.ref_combo.setEnabled(False)
        arow.addWidget(self.ref_combo)
        arow.addStretch(1)
        ebl.addLayout(arow)
        self.cb_aps_symmetric = QCheckBox(
            "APS legend symmetric about zero (white = no change)")
        self.cb_aps_symmetric.setChecked(True)
        ebl.addWidget(self.cb_aps_symmetric)
        self.epoch_box.setVisible(False)
        lay.addWidget(self.epoch_box)
        return g

    def _toggle_timeseries(self, on):
        self.epoch_box.setVisible(on)
        # the single ZTD combo is only used when not in time-series mode
        self._mapping_combos["ztd"].setEnabled(not on and self.epoch_list.count() > 0)

    def _toggle_aps(self, on):
        self.ref_combo.setEnabled(on)

    def _populate_epochs(self, headers, epoch_indices):
        """Fill the epoch checklist; pre-check ZTD-like columns."""
        self.epoch_list.clear()
        preset = set(epoch_indices)
        for i, h in enumerate(headers):
            item = QListWidgetItem(h)
            item.setData(Qt.UserRole, i)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if i in preset else Qt.Unchecked)
            self.epoch_list.addItem(item)
        self.epoch_list.itemChanged.connect(self._refresh_ref_combo)
        self._refresh_ref_combo()

    def _checked_epochs(self):
        """Return (indices, names) of the ticked epoch columns, in order."""
        idx, names = [], []
        for row in range(self.epoch_list.count()):
            it = self.epoch_list.item(row)
            if it.checkState() == Qt.Checked:
                idx.append(int(it.data(Qt.UserRole)))
                names.append(it.text())
        return idx, names

    def _refresh_ref_combo(self, *args):
        _, names = self._checked_epochs()
        current = self.ref_combo.currentText()
        self.ref_combo.clear()
        self.ref_combo.addItems(names)
        if current in names:
            self.ref_combo.setCurrentIndex(names.index(current))

    # --- DEM ------------------------------------------------------------ #
    def _build_dem_group(self):
        g = QGroupBox("2. Digital Elevation Model")
        lay = QVBoxLayout(g)

        self.dem_group = QButtonGroup(self)
        self.rb_download = QRadioButton(
            "Download automatically (AWS Terrain Tiles, no key)")
        self.rb_loaded = QRadioButton("Use a DEM already loaded in QGIS")
        self.rb_local = QRadioButton("Use a local DEM file")
        self.rb_download.setChecked(True)
        for rb in (self.rb_download, self.rb_loaded, self.rb_local):
            self.dem_group.addButton(rb)

        lay.addWidget(self.rb_download)
        drow = QHBoxLayout()
        drow.addWidget(QLabel("Target resolution (m):"))
        self.dem_res = QSpinBox()
        self.dem_res.setRange(10, 1000)
        self.dem_res.setValue(90)
        self.dem_res.setSingleStep(10)
        drow.addWidget(self.dem_res)
        drow.addWidget(QLabel("Border margin (deg):"))
        self.margin = QDoubleSpinBox()
        self.margin.setRange(0.0, 2.0)
        self.margin.setValue(0.1)
        self.margin.setSingleStep(0.05)
        drow.addWidget(self.margin)
        drow.addStretch(1)
        lay.addLayout(drow)

        # loaded raster option
        lay.addWidget(self.rb_loaded)
        lrow = QHBoxLayout()
        self.loaded_combo = QComboBox()
        self.loaded_combo.setEnabled(False)
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self._refresh_loaded_rasters)
        lrow.addWidget(self.loaded_combo)
        lrow.addWidget(refresh)
        lay.addLayout(lrow)

        # local file option
        lay.addWidget(self.rb_local)
        frow = QHBoxLayout()
        self.dem_path = QLineEdit()
        self.dem_path.setPlaceholderText("Path to a local DEM raster")
        self.dem_path.setEnabled(False)
        self.dem_btn = QPushButton("Browse…")
        self.dem_btn.setEnabled(False)
        self.dem_btn.clicked.connect(self._pick_dem)
        frow.addWidget(self.dem_path)
        frow.addWidget(self.dem_btn)
        lay.addLayout(frow)

        self.rb_download.toggled.connect(self._update_dem_mode)
        self.rb_loaded.toggled.connect(self._update_dem_mode)
        self.rb_local.toggled.connect(self._update_dem_mode)
        self._refresh_loaded_rasters()
        return g

    def _update_dem_mode(self):
        download = self.rb_download.isChecked()
        loaded = self.rb_loaded.isChecked()
        local = self.rb_local.isChecked()
        self.dem_res.setEnabled(download)
        self.margin.setEnabled(download)
        self.loaded_combo.setEnabled(loaded)
        self.dem_path.setEnabled(local)
        self.dem_btn.setEnabled(local)
        if loaded:
            self._refresh_loaded_rasters()

    def _refresh_loaded_rasters(self):
        current = self.loaded_combo.currentData()
        self.loaded_combo.clear()
        self.loaded_combo.addItem("(choose a loaded raster layer)", None)
        for lid, layer in QgsProject.instance().mapLayers().items():
            if isinstance(layer, QgsRasterLayer) and layer.isValid():
                self.loaded_combo.addItem(layer.name(), lid)
        if current is not None:
            idx = self.loaded_combo.findData(current)
            if idx >= 0:
                self.loaded_combo.setCurrentIndex(idx)

    # --- parameters ----------------------------------------------------- #
    def _build_params_group(self):
        g = QGroupBox("3. Interpolation parameters")
        form = QFormLayout(g)

        self.radius = QDoubleSpinBox()
        self.radius.setRange(1.0, 1000.0)
        self.radius.setValue(100.0)
        self.radius.setSuffix(" km")
        form.addRow("Decorrelation radius:", self.radius)

        self.max_iter = QSpinBox()
        self.max_iter.setRange(1, 50)
        self.max_iter.setValue(10)
        form.addRow("Max decomposition iterations:", self.max_iter)

        mrow = QHBoxLayout()
        self.cb_aitd = QCheckBox("AITD (augmented)")
        self.cb_aitd.setChecked(True)
        self.cb_itd = QCheckBox("ITD (original)")
        mrow.addWidget(self.cb_aitd)
        mrow.addWidget(self.cb_itd)
        mrow.addStretch(1)
        form.addRow("Method:", self._wrap(mrow))

        prow = QHBoxLayout()
        self.cb_ztd = QCheckBox("ZTD")
        self.cb_strat = QCheckBox("Stratified")
        self.cb_turb = QCheckBox("Turbulent")
        for cb in (self.cb_ztd, self.cb_strat, self.cb_turb):
            cb.setChecked(True)
            prow.addWidget(cb)
        prow.addStretch(1)
        form.addRow("Output maps:", self._wrap(prow))

        self.cb_fill = QCheckBox("Fill knots outside every window (nearest)")
        self.cb_fill.setChecked(True)
        form.addRow("", self.cb_fill)

        self.cb_sea = QCheckBox("Clamp DEM below sea level to 0 m (recommended)")
        self.cb_sea.setChecked(True)
        form.addRow("", self.cb_sea)

        self.cb_flag = QCheckBox(
            "Flag pixels above the highest station as extrapolated")
        self.cb_flag.setChecked(False)
        form.addRow("", self.cb_flag)

        crow = QHBoxLayout()
        self.cb_clip = QCheckBox("Clip output to station coverage (convex hull)")
        self.cb_clip.setChecked(False)
        self.cb_clip.toggled.connect(lambda on: self.coverage_buffer.setEnabled(on))
        crow.addWidget(self.cb_clip)
        crow.addWidget(QLabel("buffer:"))
        self.coverage_buffer = QDoubleSpinBox()
        self.coverage_buffer.setRange(0.0, 200.0)
        self.coverage_buffer.setValue(5.0)
        self.coverage_buffer.setSuffix(" km")
        self.coverage_buffer.setEnabled(False)
        crow.addWidget(self.coverage_buffer)
        crow.addStretch(1)
        form.addRow("", self._wrap(crow))
        return g

    def _wrap(self, layout):
        w = QWidget()
        w.setLayout(layout)
        return w

    # --- output --------------------------------------------------------- #
    def _build_output_group(self):
        g = QGroupBox("4. Output")
        lay = QVBoxLayout(g)
        row = QHBoxLayout()
        self.out_dir = QLineEdit()
        self.out_dir.setPlaceholderText("Output folder")
        btn = QPushButton("Browse…")
        btn.clicked.connect(self._pick_out_dir)
        row.addWidget(self.out_dir)
        row.addWidget(btn)
        lay.addLayout(row)

        prow = QHBoxLayout()
        prow.addWidget(QLabel("File prefix:"))
        self.prefix = QLineEdit("ztd")
        prow.addWidget(self.prefix)
        self.cb_load = QCheckBox("Load results into QGIS")
        self.cb_load.setChecked(True)
        prow.addWidget(self.cb_load)
        lay.addLayout(prow)
        return g

    # --- run ------------------------------------------------------------ #
    def _build_run_group(self):
        g = QGroupBox("5. Run / log")
        lay = QVBoxLayout(g)
        brow = QHBoxLayout()
        self.run_btn = QPushButton("Generate ZTD maps")
        self.run_btn.clicked.connect(self._on_run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._on_cancel)
        brow.addWidget(self.run_btn)
        brow.addWidget(self.cancel_btn)
        lay.addLayout(brow)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        lay.addWidget(self.progress)

        lay.addWidget(QLabel("Log:"))
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setMinimumHeight(160)
        lay.addWidget(self.log)
        return g

    # ------------------------------------------------------------------ #
    # File pickers / column loading
    # ------------------------------------------------------------------ #
    def _pick_ztd(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select ZTD table", "",
            "Tabular (*.csv *.txt *.tsv *.dat);;All files (*)")
        if path:
            self.ztd_path.setText(path)
            self._load_columns()

    def _pick_dem(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select DEM", "", "Rasters (*.tif *.tiff *.vrt *.asc);;All files (*)")
        if path:
            self.dem_path.setText(path)

    def _pick_out_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.out_dir.setText(path)

    def _delimiter(self):
        return {"auto": None, ",": ",", ";": ";", "tab": "\t", "space": " "}[
            self.delim.currentText()]

    def _load_columns(self):
        path = self.ztd_path.text().strip()
        if not path or not os.path.exists(path):
            self._msg("Pick a valid ZTD file first.", warn=True)
            return
        try:
            data, headers, mapping = io.read_ztd_table(
                path, delimiter=self._delimiter())
        except Exception as exc:
            try:
                import csv
                with open(path, encoding="utf-8-sig") as fh:
                    headers = next(csv.reader(fh, delimiter=self._delimiter() or ","))
                mapping = io.detect_columns(headers)
            except Exception:
                self._msg(f"Could not read file: {exc}", warn=True)
                return
        for key, cb in self._mapping_combos.items():
            cb.setEnabled(True)
            cb.clear()
            cb.addItems(headers)
            idx = mapping.get(key)
            if idx is not None and 0 <= idx < len(headers):
                cb.setCurrentIndex(idx)
        # populate the time-series epoch checklist (pre-check ZTD-like columns)
        epoch_idx = io.detect_epoch_columns(headers)
        self._populate_epochs(headers, epoch_idx)
        msg = f"Loaded {len(headers)} columns."
        if len(epoch_idx) > 1:
            msg += (f" Detected {len(epoch_idx)} ZTD epoch columns - tick "
                    "'Time series' to process them all.")
        self._append(msg + " Verify the lon/lat/height/ZTD mapping above.")

    # ------------------------------------------------------------------ #
    # Run / cancel
    # ------------------------------------------------------------------ #
    def _resolve_dem(self):
        """Return ('download'|'local', path) for the chosen DEM source."""
        if self.rb_download.isChecked():
            return "download", ""
        if self.rb_loaded.isChecked():
            lid = self.loaded_combo.currentData()
            if not lid:
                raise ValueError("Choose a loaded raster layer (or Refresh).")
            layer = QgsProject.instance().mapLayer(lid)
            if layer is None:
                raise ValueError("The selected raster layer is no longer loaded.")
            return "local", layer.source()
        # local file
        path = self.dem_path.text().strip()
        if not path or not os.path.exists(path):
            raise ValueError("Select a valid local DEM file.")
        return "local", path

    def _collect_params(self):
        ztd_path = self.ztd_path.text().strip()
        if not ztd_path or not os.path.exists(ztd_path):
            raise ValueError("Select a valid ZTD file.")
        if not any(cb.isEnabled() for cb in self._mapping_combos.values()):
            raise ValueError("Click 'Read columns' and map the columns first.")
        mapping = {k: cb.currentIndex() for k, cb in self._mapping_combos.items()}

        methods = []
        if self.cb_aitd.isChecked():
            methods.append("AITD")
        if self.cb_itd.isChecked():
            methods.append("ITD")
        if not methods:
            raise ValueError("Select at least one method (AITD / ITD).")

        products = []
        if self.cb_ztd.isChecked():
            products.append("ztd")
        if self.cb_strat.isChecked():
            products.append("stratified")
        if self.cb_turb.isChecked():
            products.append("turbulent")
        if not products:
            raise ValueError("Select at least one output map.")

        # epoch selection: time-series (multiple columns) or single ZTD column
        if self.cb_timeseries.isChecked():
            ztd_columns, epoch_names = self._checked_epochs()
            if not ztd_columns:
                raise ValueError("Tick at least one ZTD epoch column "
                                 "(or disable time series).")
        else:
            ztd_columns = [mapping["ztd"]]
            epoch_names = [self._mapping_combos["ztd"].currentText()]

        make_aps = self.cb_timeseries.isChecked() and self.cb_aps.isChecked()
        reference_epoch = self.ref_combo.currentText() if make_aps else None
        if make_aps:
            if reference_epoch not in epoch_names:
                raise ValueError("Choose a reference epoch that is also ticked.")
            if len(epoch_names) < 2:
                raise ValueError("APS maps need at least two epochs.")

        out_dir = self.out_dir.text().strip()
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(ztd_path), "ztd_output")

        dem_mode, dem_path = self._resolve_dem()

        return {
            "ztd_path": ztd_path,
            "mapping": mapping,
            "ztd_columns": ztd_columns,
            "epoch_names": epoch_names,
            "make_aps": make_aps,
            "reference_epoch": reference_epoch,
            "delimiter": self._delimiter(),
            "dem_mode": dem_mode,
            "dem_path": dem_path,
            "dem_resolution_m": float(self.dem_res.value()),
            "margin_deg": float(self.margin.value()),
            "radius_km": float(self.radius.value()),
            "max_iter": int(self.max_iter.value()),
            "methods": methods,
            "products": products,
            "fill_nearest": self.cb_fill.isChecked(),
            "sea_floor": 0.0 if self.cb_sea.isChecked() else None,
            "flag_extrapolation": self.cb_flag.isChecked(),
            "clip_coverage": self.cb_clip.isChecked(),
            "coverage_buffer_km": float(self.coverage_buffer.value()),
            "output_dir": out_dir,
            "prefix": self.prefix.text().strip() or "ztd",
            "add_dem_layer": dem_mode == "download",
        }

    def _planned_output_paths(self, params):
        """Paths the task will (over)write - used to release QGIS file locks."""
        return [os.path.normcase(os.path.abspath(p))
                for p in planned_output_files(params)]

    def _free_output_files(self, params):
        """Unload any project layers pointing at files we are about to rewrite.

        On Windows an open raster layer locks its file, so GDAL cannot overwrite
        it (the 'Deleting ... Permission denied' error). Removing the layer first
        releases the lock.
        """
        targets = set(self._planned_output_paths(params))
        proj = QgsProject.instance()
        to_remove = []
        for lid, layer in proj.mapLayers().items():
            try:
                # vector providers append e.g. "|layername=stations" - strip it
                raw = layer.source().split("|")[0]
                src = os.path.normcase(os.path.abspath(raw))
            except Exception:
                continue
            if src in targets:
                to_remove.append(lid)
        if to_remove:
            proj.removeMapLayers(to_remove)
            self._append(f"Released {len(to_remove)} previously loaded "
                         "output layer(s) before overwriting.")

    def _on_run(self):
        try:
            params = self._collect_params()
            self._free_output_files(params)
        except ValueError as exc:
            self._msg(str(exc), warn=True)
            return

        self.progress.setValue(0)
        self.log.clear()
        self._append("Starting ZTD-AITD generation…")
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)

        self.task = AitdTask(params)
        self.task.progressChanged.connect(
            lambda v: self.progress.setValue(int(v)))
        self.task.taskCompleted.connect(self._on_done)
        self.task.taskTerminated.connect(self._on_terminated)
        QgsApplication.taskManager().addTask(self.task)

    def _on_cancel(self):
        if self.task is not None:
            self.task.cancel()
            self._append("Cancelling…")

    def _on_done(self):
        self._append("Generation finished.")
        for title, (mean, std) in (self.task.summary or {}).items():
            self._append(f"  {title}: mean={mean:.4f} m, std={std*1000:.1f} mm")
        if self.cb_load.isChecked():
            self._load_outputs(self.task.outputs)
        self._reset_buttons()
        self._msg("ZTD maps generated successfully.", warn=False)

    def _on_terminated(self):
        if self.task is not None and self.task.exception is not None:
            self._append("FAILED: " + str(self.task.exception))
            self._msg("Generation failed - see the log panel.", warn=True)
        else:
            self._append("Cancelled.")
        self._reset_buttons()

    def _reset_buttons(self):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)

    # ------------------------------------------------------------------ #
    # Result loading & styling
    # ------------------------------------------------------------------ #
    def _load_outputs(self, outputs):
        proj = QgsProject.instance()
        for title, path in outputs.items():
            if path.lower().endswith((".gpkg", ".shp")):
                layer = QgsVectorLayer(path, title, "ogr")
                if not layer.isValid():
                    self._append(f"  (could not load {path})")
                    continue
                self._style_points(layer)
            else:
                layer = QgsRasterLayer(path, title)
                if not layer.isValid():
                    self._append(f"  (could not load {path})")
                    continue
                self._style(layer, title)
            proj.addMapLayer(layer)
            self._append(f"  loaded: {title}")
        self._refresh_loaded_rasters()

    def _style_points(self, layer):
        """Render GNSS stations as clear black-outlined circles, on top."""
        try:
            symbol = QgsMarkerSymbol.createSimple({
                "name": "circle",
                "color": "255,255,0,255",      # yellow fill (like the paper)
                "outline_color": "0,0,0,255",
                "outline_width": "0.4",
                "size": "2.6",
            })
            layer.renderer().setSymbol(symbol)
        except Exception:
            pass

    def _style(self, layer, title):
        ramp_name = "Spectral"
        invert = True
        low = title.lower()
        if "extrapolation mask" in low:
            self._style_mask(layer)
            return
        is_aps = "aps" in low
        if low.startswith("dem"):
            ramp_name, invert = "Greys", False
        elif is_aps:
            # difference maps: diverging ramp centred on zero
            ramp_name, invert = "RdBu", False
        try:
            stats = layer.dataProvider().bandStatistics(1)
            vmin, vmax = stats.minimumValue, stats.maximumValue
            if is_aps and self.cb_aps_symmetric.isChecked():
                m = max(abs(vmin), abs(vmax)) or 1e-6
                vmin, vmax = -m, m
            if vmin == vmax:
                vmax = vmin + 1e-6
            ramp = QgsStyle.defaultStyle().colorRamp(ramp_name)
            if ramp is None:
                return
            if invert:
                ramp.invert()
            shader_fn = QgsColorRampShader(vmin, vmax)
            shader_fn.setColorRampType(QgsColorRampShader.Interpolated)
            shader_fn.setSourceColorRamp(ramp)
            n = 8
            items = []
            for i in range(n + 1):
                v = vmin + (vmax - vmin) * i / n
                items.append(QgsColorRampShader.ColorRampItem(
                    v, ramp.color(i / n), f"{v:.3f}"))
            shader_fn.setColorRampItemList(items)
            shader = QgsRasterShader()
            shader.setRasterShaderFunction(shader_fn)
            renderer = QgsSingleBandPseudoColorRenderer(
                layer.dataProvider(), 1, shader)
            layer.setRenderer(renderer)
        except Exception:
            pass

    def _style_mask(self, layer):
        """Render the extrapolation mask as a semi-transparent red overlay.

        The raster has data (value 1) only over extrapolated terrain; everywhere
        else is NoData and therefore fully transparent.
        """
        try:
            shader_fn = QgsColorRampShader()
            shader_fn.setColorRampType(QgsColorRampShader.Interpolated)
            shader_fn.setColorRampItemList([
                QgsColorRampShader.ColorRampItem(
                    1.0, QColor(214, 40, 40),
                    "extrapolated (above highest station)")
            ])
            shader = QgsRasterShader()
            shader.setRasterShaderFunction(shader_fn)
            renderer = QgsSingleBandPseudoColorRenderer(
                layer.dataProvider(), 1, shader)
            layer.setRenderer(renderer)
            layer.setOpacity(0.5)
        except Exception:
            pass

    # ------------------------------------------------------------------ #
    def _append(self, text):
        self.log.appendPlainText(text)

    def _msg(self, text, warn=False):
        level = Qgis.Warning if warn else Qgis.Info
        self.iface.messageBar().pushMessage("ZTD-AITD", text, level=level, duration=5)
        self._append(text)
