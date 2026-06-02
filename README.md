# ZTD-AITD Map Generator — QGIS Plugin

Generate **Zenith Tropospheric Delay (ZTD) maps** from sparse GNSS observations
using the **Iterative Tropospheric Decomposition (ITD)** strategy and its
**Augmented** variant (**AITD**).

This plugin is a full, from-scratch implementation of:

> X. Song, G. Venuti, A. V. Monti-Guarnieri, M. Manzoni (2025).
> *Augmented iterative tropospheric decomposition strategy for GNSS-based
> zenith tropospheric delay map generation.*
> **Environmental Modelling and Software 194, 106669.**
> https://doi.org/10.1016/j.envsoft.2025.106669
> (Politecnico di Milano)

It reads a table of GNSS ZTDs, automatically fetches (or loads) a matching DEM,
and produces four raster products: the **ZTD map**, its **stratified** and
**turbulent** components, and the **DEM** itself.

---

## 1. Scientific background

At any point `P(x, h)` the ZTD is modelled as the sum of two parts (Eq. 1):

```
ZTD(x, h) = S(h) + T(x) + noise
```

* **Stratified** component — height-dependent, an exponential model (Eq. 2):

  ```
  S(h) = L · exp(−B · h)          (h scaled to [0, 1])
  ```

* **Turbulent** component — horizontally correlated, obtained by inverse-distance
  weighting of the residuals.

A moving **decorrelation window** of radius `r_G` selects the reference points
used at each location. Inside the window an **iterative decomposition**
alternates between fitting the stratified model and estimating the turbulent
residual until convergence (typically 4–6 iterations; Eqs. 3–7).

| Method | Where decomposition runs | Weighting | Notes |
|--------|--------------------------|-----------|-------|
| **ITD** (Sec. 2.1) | once **per DEM knot** | classic IDW `1/d²` (Eq. 5) | original method; can produce artefacts on sparse GNSS data |
| **AITD** (Sec. 2.2) | once **per observation point**, then `L, B` and residuals are interpolated to the DEM | Modified IDW `1/(1 + d²/α)` (Eq. 8) | faster; smooths `L, B`; removes interpolation peaks |

The MIDW parameter `α` is set so the weight equals `0.01` at the decorrelation
radius: `α = r_G² / 99`.

> **Exact ITD optimisation.** The iterative decomposition at a DEM knot depends
> *only* on the set of reference points inside the window, never on the knot’s
> exact position. The plugin therefore groups all knots sharing the same
> neighbour set and decomposes once per group — an exact speed-up that keeps the
> original ITD tractable on million-knot grids.

---

## 2. Features

- Reads ZTD from **CSV / TSV / whitespace-delimited** tables with automatic
  column detection (`lon, lat, height, ZTD`) and a manual column-mapping UI.
- **Automatic DEM download** — keyless, from the
  [AWS Terrain Tiles](https://registry.opendata.aws/terrain-tiles/) open dataset
  (no API key, no registration). Or use a **local DEM** file, or a **DEM already
  loaded in QGIS**.
- **Sea-level clamp** — negative elevations in the terrain tiles (sea/lake
  bathymetry) are clamped to 0 m by default, preventing unrealistic ZTD
  inflation below sea level.
- **Coverage clip** (optional) — restricts output to the convex hull of the
  stations (plus a configurable buffer), removing the unreliable border/corner
  zone *beyond* the network where horizontal extrapolation causes artifacts.
- **Extrapolation mask** (optional) — flags pixels whose elevation exceeds the
  highest GNSS station, where the stratified model is extrapolated rather than
  interpolated. Written as a separate raster (data only over flagged terrain)
  and loaded as a semi-transparent red overlay, so reliable vs uncertain zones
  are clear in figures.
- Builds the DEM bounding box automatically from the station coordinates
  (with a configurable border margin).
- Implements **both ITD and AITD**, selectable and comparable.
- **Time series / multi-epoch** — process several ZTD epoch columns in one run
  (e.g. `ZTD1, ZTD2, …`), with automatic handling of missing values per epoch,
  and optional **APS-like difference maps** (each epoch minus a reference
  epoch), as the paper does to compare against InSAR.
- Outputs **ZTD map, stratified map, turbulent map and DEM** as styled GeoTIFFs,
  plus the **GNSS stations** as a point layer (GeoPackage) for context.
- Runs on a **background thread** (`QgsTask`) with a progress bar and log — the
  QGIS UI stays responsive, and long runs can be cancelled.
- Automatic, sensible raster styling on load (Spectral ramp for delays, grey for
  the DEM).

---

## 3. Installation

The plugin already lives in your QGIS plugin folder:

```
…/QGIS3/profiles/default/python/plugins/ZTD/
```

1. Start QGIS (3.22 or newer).
2. `Plugins ▸ Manage and Install Plugins… ▸ Installed`.
3. Enable **“ZTD-AITD Map Generator”**.
4. A toolbar icon and a `Plugins ▸ ZTD-AITD` menu entry appear; click either to
   open the control window (an independent, resizable, minimizable window).

If you edit the code while QGIS is running, use the **Plugin Reloader** plugin
to reload without restarting.

### Requirements
NumPy, GDAL/OGR and (optionally) SciPy — all bundled with the standard QGIS
installation. SciPy, when present, accelerates neighbour search; the code falls
back to NumPy otherwise.

---

## 4. Input format

A delimited text file with one row per GNSS station, e.g.
[`sample_data/ztd_sample.csv`](sample_data/ztd_sample.csv):

```csv
station,lon,lat,height,ztd
MILA,9.230,45.478,140.0,2.3412
COMO,9.085,45.808,290.0,2.3105
...
```

| Column | Meaning | Unit |
|--------|---------|------|
| `lon`  | longitude | decimal degrees (WGS84) |
| `lat`  | latitude  | decimal degrees (WGS84) |
| `height` | station height | metres |
| `ztd`  | zenith tropospheric delay | metres |

Extra columns (e.g. a station name) are ignored. Column names are auto-detected;
if detection fails you can map columns manually in the panel.

**Multiple epochs.** A file may hold several ZTD columns (e.g. `ZTD1, ZTD2, …`,
one per epoch). Tick *"Time series"* and select which columns to process; blank
cells are treated as a missing station for that epoch only.

---

## 5. Usage

1. **GNSS ZTD observations** — browse to your file, pick a delimiter (or leave
   *auto*), click **Read columns** and confirm the `lon / lat / height / ZTD`
   mapping. For multi-epoch files, tick **Time series**, choose the ZTD epoch
   columns, and (optionally) enable **APS difference maps** with a reference
   epoch.
2. **DEM** — choose one of: *Download automatically* (set a target resolution,
   default 90 m, matching the paper), *Use a DEM already loaded in QGIS* (pick
   from the dropdown; press *Refresh* to update the list), or *Use a local DEM
   file*.
3. **Interpolation parameters**:
   - *Decorrelation radius* (km) — the moving-window radius `r_G` (default 100).
   - *Max iterations* — decomposition iteration cap (default 10).
   - *Method* — **AITD**, **ITD**, or both.
   - *Output maps* — ZTD, stratified, turbulent.
   - *Clamp DEM below sea level* — recommended (on by default).
   - *Flag pixels above the highest station* — adds the extrapolation-mask
     overlay (off by default).
   - *Clip output to station coverage* — masks everything outside the stations'
     convex hull plus the chosen buffer (off by default; buffer 5 km); removes
     border/corner extrapolation artifacts.
4. **Output** — choose a folder and a file prefix.
5. Click **Generate ZTD maps**. Watch the progress bar / log; results load
   automatically and styled.

---

## 6. Outputs

For each selected method and product a GeoTIFF (EPSG:4326, float32) is written:

```
<prefix>_dem.tif                 the DEM used
<prefix>_aitd_ztd.tif            AITD total ZTD map
<prefix>_aitd_stratified.tif     AITD stratified component  S(h)
<prefix>_aitd_turbulent.tif      AITD turbulent component   T(x)
<prefix>_itd_ztd.tif             ITD total ZTD map  (if ITD selected)
<prefix>_itd_stratified.tif      …
<prefix>_itd_turbulent.tif       …
<prefix>_extrapolation_mask.tif  pixels above the highest station (if enabled)
<prefix>_stations.gpkg          GNSS station points (lon, lat, height, ZTD/epoch)
```

With **multiple epochs**, each map name gains an epoch suffix
(`<prefix>_aitd_ztd_ztd1.tif`, `…_ztd2.tif`, …), and APS difference maps are
written as `<prefix>_<method>_ztd_aps_<epoch>.tif` (epoch − reference epoch,
styled with a diverging red/blue ramp). The *"APS legend symmetric about zero"*
checkbox (on by default) keeps white = no change; unticking it stretches the
ramp to the actual data range instead.

Delays are in **metres**. Knots with no DEM data (e.g. clipped sea) are NoData.

---

## 7. Code structure

```
ZTD/
├── metadata.txt              plugin manifest
├── __init__.py               classFactory entry point
├── ztd_aitd_plugin.py        main plugin (toolbar / menu / dock)
├── icon.png
├── gui/
│   └── ztd_window.py         the control window + result styling
├── core/                     pure-Python, PyQGIS-free, unit-testable
│   ├── decomposition.py      ITD + AITD (the paper’s equations)
│   ├── dem_downloader.py     AWS Terrain Tiles → mosaicked GeoTIFF
│   ├── io_utils.py           table reading, bbox, raster read/write
│   └── projection.py         lon/lat → metric UTM (for distances)
├── workers/
│   └── aitd_task.py          QgsTask orchestrating the pipeline
├── tests/
│   └── test_decomposition.py synthetic-truth validation
└── sample_data/
    └── ztd_sample.csv        32 Italian-style GNSS stations
```

The `core` package contains **no PyQGIS imports**, so the science can be tested
with plain Python + NumPy.

---

## 8. Testing

Unit test (synthetic truth field, leave-one-out validation):

```bash
python tests/test_decomposition.py
```

It checks that both ITD and AITD recover a known stratified-plus-turbulent field
with a grid RMS below 1 cm and a leave-one-out station RMS below 1.5 cm —
consistent with the 5.7–14.9 mm reported in the paper.

The `tests/_*_check.py` scripts are developer integration checks that exercise
the DEM download, projection and raster writing using QGIS’s own Python:

```powershell
& "C:\Program Files\QGIS 3.42.2\bin\python-qgis.bat" tests\_e2e_check.py
```

---

## 9. Notes and limitations

- **Distances** are computed in an auto-selected **UTM zone** (metres). For very
  large areas spanning multiple zones, accuracy near the edges may degrade
  slightly.
- **Heights** are scaled to `[0, 1]` over the combined station + DEM range, as in
  the paper; the difference between EGM96/EGM2008 height systems affects the ZTD
  by < 1 mm and is neglected.
- The **ITD** method is heavier than AITD; on very fine DEMs prefer AITD (the
  paper’s recommended method) or use a coarser DEM resolution for ITD.
- DEM tiles come from the AWS Terrain Tiles open dataset (SRTM/ASTER/others,
  depending on location). Attribute the source per its
  [licence](https://github.com/tilezen/joerd/blob/master/docs/attribution.md).

---

## 10. References

- Song, Venuti, Monti-Guarnieri, Manzoni (2025), *Environmental Modelling and
  Software* 194, 106669. https://doi.org/10.1016/j.envsoft.2025.106669
- Yu et al. (2017, 2018a, 2018b) — original ITD strategy and GACOS.
- AWS Terrain Tiles — https://registry.opendata.aws/terrain-tiles/

## Author

**Rohollah Naeijian** — rohollah.naeijian@mail.polimi.it
Developed for the **Geoinformatics Project** course, **Politecnico di Milano**.
