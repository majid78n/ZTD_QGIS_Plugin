"""Import every module of the plugin to catch syntax / import errors.

Run with QGIS's python so PyQGIS is importable.
"""
import os, sys, importlib
# make the plugin importable as a package named 'ZTD'
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
parent = os.path.dirname(plugin_dir)
sys.path.insert(0, parent)
pkg = os.path.basename(plugin_dir)

mods = [
    f"{pkg}.core.decomposition",
    f"{pkg}.core.io_utils",
    f"{pkg}.core.projection",
    f"{pkg}.core.dem_downloader",
    f"{pkg}.workers.aitd_task",
    f"{pkg}.gui.ztd_window",
    f"{pkg}.ztd_aitd_plugin",
    f"{pkg}",
]
for m in mods:
    importlib.import_module(m)
    print("OK", m)
print("ALL IMPORTS OK")
