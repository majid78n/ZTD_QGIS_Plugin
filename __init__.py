"""ZTD-AITD Map Generator - QGIS plugin entry point.

QGIS calls :func:`classFactory` to instantiate the plugin.
"""


def classFactory(iface):  # noqa: N802 (QGIS-mandated name)
    from .ztd_aitd_plugin import ZtdAitdPlugin
    return ZtdAitdPlugin(iface)
