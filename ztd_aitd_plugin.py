"""Main plugin class - registers the toolbar/menu action and the dock panel."""

from __future__ import annotations

import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .gui.ztd_window import ZtdMainWindow

PLUGIN_DIR = os.path.dirname(__file__)
MENU = "&ZTD-AITD"


class ZtdAitdPlugin:
    """QGIS plugin: ZTD map generation via (Augmented) ITD."""

    def __init__(self, iface):
        self.iface = iface
        self.action = None
        self.window = None

    # ------------------------------------------------------------------ #
    def initGui(self):  # noqa: N802 (QGIS-mandated name)
        icon_path = os.path.join(PLUGIN_DIR, "icon.png")
        icon = QIcon(icon_path) if os.path.exists(icon_path) else QIcon()
        self.action = QAction(icon, "ZTD-AITD Map Generator", self.iface.mainWindow())
        self.action.triggered.connect(self.show_window)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu(MENU, self.action)

    def unload(self):
        if self.window is not None:
            self.window.close()
            self.window.deleteLater()
            self.window = None
        if self.action is not None:
            self.iface.removePluginMenu(MENU, self.action)
            self.iface.removeToolBarIcon(self.action)
            self.action = None

    # ------------------------------------------------------------------ #
    def show_window(self):
        if self.window is None:
            # No parent -> a fully independent, minimizable top-level window.
            self.window = ZtdMainWindow(self.iface, None)
        self.window.show()
        self.window.raise_()
        self.window.activateWindow()
