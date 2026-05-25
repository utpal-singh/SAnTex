from PyQt5.QtWidgets import (
    QMainWindow, QTabWidget, QStatusBar, QMenuBar, QAction, QFileDialog,
    QMessageBox, QWidget
)
from PyQt5.QtCore import Qt

from backend.ebsd_backend import EBSDBackend
from backend.material_backend import MaterialBackend
from backend.anisotropy_backend import AnisotropyBackend
from backend.grains_backend import GrainsBackend
from backend.odf_backend import ODFBackend


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SAnTex — Seismic Anisotropy Toolbox")
        self.setMinimumSize(1200, 800)

        # Shared backend instances
        self.ebsd_backend       = EBSDBackend()
        self.material_backend   = MaterialBackend()
        self.anisotropy_backend = AnisotropyBackend()
        self.grains_backend     = GrainsBackend()
        self.odf_backend        = ODFBackend()

        self._build_menu()
        self._build_tabs()
        self._build_statusbar()

    # ------------------------------------------------------------------
    def _build_menu(self):
        menubar: QMenuBar = self.menuBar()

        file_menu = menubar.addMenu("File")
        open_act = QAction("Open EBSD file…", self)
        open_act.setShortcut("Ctrl+O")
        open_act.triggered.connect(self._open_ebsd)
        file_menu.addAction(open_act)
        file_menu.addSeparator()
        quit_act = QAction("Quit", self)
        quit_act.setShortcut("Ctrl+Q")
        quit_act.triggered.connect(self.close)
        file_menu.addAction(quit_act)

        help_menu = menubar.addMenu("Help")
        about_act = QAction("About SAnTex", self)
        about_act.triggered.connect(self._show_about)
        help_menu.addAction(about_act)

    def _build_tabs(self):
        from frontend.tabs.material_tab import MaterialTab
        from frontend.tabs.anisotropy_tab import AnisotropyTab
        from frontend.tabs.ebsd_tab import EBSDTab
        from frontend.tabs.modal_rock_tab import ModalRockTab
        from frontend.tabs.grains_tab import GrainsTab
        from frontend.tabs.odf_tab import ODFTab

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)

        self.material_tab   = MaterialTab(self.material_backend, self.anisotropy_backend)
        self.anisotropy_tab = AnisotropyTab(self.anisotropy_backend)
        self.ebsd_tab       = EBSDTab(self.ebsd_backend, self.material_backend,
                                      self.anisotropy_backend)
        self.modal_rock_tab = ModalRockTab(self.material_backend, self.anisotropy_backend)
        self.grains_tab     = GrainsTab(self.grains_backend, self.ebsd_backend,
                                        ebsd_tab=self.ebsd_tab)
        self.odf_tab        = ODFTab(self.odf_backend, self.ebsd_backend)

        self.tabs.addTab(self.material_tab,   "Material Database")
        self.tabs.addTab(self.anisotropy_tab, "Anisotropy")
        self.tabs.addTab(self.ebsd_tab,       "EBSD")
        self.tabs.addTab(self.modal_rock_tab, "Modal Rock")
        self.tabs.addTab(self.grains_tab,     "Grains")
        self.tabs.addTab(self.odf_tab,        "ODF & PF")

        self.setCentralWidget(self.tabs)

    def _build_statusbar(self):
        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self.status.showMessage("Ready")

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _open_ebsd(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open EBSD file", "",
            "EBSD files (*.ctf *.ang *.osc *.txt);;All files (*)"
        )
        if not path:
            return
        try:
            self.ebsd_backend.load(path)
            n = self.ebsd_backend.n_points
            self.status.showMessage(f"Loaded {path}  —  {n} points")
            self.ebsd_tab.on_file_loaded()
            self.tabs.setCurrentWidget(self.ebsd_tab)
        except Exception as e:
            QMessageBox.critical(self, "Load error", str(e))

    def _show_about(self):
        QMessageBox.information(
            self, "About SAnTex",
            "SAnTex — Seismic Anisotropy Toolbox\n\n"
            "PyQt5 GUI for computing and visualising seismic anisotropy\n"
            "from single-crystal elastic constants and EBSD orientation data.\n\n"
            "Backend: santex Python library\n"
            "3-D visualisation: PyVista / VTK\n"
            "2-D plots: Matplotlib"
        )
