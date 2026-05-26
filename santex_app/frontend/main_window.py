from PyQt5.QtWidgets import (
    QMainWindow, QTabWidget, QStatusBar, QMenuBar, QAction, QFileDialog,
    QMessageBox, QWidget
)
from PyQt5.QtCore import Qt

from backend.ebsd_backend  import EBSDBackend
from backend.ebsd_project  import EBSDProject
from backend.material_backend  import MaterialBackend
from backend.anisotropy_backend import AnisotropyBackend
from backend.grains_backend import GrainsBackend
from backend.odf_backend import ODFBackend


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SAnTex — Seismic Anisotropy Toolbox")
        self.setMinimumSize(1200, 800)

        # ── Shared backend instances ──────────────────────────────────
        self.project            = EBSDProject(max_in_ram=2)
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

        # ── File menu ──────────────────────────────────────────────────
        file_menu = menubar.addMenu("File")

        # Project submenu
        proj_menu = file_menu.addMenu("Project")
        new_proj_act   = QAction("New project",    self); new_proj_act.setShortcut("Ctrl+Shift+N")
        open_proj_act  = QAction("Open project…",  self); open_proj_act.setShortcut("Ctrl+Shift+O")
        save_proj_act  = QAction("Save project",   self); save_proj_act.setShortcut("Ctrl+Shift+S")
        saveas_proj_act = QAction("Save project as…", self)
        new_proj_act.triggered.connect(lambda: self.project_tab._new_project())
        open_proj_act.triggered.connect(lambda: self.project_tab._open_project())
        save_proj_act.triggered.connect(lambda: self.project_tab._save_project())
        saveas_proj_act.triggered.connect(lambda: self.project_tab._saveas_project())
        for act in (new_proj_act, open_proj_act, save_proj_act, saveas_proj_act):
            proj_menu.addAction(act)

        file_menu.addSeparator()

        # Add EBSD files to project (shortcut)
        add_act = QAction("Add EBSD file(s) to project…", self)
        add_act.setShortcut("Ctrl+O")
        add_act.triggered.connect(self._add_ebsd_files)
        file_menu.addAction(add_act)

        # Quick-open (single file, load into EBSD tab directly)
        quick_act = QAction("Quick open EBSD file…", self)
        quick_act.setShortcut("Ctrl+Shift+Q")
        quick_act.triggered.connect(self._quick_open_ebsd)
        file_menu.addAction(quick_act)

        file_menu.addSeparator()
        quit_act = QAction("Quit", self)
        quit_act.setShortcut("Ctrl+Q")
        quit_act.triggered.connect(self.close)
        file_menu.addAction(quit_act)

        # ── Help menu ──────────────────────────────────────────────────
        help_menu = menubar.addMenu("Help")
        about_act = QAction("About SAnTex", self)
        about_act.triggered.connect(self._show_about)
        help_menu.addAction(about_act)

    def _build_tabs(self):
        from frontend.tabs.material_tab   import MaterialTab
        from frontend.tabs.anisotropy_tab import AnisotropyTab
        from frontend.tabs.ebsd_tab       import EBSDTab
        from frontend.tabs.modal_rock_tab import ModalRockTab
        from frontend.tabs.grains_tab     import GrainsTab
        from frontend.tabs.odf_tab        import ODFTab
        from frontend.tabs.project_tab    import ProjectTab

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)

        # Project tab (first — it's the entry point for multi-file work)
        self.project_tab = ProjectTab(
            self.project, self.ebsd_backend, parent=self
        )
        self.project_tab.open_dataset.connect(self._open_dataset_in_ebsd_tab)
        self.project_tab.project_changed.connect(self._on_project_changed)

        self.material_tab   = MaterialTab(self.material_backend, self.anisotropy_backend)
        self.anisotropy_tab = AnisotropyTab(self.anisotropy_backend)
        self.ebsd_tab       = EBSDTab(
            self.ebsd_backend, self.material_backend, self.anisotropy_backend
        )
        self.modal_rock_tab = ModalRockTab(self.material_backend, self.anisotropy_backend)
        self.grains_tab     = GrainsTab(
            self.grains_backend, self.ebsd_backend, ebsd_tab=self.ebsd_tab
        )
        self.odf_tab        = ODFTab(self.odf_backend, self.ebsd_backend)

        self.tabs.addTab(self.project_tab,    "📁 Project")
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
        self.status.showMessage("Ready  —  Open a project or add EBSD files to get started.")

    # ------------------------------------------------------------------
    # File loading slots
    # ------------------------------------------------------------------

    def _add_ebsd_files(self):
        """Add one or more CTF files to the project (no immediate full load)."""
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add EBSD files to project", "",
            "EBSD files (*.ctf *.ang *.osc);;All files (*)",
        )
        if not paths:
            return
        added = 0
        for path in paths:
            try:
                self.project.add_file(path)
                added += 1
            except Exception as e:
                QMessageBox.warning(self, "Add file", f"Could not add {path}:\n{e}")
        if added:
            self.project_tab._refresh_all()
            self.project_tab._update_mem_label()
            self.project_tab.project_changed.emit()
            # Switch to project tab so user can see the list
            self.tabs.setCurrentWidget(self.project_tab)
            self.status.showMessage(
                f"Added {added} file(s) to project  —  double-click a row to open in EBSD tab."
            )

    def _quick_open_ebsd(self):
        """Open a single EBSD file directly into the EBSD tab (bypasses project LRU)."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Quick-open EBSD file", "",
            "EBSD files (*.ctf *.ang *.osc);;All files (*)",
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

    def _open_dataset_in_ebsd_tab(self, dataset_id: str):
        """Load *dataset_id* from the project cache into the EBSD tab."""
        meta = self.project.get_meta(dataset_id)
        if meta is None:
            QMessageBox.warning(self, "Open dataset", "Dataset not found in project.")
            return
        try:
            self.ebsd_backend.load_from_project(self.project, dataset_id)
            n = self.ebsd_backend.n_points
            self.status.showMessage(
                f"Loaded '{meta.name}'  —  {n:,} points"
                + ("  [from cache]" if self.project.is_cached(dataset_id) else "  [from disk]")
            )
            self.ebsd_tab.on_file_loaded()
            self.tabs.setCurrentWidget(self.ebsd_tab)
            # Refresh cache badges in the project panel
            self.project_tab.refresh_cache_badges()
        except Exception as e:
            QMessageBox.critical(self, "Load error", str(e))

    def _on_project_changed(self):
        """Called whenever datasets are added, removed, or toggled."""
        total = len(self.project)
        active = len(self.project.active_datasets())
        self.status.showMessage(
            f"Project: {total} dataset(s)  ·  {active} active"
        )

    # ------------------------------------------------------------------
    def _show_about(self):
        QMessageBox.information(
            self, "About SAnTex",
            "SAnTex — Seismic Anisotropy Toolbox\n\n"
            "PyQt5 GUI for computing and visualising seismic anisotropy\n"
            "from single-crystal elastic constants and EBSD orientation data.\n\n"
            "Backend: santex Python library\n"
            "3-D visualisation: PyVista / VTK\n"
            "2-D plots: Plotly"
        )
