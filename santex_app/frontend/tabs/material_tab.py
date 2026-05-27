"""
MaterialTab — browse the mineral database, pick P/T conditions,
view the Voigt stiffness matrix, scalar anisotropy metrics,
2-D stereonets for *all* velocity scalars, and the 3-D velocity surface.
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QTableWidget, QTableWidgetItem, QAbstractItemView,
    QDoubleSpinBox, QPushButton, QSplitter, QHeaderView,
    QComboBox, QFormLayout, QSizePolicy, QTabWidget,
    QProgressBar,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES
from frontend.widgets.pyvista_widget import PyVistaWidget
from frontend.tabs.ebsd_tab import _PlotOptions, _SCALAR_LABELS
from frontend.tabs._stereonet import make_stereonet_figure
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Worker thread — velocity surface computation can be slow at high resolution
# ---------------------------------------------------------------------------

class _SurfaceWorker(QThread):
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, backend, n_theta: int = 60, n_phi: int = 120):
        super().__init__()
        self.backend = backend
        self.n_theta = n_theta
        self.n_phi   = n_phi

    def run(self):
        try:
            surf = self.backend.compute_velocity_surface(self.n_theta, self.n_phi)
            if surf is None:
                self.error.emit("compute_velocity_surface returned None.")
                return
            self.finished.emit(surf)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# MaterialTab
# ---------------------------------------------------------------------------

class MaterialTab(QWidget):
    """Browse the mineral database, pick P/T conditions, view elastic tensor
    and all anisotropy metrics — stereonet + 3-D surface."""

    def __init__(self, material_backend, anisotropy_backend, parent=None):
        super().__init__(parent)
        self.mb = material_backend
        self.ab = anisotropy_backend
        self._surf_data: dict | None = None
        self._worker: _SurfaceWorker | None = None
        self._build_ui()
        self._populate_table()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ── Left: material list + PT controls ─────────────────────────
        left = QWidget()
        lv = QVBoxLayout(left)

        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(
            ["Phase", "Crystal System", "Primary Phase", "ρ (g/cm³)"]
        )
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.selectionModel().selectionChanged.connect(self._on_select)
        lv.addWidget(QLabel("Mineral database:"))
        lv.addWidget(self.table)

        pt_box = QGroupBox("Conditions")
        pt_form = QFormLayout(pt_box)

        self.pressure_spin = QDoubleSpinBox()
        self.pressure_spin.setRange(0, 300)
        self.pressure_spin.setSuffix(" GPa")
        self.pressure_spin.setDecimals(2)
        self.pressure_spin.setValue(0.0)
        pt_form.addRow("Pressure:", self.pressure_spin)

        self.temp_spin = QDoubleSpinBox()
        self.temp_spin.setRange(0, 3000)
        self.temp_spin.setSuffix(" K")
        self.temp_spin.setDecimals(1)
        self.temp_spin.setValue(300.0)
        pt_form.addRow("Temperature:", self.temp_spin)

        lv.addWidget(pt_box)

        self.calc_btn = QPushButton("Compute anisotropy")
        self.calc_btn.clicked.connect(self._compute)
        lv.addWidget(self.calc_btn)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        lv.addWidget(self.progress)

        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        lv.addWidget(self.info_label)

        splitter.addWidget(left)

        # ── Right: Voigt matrix + metrics + visualisation ──────────────
        right = QWidget()
        rv = QVBoxLayout(right)

        # Voigt matrix
        rv.addWidget(QLabel("Voigt stiffness matrix (GPa):"))
        self.voigt_table = QTableWidget(6, 6)
        self.voigt_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.voigt_table.setMaximumHeight(200)
        self.voigt_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.voigt_table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        rv.addWidget(self.voigt_table)

        # Anisotropy metrics table
        rv.addWidget(QLabel("Anisotropy metrics:"))
        self.anis_table = QTableWidget()
        self.anis_table.setColumnCount(2)
        self.anis_table.setHorizontalHeaderLabels(["Quantity", "Value"])
        self.anis_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.anis_table.setMaximumHeight(160)
        self.anis_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        rv.addWidget(self.anis_table)

        # Plot options — includes scalar selector (vp / vs1 / vs2 / avs …)
        self._sn_opts = _PlotOptions(
            "Plot options",
            default_cmap="RdBu_r",
            default_pt_size=2,
            show_scalar=True,   # ← enables vp/vs1/vs2/avs/vpvs1/vpvs2 combo
        )
        self._sn_opts.scalar_combo.currentIndexChanged.connect(
            lambda _: self._replot_if_ready()
        )
        self._sn_opts.style_combo.currentTextChanged.connect(
            lambda _: self._replot_if_ready()
        )
        self._sn_opts.cmap_combo.currentTextChanged.connect(
            lambda _: self._replot_if_ready()
        )
        rv.addWidget(self._sn_opts)

        # Visualisation tabs: 2-D stereonet  |  3-D surface
        self.vis_tabs = QTabWidget()

        # 2-D stereonet
        sn2d = QWidget()
        sn2d_v = QVBoxLayout(sn2d)
        sn2d_v.setContentsMargins(0, 0, 0, 0)
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        self._browser_btn = QPushButton("🌐 Open in browser")
        self._browser_btn.setToolTip("Open this stereonet in the system browser")
        btn_row.addWidget(self._browser_btn)
        sn2d_v.addLayout(btn_row)
        self.mpl = PlotlyWidget()
        self._browser_btn.clicked.connect(self.mpl.open_in_browser)
        sn2d_v.addWidget(self.mpl)
        self.vis_tabs.addTab(sn2d, "2-D stereonet")

        # 3-D surface
        self.pv3d = PyVistaWidget()
        self.vis_tabs.addTab(self.pv3d, "3-D surface")

        rv.addWidget(self.vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([380, 820])

    # ------------------------------------------------------------------
    def _populate_table(self):
        phases = self.mb.get_all_phases()
        self.table.setRowCount(len(phases))
        for row, p in enumerate(phases):
            self.table.setItem(row, 0, QTableWidgetItem(str(p["phase"])))
            self.table.setItem(row, 1, QTableWidgetItem(str(p["crystal_system"])))
            self.table.setItem(row, 2, QTableWidgetItem(str(p["primary_phase"])))
            try:
                rho_str = f"{float(p['density']):.3f}"
            except (TypeError, ValueError):
                rho_str = str(p["density"])
            self.table.setItem(row, 3, QTableWidgetItem(rho_str))

    def _selected_phase(self) -> str | None:
        rows = self.table.selectedItems()
        if not rows:
            return None
        return self.table.item(self.table.currentRow(), 0).text()

    def _on_select(self):
        phase = self._selected_phase()
        if phase is None:
            return
        rho = self.mb.get_density(phase)
        self.info_label.setText(
            f"<b>{phase}</b>  |  ρ = {rho:.1f} kg/m³" if rho else f"<b>{phase}</b>"
        )

    # ------------------------------------------------------------------
    # Compute — load tensor → metrics → kick off async surface computation
    # ------------------------------------------------------------------

    def _compute(self):
        phase = self._selected_phase()
        if phase is None:
            return

        p_gpa = self.pressure_spin.value()
        t_k   = self.temp_spin.value()

        cij = self.mb.get_voigt_matrix_gpa(phase, p_gpa, t_k)
        rho = self.mb.get_density(phase)

        if cij is None or rho is None:
            self.info_label.setText(
                "Could not retrieve elastic constants for this phase / P-T range."
            )
            return

        self._fill_voigt_table(cij)

        # Load into the backend (Pa)
        self.ab.set_from_voigt(cij * 1e9, rho)

        # Scalar metrics (fast — synchronous)
        anis = self.ab.compute_anisotropy_values()
        if anis:
            self._fill_anis_table(anis)

        # Velocity surface (potentially slow — async worker)
        self._surf_data = None
        self.calc_btn.setEnabled(False)
        self.progress.setVisible(True)

        self._worker = _SurfaceWorker(self.ab, n_theta=60, n_phi=120)
        self._worker.finished.connect(self._on_surface_done)
        self._worker.error.connect(self._on_surface_error)
        self._worker.start()

    def _on_surface_done(self, surf: dict):
        self._surf_data = surf
        self.progress.setVisible(False)
        self.calc_btn.setEnabled(True)
        self._replot()

    def _on_surface_error(self, msg: str):
        self.progress.setVisible(False)
        self.calc_btn.setEnabled(True)
        self.info_label.setText(f"Surface error: {msg}")
        # Still try to draw a basic stereonet from the tensor already loaded
        self._replot_stereonet_only()

    # ------------------------------------------------------------------
    # Replot helpers
    # ------------------------------------------------------------------

    def _replot_if_ready(self):
        """Called when display options change — only redraws if data exists."""
        if self._surf_data is not None:
            self._replot()
        elif self.ab.cijkl is not None:
            # Surface not yet computed, but tensor is loaded → at least stereonet
            self._replot_stereonet_only()

    def _replot(self):
        """Full redraw: 2-D stereonet + 3-D surface."""
        self._replot_stereonet_only()
        scalar = self._sn_opts.scalar
        if self._surf_data is not None:
            self.pv3d.plot_velocity_surface(self._surf_data, scalar=scalar)

    def _replot_stereonet_only(self):
        opts   = self._sn_opts
        scalar = opts.scalar

        if opts.style == "Scatter (dots)":
            plot_data = self.ab.compute_stereonet_data(step_deg=2.0)
        else:
            plot_data = self.ab.compute_stereonet_grid(grid_size=300)

        if plot_data is None:
            return

        label = _SCALAR_LABELS.get(scalar, scalar.upper())
        fig = make_stereonet_figure(
            data=plot_data,
            scalar=scalar,
            style=opts.style,
            colorscale=opts.colorscale,
            vmin=opts.vmin,
            vmax=opts.vmax,
            show_colorbar=opts.show_colorbar,
            pt_size=opts.pt_size,
            title=f"{label} — upper-hemisphere stereonet",
        )
        self.mpl.show_figure(fig)
        self.vis_tabs.setCurrentIndex(0)   # switch to stereonet tab

    # ------------------------------------------------------------------
    # Table helpers
    # ------------------------------------------------------------------

    def _fill_voigt_table(self, cij: np.ndarray):
        for i in range(6):
            for j in range(6):
                self.voigt_table.setItem(
                    i, j, QTableWidgetItem(f"{cij[i, j]:.2f}")
                )

    def _fill_anis_table(self, anis: dict):
        items = [
            (k, v) for k, v in anis.items()
            if not isinstance(v, np.ndarray)
        ]
        self.anis_table.setRowCount(len(items))
        for row, (k, v) in enumerate(items):
            self.anis_table.setItem(row, 0, QTableWidgetItem(str(k)))
            val_str = f"{v:.4f}" if isinstance(v, float) else str(v)
            self.anis_table.setItem(row, 1, QTableWidgetItem(val_str))
