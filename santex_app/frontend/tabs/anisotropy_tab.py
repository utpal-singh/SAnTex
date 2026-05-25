"""
AnisotropyTab — 6×6 Voigt matrix input + 2-D stereonet + 3-D surface.

Stereonet rendered with Plotly (QWebEngineView).
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QDoubleSpinBox, QSpinBox, QPushButton, QSplitter,
    QComboBox, QFormLayout, QTabWidget, QSizePolicy,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QProgressBar, QTextEdit,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES
from frontend.widgets.pyvista_widget import PyVistaWidget
from frontend.tabs.ebsd_tab import _PlotOptions

import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Worker thread
# ---------------------------------------------------------------------------

class _ComputeWorker(QThread):
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, backend, n_theta, n_phi):
        super().__init__()
        self.backend = backend
        self.n_theta = n_theta
        self.n_phi   = n_phi

    def run(self):
        try:
            surf = self.backend.compute_velocity_surface(self.n_theta, self.n_phi)
            if surf is None:
                self.error.emit("No elastic tensor loaded.")
                return
            self.finished.emit(surf)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Main tab
# ---------------------------------------------------------------------------

class AnisotropyTab(QWidget):
    def __init__(self, anisotropy_backend, parent=None):
        super().__init__(parent)
        self.ab = anisotropy_backend
        self._worker = None
        self._surf_data = None
        self._build_ui()

    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ---- Left: tensor input + controls ----
        left = QWidget()
        lv = QVBoxLayout(left)

        lv.addWidget(QLabel("6×6 Voigt stiffness matrix (GPa):"))
        self.voigt_input = QTableWidget(6, 6)
        self.voigt_input.setMaximumHeight(240)
        self.voigt_input.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.voigt_input.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._clear_voigt_input()
        lv.addWidget(self.voigt_input)

        density_row = QHBoxLayout()
        density_row.addWidget(QLabel("Density (kg/m³):"))
        self.density_spin = QDoubleSpinBox()
        self.density_spin.setRange(100, 20000)
        self.density_spin.setValue(3300)
        self.density_spin.setDecimals(1)
        density_row.addWidget(self.density_spin)
        lv.addLayout(density_row)

        res_box = QGroupBox("Velocity surface resolution")
        res_form = QFormLayout(res_box)
        self.n_theta_spin = QSpinBox()
        self.n_theta_spin.setRange(10, 300)
        self.n_theta_spin.setValue(60)
        res_form.addRow("n θ:", self.n_theta_spin)
        self.n_phi_spin = QSpinBox()
        self.n_phi_spin.setRange(10, 600)
        self.n_phi_spin.setValue(120)
        res_form.addRow("n φ:", self.n_phi_spin)
        lv.addWidget(res_box)

        self.compute_btn = QPushButton("Compute anisotropy")
        self.compute_btn.clicked.connect(self._start_compute)
        lv.addWidget(self.compute_btn)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        lv.addWidget(self.progress)

        lv.addWidget(QLabel("Scalar metrics:"))
        self.metrics_text = QTextEdit()
        self.metrics_text.setReadOnly(True)
        self.metrics_text.setMaximumHeight(180)
        lv.addWidget(self.metrics_text)

        lv.addStretch()
        splitter.addWidget(left)

        # ---- Right: visualisation ----
        right = QWidget()
        rv = QVBoxLayout(right)

        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel("Display:"))
        self.scalar_combo = QComboBox()
        self.scalar_combo.addItems(["vp", "vs1", "vs2", "avs", "vpvs1", "vpvs2"])
        self.scalar_combo.currentTextChanged.connect(self._replot)
        ctrl_row.addWidget(self.scalar_combo)
        ctrl_row.addStretch()
        rv.addLayout(ctrl_row)

        self._sn_opts = _PlotOptions("Stereonet options", default_cmap="RdBu_r",
                                     default_pt_size=3)
        rv.addWidget(self._sn_opts)

        self.vis_tabs = QTabWidget()
        self.mpl_2d = PlotlyWidget()
        self.vis_tabs.addTab(self.mpl_2d, "2-D stereonet")
        self.pv3d = PyVistaWidget()
        self.vis_tabs.addTab(self.pv3d, "3-D surface")
        rv.addWidget(self.vis_tabs)

        splitter.addWidget(right)
        splitter.setSizes([380, 720])

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_tensor(self, cij_gpa: np.ndarray, density: float):
        for i in range(6):
            for j in range(6):
                self.voigt_input.setItem(i, j, QTableWidgetItem(f"{cij_gpa[i,j]:.2f}"))
        self.density_spin.setValue(density)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _clear_voigt_input(self):
        for i in range(6):
            for j in range(6):
                self.voigt_input.setItem(i, j, QTableWidgetItem("0.0"))

    def _read_voigt_input(self) -> np.ndarray | None:
        try:
            cij = np.zeros((6, 6))
            for i in range(6):
                for j in range(6):
                    item = self.voigt_input.item(i, j)
                    cij[i, j] = float(item.text()) if item else 0.0
            return cij
        except ValueError:
            return None

    def _start_compute(self):
        cij_gpa = self._read_voigt_input()
        if cij_gpa is None:
            return
        self.ab.set_from_voigt(cij_gpa * 1e9, self.density_spin.value())
        self.compute_btn.setEnabled(False)
        self.progress.setVisible(True)
        self._worker = _ComputeWorker(self.ab, self.n_theta_spin.value(),
                                      self.n_phi_spin.value())
        self._worker.finished.connect(self._on_compute_done)
        self._worker.error.connect(self._on_compute_error)
        self._worker.start()

    def _on_compute_done(self, surf: dict):
        self._surf_data = surf
        self.progress.setVisible(False)
        self.compute_btn.setEnabled(True)
        self._show_metrics()
        self._replot()

    def _on_compute_error(self, msg: str):
        self.progress.setVisible(False)
        self.compute_btn.setEnabled(True)
        self.metrics_text.setPlainText(f"Error: {msg}")

    def _show_metrics(self):
        anis = self.ab.compute_anisotropy_values()
        if not anis:
            return
        lines = [f"{k}: {v:.4f}" if isinstance(v, float) else f"{k}: {v}"
                 for k, v in anis.items()]
        self.metrics_text.setPlainText("\n".join(lines))

    def _replot(self):
        if self._surf_data is None:
            return
        scalar = self.scalar_combo.currentText()
        self._plot_2d(scalar)
        self._plot_3d(scalar)

    def _plot_2d(self, scalar: str):
        data = self.ab.compute_stereonet_data(step_deg=3.0)
        if data is None:
            return
        opts  = self._sn_opts
        vals  = data[scalar]
        unit  = "km/s" if scalar in ("vp", "vs1", "vs2") else "%"
        if unit == "km/s":
            vals = vals / 1000.0

        vmin = opts.vmin if opts.vmin is not None else float(np.nanmin(vals))
        vmax = opts.vmax if opts.vmax is not None else float(np.nanmax(vals))
        cb   = dict(title=f"{scalar} ({unit})") if opts.show_colorbar else None

        fig = go.Figure()
        fig.add_shape(type="circle", xref="x", yref="y",
                      x0=-1, y0=-1, x1=1, y1=1,
                      line=dict(color="black", width=1))
        fig.add_trace(go.Scattergl(
            x=data["x"], y=data["y"],
            mode="markers",
            marker=dict(color=vals, colorscale=opts.colorscale,
                        cmin=vmin, cmax=vmax, size=opts.pt_size,
                        colorbar=cb, showscale=opts.show_colorbar),
            hovertemplate=f"{scalar}=%{{marker.color:.3f}} {unit}<extra></extra>",
            name=scalar,
        ))
        fig.update_layout(
            title=f"{scalar} — upper-hemisphere stereonet",
            xaxis=dict(range=[-1.1, 1.1], showgrid=False, zeroline=False,
                       showticklabels=False),
            yaxis=dict(range=[-1.1, 1.1], scaleanchor="x", showgrid=False,
                       zeroline=False, showticklabels=False),
        )
        self.mpl_2d.show_figure(fig)

    def _plot_3d(self, scalar: str):
        if self._surf_data:
            self.pv3d.plot_velocity_surface(self._surf_data, scalar=scalar)
