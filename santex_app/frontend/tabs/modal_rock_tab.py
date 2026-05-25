import numpy as np
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QDoubleSpinBox, QPushButton, QSplitter, QFormLayout,
    QComboBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QSpinBox, QSizePolicy, QTabWidget,
    QMessageBox, QTextEdit
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from frontend.widgets.plotly_widget import PlotlyWidget
from frontend.tabs.ebsd_tab import _PlotOptions
import plotly.graph_objects as go
from frontend.widgets.pyvista_widget import PyVistaWidget


# ---------------------------------------------------------------------------
# Worker threads
# ---------------------------------------------------------------------------

class _ModalWorker(QThread):
    finished = pyqtSignal(object, float)
    error    = pyqtSignal(str)

    def __init__(self, mb, minerals, fractions, pressure, temp):
        super().__init__()
        self.mb = mb
        self.minerals = minerals
        self.fractions = fractions
        self.pressure = pressure
        self.temp = temp

    def run(self):
        try:
            cij, rho = self.mb.modal_rock(
                self.minerals, self.fractions, self.pressure, self.temp)
            if cij is None:
                self.error.emit("modal_rock returned None.")
                return
            self.finished.emit(cij, rho)
        except Exception as e:
            self.error.emit(str(e))


class _HSWorker(QThread):
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, mb, phase_ids, fractions, temp_c, pressure):
        super().__init__()
        self.mb = mb
        self.phase_ids = phase_ids
        self.fractions = fractions
        self.temp_c = temp_c
        self.pressure = pressure

    def run(self):
        try:
            result = self.mb.hashin_shtrikman(
                self.phase_ids, self.fractions, self.temp_c, self.pressure)
            if result is None:
                self.error.emit("Hashin-Shtrikman returned None.")
                return
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Modal Rock tab
# ---------------------------------------------------------------------------

class ModalRockTab(QWidget):
    """
    Two sub-tabs:
      1. Anisotropic modal rock (Voigt average of several anisotropic minerals)
      2. Isotropic Hashin-Shtrikman bounds (uses Isotropy / Birch-Murnaghan EOS)
    """

    def __init__(self, material_backend, anisotropy_backend, parent=None):
        super().__init__(parent)
        self.mb = material_backend
        self.ab = anisotropy_backend
        self._worker = None
        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QVBoxLayout(self)
        self.sub_tabs = QTabWidget()
        root.addWidget(self.sub_tabs)

        self.sub_tabs.addTab(self._build_aniso_panel(), "Anisotropic (Voigt modal)")
        self.sub_tabs.addTab(self._build_hs_panel(),    "Isotropic (Hashin-Shtrikman)")

    # ------------------------------------------------------------------
    # Sub-panel 1: anisotropic modal
    # ------------------------------------------------------------------

    def _build_aniso_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        # Left: mineral rows
        left = QWidget()
        left_layout = QVBoxLayout(left)

        left_layout.addWidget(QLabel("Mineral assemblage:"))
        self.aniso_table = QTableWidget(0, 3)
        self.aniso_table.setHorizontalHeaderLabels(["Mineral", "Volume fraction", ""])
        self.aniso_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.aniso_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.aniso_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        left_layout.addWidget(self.aniso_table)

        add_btn = QPushButton("+ Add mineral")
        add_btn.clicked.connect(self._add_aniso_row)
        left_layout.addWidget(add_btn)

        pt_box = QGroupBox("Conditions")
        pt_form = QFormLayout(pt_box)
        self.aniso_pressure = QDoubleSpinBox()
        self.aniso_pressure.setRange(0, 300)
        self.aniso_pressure.setSuffix(" GPa")
        self.aniso_pressure.setDecimals(2)
        pt_form.addRow("Pressure:", self.aniso_pressure)

        self.aniso_temp = QDoubleSpinBox()
        self.aniso_temp.setRange(0, 3000)
        self.aniso_temp.setSuffix(" K")
        self.aniso_temp.setValue(300)
        self.aniso_temp.setDecimals(1)
        pt_form.addRow("Temperature:", self.aniso_temp)
        left_layout.addWidget(pt_box)

        self.aniso_compute_btn = QPushButton("Compute modal rock")
        self.aniso_compute_btn.clicked.connect(self._start_modal)
        left_layout.addWidget(self.aniso_compute_btn)

        self.aniso_result_text = QTextEdit()
        self.aniso_result_text.setReadOnly(True)
        self.aniso_result_text.setMaximumHeight(120)
        left_layout.addWidget(self.aniso_result_text)
        left_layout.addStretch()
        splitter.addWidget(left)

        # Right: stereonet options + plot + 3-D
        right = QWidget()
        right_layout = QVBoxLayout(right)
        self._aniso_sn_opts = _PlotOptions("Stereonet options", default_cmap="RdBu_r",
                                            default_pt_size=2)
        right_layout.addWidget(self._aniso_sn_opts)
        self.aniso_vis_tabs = QTabWidget()
        self.aniso_mpl  = PlotlyWidget()
        self.aniso_pv3d = PyVistaWidget()
        self.aniso_vis_tabs.addTab(self.aniso_mpl, "Stereonet")
        self.aniso_vis_tabs.addTab(self.aniso_pv3d, "3-D surface")
        right_layout.addWidget(self.aniso_vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([380, 720])
        return panel

    def _add_aniso_row(self):
        names = self.mb.get_phase_names()
        row = self.aniso_table.rowCount()
        self.aniso_table.insertRow(row)

        combo = QComboBox()
        combo.addItems(names)
        self.aniso_table.setCellWidget(row, 0, combo)

        frac_spin = QDoubleSpinBox()
        frac_spin.setRange(0, 1)
        frac_spin.setDecimals(3)
        frac_spin.setValue(1.0 / max(1, row + 1))
        frac_spin.setSingleStep(0.05)
        self.aniso_table.setCellWidget(row, 1, frac_spin)

        del_btn = QPushButton("✕")
        del_btn.clicked.connect(lambda _, r=row: self.aniso_table.removeRow(r))
        self.aniso_table.setCellWidget(row, 2, del_btn)

    def _start_modal(self):
        minerals, fractions = [], []
        for row in range(self.aniso_table.rowCount()):
            combo = self.aniso_table.cellWidget(row, 0)
            spin  = self.aniso_table.cellWidget(row, 1)
            if combo and spin:
                minerals.append(combo.currentText())
                fractions.append(spin.value())

        if not minerals:
            return

        total = sum(fractions)
        if abs(total - 1.0) > 0.01:
            QMessageBox.warning(self, "Fractions", f"Volume fractions sum to {total:.3f}, not 1.0.")
            return

        self.aniso_compute_btn.setEnabled(False)
        self._worker = _ModalWorker(
            self.mb, minerals, fractions,
            self.aniso_pressure.value(), self.aniso_temp.value()
        )
        self._worker.finished.connect(self._on_modal_done)
        self._worker.error.connect(self._on_modal_error)
        self._worker.start()

    def _on_modal_done(self, cij_gpa: np.ndarray, rho: float):
        self.aniso_compute_btn.setEnabled(True)
        self.ab.set_from_voigt(cij_gpa * 1e9, rho)
        anis = self.ab.compute_anisotropy_values()
        if anis:
            lines = [f"{k}: {v:.4f}" for k, v in anis.items() if isinstance(v, float)]
            self.aniso_result_text.setPlainText(f"ρ = {rho:.1f} kg/m³\n" + "\n".join(lines))

        data = self.ab.compute_stereonet_data()
        if data:
            import numpy as _np
            opts = self._aniso_sn_opts
            vals = data["vp"] / 1000.0
            vmin = opts.vmin if opts.vmin is not None else float(_np.nanmin(vals))
            vmax = opts.vmax if opts.vmax is not None else float(_np.nanmax(vals))
            cb   = dict(title="Vp (km/s)") if opts.show_colorbar else None
            fig  = go.Figure()
            fig.add_shape(type="circle", xref="x", yref="y",
                          x0=-1, y0=-1, x1=1, y1=1,
                          line=dict(color="black", width=1))
            fig.add_trace(go.Scattergl(
                x=data["x"], y=data["y"],
                mode="markers",
                marker=dict(color=vals, colorscale=opts.colorscale,
                            cmin=vmin, cmax=vmax, size=opts.pt_size,
                            colorbar=cb, showscale=opts.show_colorbar),
            ))
            fig.update_layout(
                title="Vp — modal rock",
                xaxis=dict(range=[-1.1, 1.1], showgrid=False, zeroline=False,
                           showticklabels=False),
                yaxis=dict(range=[-1.1, 1.1], scaleanchor="x", showgrid=False,
                           zeroline=False, showticklabels=False),
            )
            self.aniso_mpl.show_figure(fig)

        surf = self.ab.compute_velocity_surface()
        if surf:
            self.aniso_pv3d.plot_velocity_surface(surf, scalar="vp")

    def _on_modal_error(self, msg: str):
        self.aniso_compute_btn.setEnabled(True)
        QMessageBox.critical(self, "Error", msg)

    # ------------------------------------------------------------------
    # Sub-panel 2: Hashin-Shtrikman
    # ------------------------------------------------------------------

    def _build_hs_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        # Left
        left = QWidget()
        left_layout = QVBoxLayout(left)

        left_layout.addWidget(QLabel("Phase assemblage (isotropic EOS):"))
        self.hs_table = QTableWidget(0, 3)
        self.hs_table.setHorizontalHeaderLabels(["Phase", "Volume fraction", ""])
        self.hs_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.hs_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.hs_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        left_layout.addWidget(self.hs_table)

        add_hs_btn = QPushButton("+ Add phase")
        add_hs_btn.clicked.connect(self._add_hs_row)
        left_layout.addWidget(add_hs_btn)

        hs_pt_box = QGroupBox("Conditions")
        hs_pt_form = QFormLayout(hs_pt_box)
        self.hs_pressure = QDoubleSpinBox()
        self.hs_pressure.setRange(0, 300)
        self.hs_pressure.setSuffix(" GPa")
        self.hs_pressure.setDecimals(2)
        hs_pt_form.addRow("Pressure:", self.hs_pressure)

        self.hs_temp = QDoubleSpinBox()
        self.hs_temp.setRange(-273, 3000)
        self.hs_temp.setSuffix(" °C")
        self.hs_temp.setValue(25)
        self.hs_temp.setDecimals(1)
        hs_pt_form.addRow("Temperature:", self.hs_temp)
        left_layout.addWidget(hs_pt_box)

        self.hs_compute_btn = QPushButton("Compute HS bounds")
        self.hs_compute_btn.clicked.connect(self._start_hs)
        left_layout.addWidget(self.hs_compute_btn)

        self.hs_result_text = QTextEdit()
        self.hs_result_text.setReadOnly(True)
        left_layout.addWidget(self.hs_result_text)
        left_layout.addStretch()
        splitter.addWidget(left)

        # Right: plot
        right = QWidget()
        right_layout = QVBoxLayout(right)
        self.hs_mpl = PlotlyWidget()
        right_layout.addWidget(self.hs_mpl)
        splitter.addWidget(right)
        splitter.setSizes([380, 720])
        return panel

    def _add_hs_row(self):
        iso_phases = self.mb.get_isotropic_phases()   # [(phase_id, name), ...]
        row = self.hs_table.rowCount()
        self.hs_table.insertRow(row)

        combo = QComboBox()
        for pid, pname in iso_phases:
            combo.addItem(pname, userData=pid)
        self.hs_table.setCellWidget(row, 0, combo)

        frac_spin = QDoubleSpinBox()
        frac_spin.setRange(0, 1)
        frac_spin.setDecimals(3)
        frac_spin.setValue(1.0 / max(1, row + 1))
        frac_spin.setSingleStep(0.05)
        self.hs_table.setCellWidget(row, 1, frac_spin)

        del_btn = QPushButton("✕")
        del_btn.clicked.connect(lambda _, r=row: self.hs_table.removeRow(r))
        self.hs_table.setCellWidget(row, 2, del_btn)

    def _start_hs(self):
        phase_ids, fractions = [], []
        for row in range(self.hs_table.rowCount()):
            combo = self.hs_table.cellWidget(row, 0)
            spin  = self.hs_table.cellWidget(row, 1)
            if combo and spin:
                phase_ids.append(combo.currentData())
                fractions.append(spin.value())

        if not phase_ids:
            return

        self.hs_compute_btn.setEnabled(False)
        self._hs_worker = _HSWorker(
            self.mb, phase_ids, fractions,
            self.hs_temp.value(), self.hs_pressure.value()
        )
        self._hs_worker.finished.connect(self._on_hs_done)
        self._hs_worker.error.connect(self._on_hs_error)
        self._hs_worker.start()

    def _on_hs_done(self, result: dict):
        self.hs_compute_btn.setEnabled(True)
        medium = result.get("medium", {})
        upper  = result.get("upper",  {})
        lower  = result.get("lower",  {})

        lines = ["=== Hashin-Shtrikman bounds ===\n"]
        for label, d in [("VRH medium", medium), ("HS upper", upper), ("HS lower", lower)]:
            lines.append(f"--- {label} ---")
            for k, v in d.items():
                if isinstance(v, (int, float)):
                    lines.append(f"  {k}: {v:.4f}")
            lines.append("")
        self.hs_result_text.setPlainText("\n".join(lines))

        # Grouped bar chart: Vp / Vs for lower / VRH / upper
        labels, vp_vals, vs_vals = [], [], []
        for lbl, d in [("Lower", lower), ("VRH", medium), ("Upper", upper)]:
            vp = d.get("Vp", d.get("vp", None))
            vs = d.get("Vs", d.get("vs", None))
            if vp is not None:
                labels.append(lbl)
                vp_vals.append(vp)
                if vs is not None:
                    vs_vals.append(vs)
        fig = go.Figure()
        if vp_vals:
            fig.add_trace(go.Bar(name="Vp", x=labels, y=vp_vals,
                                 marker_color="#4e79a7"))
        if vs_vals and len(vs_vals) == len(labels):
            fig.add_trace(go.Bar(name="Vs", x=labels, y=vs_vals,
                                 marker_color="#f28e2b"))
        fig.update_layout(barmode="group",
                          title="Hashin-Shtrikman bounds",
                          yaxis_title="Velocity (km/s)")
        self.hs_mpl.show_figure(fig)

    def _on_hs_error(self, msg: str):
        self.hs_compute_btn.setEnabled(True)
        QMessageBox.critical(self, "Error", msg)
