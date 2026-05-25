import numpy as np
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QTableWidget, QTableWidgetItem, QAbstractItemView,
    QDoubleSpinBox, QPushButton, QSplitter, QHeaderView,
    QComboBox, QFormLayout, QSizePolicy
)
from PyQt5.QtCore import Qt

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES
from frontend.tabs.ebsd_tab import _PlotOptions
import plotly.graph_objects as go


class MaterialTab(QWidget):
    """Browse the mineral database, pick P/T conditions, view elastic tensor."""

    def __init__(self, material_backend, anisotropy_backend, parent=None):
        super().__init__(parent)
        self.mb = material_backend
        self.ab = anisotropy_backend
        self._build_ui()
        self._populate_table()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ---- Left: material list + PT controls ----
        left = QWidget()
        left_layout = QVBoxLayout(left)

        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Phase", "Crystal System", "Primary Phase", "ρ (g/cm³)"])
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.selectionModel().selectionChanged.connect(self._on_select)
        left_layout.addWidget(QLabel("Mineral database:"))
        left_layout.addWidget(self.table)

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

        left_layout.addWidget(pt_box)

        self.calc_btn = QPushButton("Compute anisotropy")
        self.calc_btn.clicked.connect(self._compute)
        left_layout.addWidget(self.calc_btn)

        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        left_layout.addWidget(self.info_label)

        splitter.addWidget(left)

        # ---- Right: Voigt matrix + anisotropy metrics ----
        right = QWidget()
        right_layout = QVBoxLayout(right)

        right_layout.addWidget(QLabel("Voigt stiffness matrix (GPa):"))
        self.voigt_table = QTableWidget(6, 6)
        self.voigt_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.voigt_table.setMaximumHeight(220)
        self.voigt_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.voigt_table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        right_layout.addWidget(self.voigt_table)

        right_layout.addWidget(QLabel("Anisotropy values:"))
        self.anis_table = QTableWidget()
        self.anis_table.setColumnCount(2)
        self.anis_table.setHorizontalHeaderLabels(["Quantity", "Value"])
        self.anis_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.anis_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        right_layout.addWidget(self.anis_table)

        self._sn_opts = _PlotOptions("Stereonet options", default_cmap="RdBu_r",
                                     default_pt_size=2)
        right_layout.addWidget(self._sn_opts)
        self.mpl = PlotlyWidget()
        right_layout.addWidget(self.mpl)

        splitter.addWidget(right)
        splitter.setSizes([400, 700])

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

    def _compute(self):
        phase = self._selected_phase()
        if phase is None:
            return
        p_gpa = self.pressure_spin.value()
        t_k   = self.temp_spin.value()

        cij = self.mb.get_voigt_matrix_gpa(phase, p_gpa, t_k)
        rho = self.mb.get_density(phase)

        if cij is None or rho is None:
            self.info_label.setText("Could not retrieve elastic constants.")
            return

        self._fill_voigt_table(cij)

        # Convert GPa → Pa for anisotropy backend
        cij_pa = cij * 1e9
        self.ab.set_from_voigt(cij_pa, rho)
        anis = self.ab.compute_anisotropy_values()
        if anis:
            self._fill_anis_table(anis)
            self._plot_stereonet(anis)

    def _fill_voigt_table(self, cij: np.ndarray):
        for i in range(6):
            for j in range(6):
                self.voigt_table.setItem(i, j, QTableWidgetItem(f"{cij[i,j]:.2f}"))

    def _fill_anis_table(self, anis: dict):
        items = [(k, v) for k, v in anis.items() if not isinstance(v, np.ndarray)]
        self.anis_table.setRowCount(len(items))
        for row, (k, v) in enumerate(items):
            self.anis_table.setItem(row, 0, QTableWidgetItem(str(k)))
            self.anis_table.setItem(row, 1, QTableWidgetItem(f"{v:.4f}" if isinstance(v, float) else str(v)))

    def _plot_stereonet(self, anis: dict):
        import numpy as _np
        data = self.ab.compute_stereonet_data(step_deg=3.0)
        if data is None:
            return
        opts = self._sn_opts
        vals = data["vp"] / 1000.0
        vmin = opts.vmin if opts.vmin is not None else float(_np.nanmin(vals))
        vmax = opts.vmax if opts.vmax is not None else float(_np.nanmax(vals))
        cb   = dict(title="Vp (km/s)") if opts.show_colorbar else None
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
            hovertemplate="Vp=%{marker.color:.3f} km/s<extra></extra>",
        ))
        fig.update_layout(
            title="Vp stereonet (upper hemisphere)",
            xaxis=dict(range=[-1.1, 1.1], showgrid=False, zeroline=False,
                       showticklabels=False),
            yaxis=dict(range=[-1.1, 1.1], scaleanchor="x", showgrid=False,
                       zeroline=False, showticklabels=False),
        )
        self.mpl.show_figure(fig)
