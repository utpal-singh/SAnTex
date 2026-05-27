"""
CrystalTab — interactive 3-D crystal geometry viewer (PyVista).

Shows the unit cell as a wireframe parallelepiped, optional Miller-plane
polygons coloured by index family, and direction arrows with [uvw] labels.

Minerals included:
  Forsterite  Mg₂SiO₄  Orthorhombic  Pbnm (62)
  Enstatite   MgSiO₃   Orthorhombic  Pbca (61)
  Diopside    CaMgSi₂O₆ Monoclinic   C2/c (15)
  Quartz      SiO₂     Trigonal      P3₁21 (152)
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QFormLayout,
    QGroupBox, QLabel, QComboBox, QPushButton, QCheckBox,
    QScrollArea, QSplitter, QLineEdit, QColorDialog,
    QDoubleSpinBox, QSpinBox, QHBoxLayout, QFrame, QSlider,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor

# ---------------------------------------------------------------------------
# Crystallographic data
# ---------------------------------------------------------------------------

# Colours for Miller plane families (cycle through these)
_PLANE_COLORS = [
    "#e74c3c", "#3498db", "#2ecc71", "#f39c12",
    "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
]

CRYSTALS: dict[str, dict] = {
    "Forsterite (Mg₂SiO₄)": {
        "formula": "Mg₂SiO₄",
        "system":  "Orthorhombic",
        "a": 4.756, "b": 10.195, "c": 5.981,
        "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
        "spacegroup": "Pbnm  (No. 62)",
        "color": "#3498DB",
        "common_planes": [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 1, 0), (1, 0, 1), (0, 1, 1),
            (1, 1, 1), (2, 1, 0),
        ],
        "common_dirs": [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [1, 1, 1],
        ],
    },
    "Enstatite (MgSiO₃)": {
        "formula": "MgSiO₃",
        "system":  "Orthorhombic",
        "a": 18.228, "b": 8.815, "c": 5.179,
        "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
        "spacegroup": "Pbca  (No. 61)",
        "color": "#E67E22",
        "common_planes": [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
        ],
        "common_dirs": [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 1, 1],
        ],
    },
    "Diopside (CaMgSi₂O₆)": {
        "formula": "CaMgSi₂O₆",
        "system":  "Monoclinic",
        "a": 9.746, "b": 8.899, "c": 5.251,
        "alpha": 90.0, "beta": 105.79, "gamma": 90.0,
        "spacegroup": "C2/c  (No. 15)",
        "color": "#2ECC71",
        "common_planes": [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 1, 0), (1, -1, 0), (1, 0, 1),
            (0, 1, 1), (1, 1, 1),
        ],
        "common_dirs": [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [1, 1, 1],
        ],
    },
    "Quartz (SiO₂)": {
        "formula": "SiO₂",
        "system":  "Trigonal",
        "a": 4.913, "b": 4.913, "c": 5.405,
        "alpha": 90.0, "beta": 90.0, "gamma": 120.0,
        "spacegroup": "P3₁21  (No. 152)",
        "color": "#9B59B6",
        "common_planes": [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 0, 1), (1, 1, 0), (1, -1, 0),
            (1, 0, -1), (1, 1, 1),
        ],
        "common_dirs": [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [1, 1, 1],
        ],
    },
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _lattice_matrix(a: float, b: float, c: float,
                    alpha: float, beta: float, gamma: float) -> np.ndarray:
    """3×3 column matrix  M  such that  M @ [u,v,w]  gives Cartesian coords.

    Convention: a along X, b in XY plane.
    """
    al = np.radians(alpha)
    be = np.radians(beta)
    ga = np.radians(gamma)

    cos_a, cos_b, cos_g = np.cos(al), np.cos(be), np.cos(ga)
    sin_g = np.sin(ga)

    ax = a
    bx = b * cos_g
    by = b * sin_g
    cx = c * cos_b
    cy = c * (cos_a - cos_b * cos_g) / sin_g
    v  = 1.0 - cos_a**2 - cos_b**2 - cos_g**2 + 2.0 * cos_a * cos_b * cos_g
    cz = c * np.sqrt(max(0.0, v)) / sin_g

    return np.array([[ax, bx, cx],
                     [0., by, cy],
                     [0., 0., cz]])


def _cell_vertices(M: np.ndarray) -> np.ndarray:
    """8 Cartesian vertices of the unit cell parallelepiped."""
    frac = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
        [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
    ], dtype=float)
    return (M @ frac.T).T          # shape (8, 3)


_CELL_EDGES = [
    (0, 1), (2, 3), (4, 5), (6, 7),   # along a
    (0, 2), (1, 3), (4, 6), (5, 7),   # along b
    (0, 4), (1, 5), (2, 6), (3, 7),   # along c
]

_CELL_FACES_FRAC = [
    # (u,v,w) fractional-coord quads for each face: ±a, ±b, ±c
    [(0,0,0),(0,1,0),(0,1,1),(0,0,1)],   # u=0  ≡ (100)
    [(1,0,0),(1,1,0),(1,1,1),(1,0,1)],   # u=1
    [(0,0,0),(1,0,0),(1,0,1),(0,0,1)],   # v=0  ≡ (010)
    [(0,1,0),(1,1,0),(1,1,1),(0,1,1)],   # v=1
    [(0,0,0),(1,0,0),(1,1,0),(0,1,0)],   # w=0  ≡ (001)
    [(0,0,1),(1,0,1),(1,1,1),(0,1,1)],   # w=1
]


def _miller_plane_polygon(h: int, k: int, l: int,
                          M: np.ndarray) -> "np.ndarray | None":
    """Find the ordered polygon formed by plane (hkl) clipped to the unit cell.

    Returns an (N, 3) float array of Cartesian vertices, or None if the plane
    misses the cell interior.

    The intercept form is:  h·u + k·v + l·w = 1
    (translated along the normal so it passes through the cell if possible).
    """
    # Edge endpoints in fractional coordinates
    edges_frac = [
        ([0,0,0],[1,0,0]), ([0,1,0],[1,1,0]),
        ([0,0,1],[1,0,1]), ([0,1,1],[1,1,1]),
        ([0,0,0],[0,1,0]), ([1,0,0],[1,1,0]),
        ([0,0,1],[0,1,1]), ([1,0,1],[1,1,1]),
        ([0,0,0],[0,0,1]), ([1,0,0],[1,0,1]),
        ([0,1,0],[0,1,1]), ([1,1,0],[1,1,1]),
    ]
    hkl = np.array([h, k, l], dtype=float)

    # Try to find a good intercept value: n such that the plane cuts the cell
    # nmax = h+k+l if all positive; try n=1 first then n=h+k+l if needed
    def _intersect(n_val: float):
        pts = []
        for p0f, p1f in edges_frac:
            p0 = np.array(p0f, dtype=float)
            p1 = np.array(p1f, dtype=float)
            d  = p1 - p0
            den = float(hkl @ d)
            if abs(den) < 1e-10:
                continue
            t = (n_val - float(hkl @ p0)) / den
            if -1e-9 <= t <= 1 + 1e-9:
                t = float(np.clip(t, 0.0, 1.0))
                pf = p0 + t * d
                if all(-1e-9 <= pf[i] <= 1 + 1e-9 for i in range(3)):
                    pts.append(M @ pf)
        return pts

    # Determine a good n value (place the plane near the cell centre)
    pts = []
    for n_try in (1, 2, 0.5, 1.5, 0.25, 0.75):
        pts = _intersect(n_try)
        if len(pts) >= 3:
            break
    if len(pts) < 3:
        return None

    pts = np.array(pts)
    # Remove near-duplicates
    unique = [pts[0]]
    for p in pts[1:]:
        if all(np.linalg.norm(p - u) > 1e-6 for u in unique):
            unique.append(p)
    pts = np.array(unique)
    if len(pts) < 3:
        return None

    # Order points by angle around centroid (Graham scan in local 2D)
    centroid = pts.mean(axis=0)
    normal_cart = np.linalg.solve(M.T, hkl)          # M^{-T} @ hkl
    nlen = np.linalg.norm(normal_cart)
    if nlen < 1e-10:
        return None
    normal_cart /= nlen

    u_ax = pts[0] - centroid
    u_len = np.linalg.norm(u_ax)
    if u_len < 1e-10:
        return None
    u_ax /= u_len
    v_ax = np.cross(normal_cart, u_ax)

    angles = [np.arctan2(float((p - centroid) @ v_ax),
                          float((p - centroid) @ u_ax)) for p in pts]
    order = np.argsort(angles)
    return pts[order]


def _direction_arrow(u: int, v: int, w: int, M: np.ndarray,
                     scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (start, direction_unit) Cartesian vectors for [uvw]."""
    vec = M @ np.array([u, v, w], dtype=float)
    length = np.linalg.norm(vec)
    if length < 1e-10:
        return np.zeros(3), np.array([1, 0, 0])
    return np.zeros(3), vec / length * scale * length


def _plane_label(h: int, k: int, l: int) -> str:
    def _idx(i: int) -> str:
        if i < 0:
            return f"{abs(i)}̅"   # overline for negative
        return str(i)
    return f"({_idx(h)}{_idx(k)}{_idx(l)})"


def _dir_label(u: int, v: int, w: int) -> str:
    def _idx(i: int) -> str:
        if i < 0:
            return f"{abs(i)}̅"
        return str(i)
    return f"[{_idx(u)}{_idx(v)}{_idx(w)}]"


def _parse_hkl(text: str):
    """Parse 'h k l' or 'h,k,l' into (int,int,int). Raises ValueError."""
    import re
    nums = re.findall(r"-?\d+", text.strip())
    if len(nums) != 3:
        raise ValueError("Enter exactly three integers, e.g.  1 1 0")
    return tuple(int(n) for n in nums)


# ---------------------------------------------------------------------------
# PyVista render function
# ---------------------------------------------------------------------------

def _build_crystal_scene(
        plotter,
        crystal_name: str,
        show_faces: list[tuple[int, int, int]],
        show_dirs:  list[list[int]],
        show_cell:  bool = True,
        face_opacity: float = 0.45,
        show_labels: bool = True,
        bg_color: str = "#1e1e2e",
        cell_color: str = "#ffffff",
        axes_scale: float = 1.0,
):
    """Clear *plotter* and render the crystal unit cell scene."""
    try:
        import pyvista as pv
    except ImportError:
        return

    data = CRYSTALS[crystal_name]
    M = _lattice_matrix(
        data["a"], data["b"], data["c"],
        data["alpha"], data["beta"], data["gamma"],
    )
    pts = _cell_vertices(M)

    plotter.clear()
    plotter.set_background(bg_color)

    # ── Unit cell wireframe ──────────────────────────────────────────────
    if show_cell:
        for i0, i1 in _CELL_EDGES:
            line = pv.Line(pts[i0], pts[i1])
            plotter.add_mesh(line, color=cell_color, line_width=2,
                             render_lines_as_tubes=True)

        # Lattice vector labels at the cell corners
        if show_labels:
            a_end = M @ np.array([1, 0, 0])
            b_end = M @ np.array([0, 1, 0])
            c_end = M @ np.array([0, 0, 1])
            for vec, lbl, col in [
                (a_end, "a", "#e74c3c"),
                (b_end, "b", "#2ecc71"),
                (c_end, "c", "#3498db"),
            ]:
                plotter.add_point_labels(
                    [vec * 1.08], [lbl],
                    font_size=14, bold=True, text_color=col,
                    show_points=False, always_visible=True,
                )

    # ── Lattice-vector arrows ────────────────────────────────────────────
    origin = np.zeros(3)
    a_max = max(data["a"], data["b"], data["c"])
    arr_scale = a_max * 0.22

    for vec_frac, col in [
        ([1, 0, 0], "#e74c3c"),
        ([0, 1, 0], "#2ecc71"),
        ([0, 0, 1], "#3498db"),
    ]:
        v_cart = M @ np.array(vec_frac, dtype=float)
        v_len  = np.linalg.norm(v_cart)
        arrow  = pv.Arrow(
            start=origin,
            direction=v_cart / v_len,
            scale=v_len + arr_scale,
            shaft_radius=0.015,
            tip_radius=0.05,
            tip_length=0.12,
        )
        plotter.add_mesh(arrow, color=col)

    # ── Miller-plane polygons ────────────────────────────────────────────
    label_pts, label_texts = [], []
    for ci, (h, k, l) in enumerate(show_faces):
        poly_pts = _miller_plane_polygon(h, k, l, M)
        if poly_pts is None or len(poly_pts) < 3:
            continue
        n_face = len(poly_pts)
        faces  = np.zeros(n_face + 1, dtype=int)
        faces[0] = n_face
        faces[1:] = np.arange(n_face)
        mesh = pv.PolyData(poly_pts.astype(np.float32), faces=faces)
        col  = _PLANE_COLORS[ci % len(_PLANE_COLORS)]
        plotter.add_mesh(mesh, color=col, opacity=face_opacity,
                         show_edges=True, edge_color=col, line_width=1.5)
        if show_labels:
            centroid = poly_pts.mean(axis=0)
            # Push label slightly outside cell
            normal_dir = np.linalg.solve(M.T, np.array([h, k, l], float))
            n_len = np.linalg.norm(normal_dir)
            if n_len > 1e-10:
                centroid += normal_dir / n_len * a_max * 0.08
            label_pts.append(centroid)
            label_texts.append(_plane_label(h, k, l))

    if show_labels and label_pts:
        plotter.add_point_labels(
            label_pts, label_texts,
            font_size=12, bold=True, text_color="#f5f5f5",
            shape="rounded_rect", shape_color="#333333",
            shape_opacity=0.75, show_points=False, always_visible=True,
        )

    # ── Direction arrows ─────────────────────────────────────────────────
    dir_label_pts, dir_label_texts = [], []
    dir_colors = ["#f39c12", "#e74c3c", "#9b59b6", "#1abc9c",
                  "#34495e", "#e67e22", "#f1c40f", "#16a085"]
    for di, uvw in enumerate(show_dirs):
        u, v, w = uvw
        vec = M @ np.array([u, v, w], dtype=float)
        v_len = np.linalg.norm(vec)
        if v_len < 1e-10:
            continue
        col = dir_colors[di % len(dir_colors)]
        tip_start = origin + vec * 0.5    # arrow from centre of cell
        arrow = pv.Arrow(
            start=tip_start - vec * 0.5,
            direction=vec / v_len,
            scale=v_len * 1.15,
            shaft_radius=0.02,
            tip_radius=0.065,
            tip_length=0.15,
        )
        plotter.add_mesh(arrow, color=col)
        if show_labels:
            dir_label_pts.append(tip_start + vec / v_len * v_len * 0.65)
            dir_label_texts.append(_dir_label(u, v, w))

    if show_labels and dir_label_pts:
        plotter.add_point_labels(
            dir_label_pts, dir_label_texts,
            font_size=11, bold=True, text_color="#f39c12",
            shape="rounded_rect", shape_color="#222222",
            shape_opacity=0.7, show_points=False, always_visible=True,
        )

    # ── Axes widget ──────────────────────────────────────────────────────
    plotter.add_axes(
        xlabel="a", ylabel="b", zlabel="c",
        line_width=3,
    )
    plotter.reset_camera()


# ---------------------------------------------------------------------------
# Crystal Tab widget
# ---------------------------------------------------------------------------

class CrystalTab(QWidget):
    """Interactive 3-D crystal unit-cell viewer."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._face_checks: dict[tuple, QCheckBox] = {}
        self._dir_checks:  dict[tuple, QCheckBox] = {}
        self._bg_color    = "#1e1e2e"
        self._cell_color  = "#ffffff"
        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ── Left: scrollable controls ──────────────────────────────────
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(340)
        left_w = QWidget()
        lv = QVBoxLayout(left_w)
        lv.setContentsMargins(4, 4, 4, 4)
        scroll.setWidget(left_w)

        # Crystal selector
        sel_box = QGroupBox("Crystal")
        sel_form = QFormLayout(sel_box)
        self.crystal_combo = QComboBox()
        self.crystal_combo.addItems(list(CRYSTALS.keys()))
        self.crystal_combo.currentTextChanged.connect(self._on_crystal_changed)
        sel_form.addRow("Mineral:", self.crystal_combo)

        self.info_label = QLabel()
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("color:#555; font-size:10px;")
        sel_form.addRow(self.info_label)
        lv.addWidget(sel_box)

        # Lattice parameters (read-only display)
        self.params_box = QGroupBox("Lattice parameters")
        self.params_label = QLabel()
        self.params_label.setWordWrap(True)
        self.params_label.setStyleSheet("font-family:monospace; font-size:10px;")
        pbl = QVBoxLayout(self.params_box)
        pbl.addWidget(self.params_label)
        lv.addWidget(self.params_box)

        # Miller planes
        self.planes_box = QGroupBox("Miller planes  (hkl)")
        pv2 = QVBoxLayout(self.planes_box)
        self._planes_check_area = QWidget()
        self._planes_check_layout = QVBoxLayout(self._planes_check_area)
        self._planes_check_layout.setSpacing(2)
        pv2.addWidget(self._planes_check_area)

        # Custom plane input
        cp_row = QHBoxLayout()
        self._custom_plane_edit = QLineEdit()
        self._custom_plane_edit.setPlaceholderText("h k l  e.g. 2 1 1")
        add_plane_btn = QPushButton("Add")
        add_plane_btn.setFixedWidth(42)
        add_plane_btn.clicked.connect(self._add_custom_plane)
        cp_row.addWidget(self._custom_plane_edit)
        cp_row.addWidget(add_plane_btn)
        pv2.addLayout(cp_row)
        lv.addWidget(self.planes_box)

        # Miller directions
        self.dirs_box = QGroupBox("Directions  [uvw]")
        dv2 = QVBoxLayout(self.dirs_box)
        self._dirs_check_area = QWidget()
        self._dirs_check_layout = QVBoxLayout(self._dirs_check_area)
        self._dirs_check_layout.setSpacing(2)
        dv2.addWidget(self._dirs_check_area)

        cd_row = QHBoxLayout()
        self._custom_dir_edit = QLineEdit()
        self._custom_dir_edit.setPlaceholderText("u v w  e.g. 1 1 0")
        add_dir_btn = QPushButton("Add")
        add_dir_btn.setFixedWidth(42)
        add_dir_btn.clicked.connect(self._add_custom_dir)
        cd_row.addWidget(self._custom_dir_edit)
        cd_row.addWidget(add_dir_btn)
        dv2.addLayout(cd_row)
        lv.addWidget(self.dirs_box)

        # Display options
        disp_box = QGroupBox("Display options")
        disp_form = QFormLayout(disp_box)
        self.cell_check  = QCheckBox("Show cell wireframe")
        self.cell_check.setChecked(True)
        self.label_check = QCheckBox("Show labels")
        self.label_check.setChecked(True)

        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0.05, 1.0)
        self.opacity_spin.setSingleStep(0.05)
        self.opacity_spin.setValue(0.45)
        self.opacity_spin.setDecimals(2)

        self.bg_btn   = QPushButton()
        self.cell_btn = QPushButton()
        self._set_btn_color(self.bg_btn,   self._bg_color)
        self._set_btn_color(self.cell_btn, self._cell_color)
        self.bg_btn.setFixedHeight(22)
        self.cell_btn.setFixedHeight(22)
        self.bg_btn.clicked.connect(self._pick_bg_color)
        self.cell_btn.clicked.connect(self._pick_cell_color)

        disp_form.addRow(self.cell_check)
        disp_form.addRow(self.label_check)
        disp_form.addRow("Face opacity:", self.opacity_spin)
        disp_form.addRow("Background:", self.bg_btn)
        disp_form.addRow("Cell edges:", self.cell_btn)
        lv.addWidget(disp_box)

        # Render button
        render_btn = QPushButton("⚙  Render")
        render_btn.setMinimumHeight(34)
        render_btn.setStyleSheet("font-weight:bold; font-size:13px;")
        render_btn.clicked.connect(self._render)
        lv.addWidget(render_btn)
        lv.addStretch()

        splitter.addWidget(scroll)

        # ── Right: PyVista viewer ──────────────────────────────────────
        self._pv_widget = _make_pv_widget(self)
        splitter.addWidget(self._pv_widget)
        splitter.setSizes([340, 800])

        # Populate first crystal
        self._on_crystal_changed(self.crystal_combo.currentText())

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _set_btn_color(self, btn: QPushButton, hex_color: str):
        btn.setText(hex_color)
        btn.setStyleSheet(
            f"background-color:{hex_color};"
            f"color:{'#000' if _lightness(hex_color) > 0.5 else '#fff'};"
            f"border:1px solid #888;"
        )

    def _pick_bg_color(self):
        c = QColorDialog.getColor(QColor(self._bg_color), self, "Background colour")
        if c.isValid():
            self._bg_color = c.name()
            self._set_btn_color(self.bg_btn, self._bg_color)

    def _pick_cell_color(self):
        c = QColorDialog.getColor(QColor(self._cell_color), self, "Cell edge colour")
        if c.isValid():
            self._cell_color = c.name()
            self._set_btn_color(self.cell_btn, self._cell_color)

    # ------------------------------------------------------------------
    # Crystal changed
    # ------------------------------------------------------------------

    def _on_crystal_changed(self, name: str):
        data = CRYSTALS.get(name, {})
        self.info_label.setText(
            f"{data.get('formula','')}  ·  {data.get('system','')}  "
            f"·  {data.get('spacegroup','')}"
        )
        a, b, c = data.get("a", 0), data.get("b", 0), data.get("c", 0)
        al, be, ga = data.get("alpha", 90), data.get("beta", 90), data.get("gamma", 90)
        self.params_label.setText(
            f"a = {a:.3f} Å    b = {b:.3f} Å    c = {c:.3f} Å\n"
            f"α = {al:.2f}°  β = {be:.2f}°  γ = {ga:.2f}°"
        )
        self._rebuild_plane_checks(name)
        self._rebuild_dir_checks(name)
        self._render()

    # ------------------------------------------------------------------
    # Rebuild check-box panels
    # ------------------------------------------------------------------

    def _clear_layout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    def _rebuild_plane_checks(self, name: str):
        self._clear_layout(self._planes_check_layout)
        self._face_checks.clear()
        data = CRYSTALS.get(name, {})
        for i, hkl in enumerate(data.get("common_planes", [])):
            h, k, l = hkl
            lbl  = _plane_label(h, k, l)
            col  = _PLANE_COLORS[i % len(_PLANE_COLORS)]
            chk  = QCheckBox()
            chk.setChecked(i < 3)      # first three pre-ticked
            row_w = QWidget()
            row_l = QHBoxLayout(row_w)
            row_l.setContentsMargins(0, 0, 0, 0)
            swatch = QLabel("  ")
            swatch.setFixedWidth(14)
            swatch.setStyleSheet(f"background:{col}; border-radius:3px;")
            row_l.addWidget(chk)
            row_l.addWidget(swatch)
            row_l.addWidget(QLabel(lbl))
            row_l.addStretch()
            self._planes_check_layout.addWidget(row_w)
            self._face_checks[(h, k, l)] = chk

    def _rebuild_dir_checks(self, name: str):
        self._clear_layout(self._dirs_check_layout)
        self._dir_checks.clear()
        data = CRYSTALS.get(name, {})
        dir_colors = ["#f39c12", "#e74c3c", "#9b59b6",
                      "#1abc9c", "#e67e22", "#f1c40f"]
        for i, uvw in enumerate(data.get("common_dirs", [])):
            u, v, w = uvw
            lbl = _dir_label(u, v, w)
            col = dir_colors[i % len(dir_colors)]
            chk = QCheckBox()
            chk.setChecked(i < 3)
            row_w = QWidget()
            row_l = QHBoxLayout(row_w)
            row_l.setContentsMargins(0, 0, 0, 0)
            swatch = QLabel("→")
            swatch.setStyleSheet(f"color:{col}; font-weight:bold; font-size:14px;")
            swatch.setFixedWidth(18)
            row_l.addWidget(chk)
            row_l.addWidget(swatch)
            row_l.addWidget(QLabel(lbl))
            row_l.addStretch()
            self._dirs_check_layout.addWidget(row_w)
            self._dir_checks[tuple(uvw)] = chk

    # ------------------------------------------------------------------
    # Custom plane / direction
    # ------------------------------------------------------------------

    def _add_custom_plane(self):
        txt = self._custom_plane_edit.text().strip()
        if not txt:
            return
        try:
            h, k, l = _parse_hkl(txt)
        except ValueError as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Invalid input", str(e))
            return
        if (h, k, l) in self._face_checks:
            return   # already listed
        lbl = _plane_label(h, k, l)
        ci  = len(self._face_checks)
        col = _PLANE_COLORS[ci % len(_PLANE_COLORS)]
        chk = QCheckBox()
        chk.setChecked(True)
        row_w = QWidget()
        row_l = QHBoxLayout(row_w)
        row_l.setContentsMargins(0, 0, 0, 0)
        swatch = QLabel("  ")
        swatch.setFixedWidth(14)
        swatch.setStyleSheet(f"background:{col}; border-radius:3px;")
        row_l.addWidget(chk)
        row_l.addWidget(swatch)
        row_l.addWidget(QLabel(lbl))
        row_l.addStretch()
        self._planes_check_layout.addWidget(row_w)
        self._face_checks[(h, k, l)] = chk
        self._custom_plane_edit.clear()

    def _add_custom_dir(self):
        txt = self._custom_dir_edit.text().strip()
        if not txt:
            return
        try:
            u, v, w = _parse_hkl(txt)
        except ValueError as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Invalid input", str(e))
            return
        key = (u, v, w)
        if key in self._dir_checks:
            return
        lbl = _dir_label(u, v, w)
        dir_colors = ["#f39c12","#e74c3c","#9b59b6","#1abc9c","#e67e22","#f1c40f"]
        col = dir_colors[len(self._dir_checks) % len(dir_colors)]
        chk = QCheckBox()
        chk.setChecked(True)
        row_w = QWidget()
        row_l = QHBoxLayout(row_w)
        row_l.setContentsMargins(0, 0, 0, 0)
        swatch = QLabel("→")
        swatch.setStyleSheet(f"color:{col}; font-weight:bold; font-size:14px;")
        swatch.setFixedWidth(18)
        row_l.addWidget(chk)
        row_l.addWidget(swatch)
        row_l.addWidget(QLabel(lbl))
        row_l.addStretch()
        self._dirs_check_layout.addWidget(row_w)
        self._dir_checks[key] = chk
        self._custom_dir_edit.clear()

    # ------------------------------------------------------------------
    # Render
    # ------------------------------------------------------------------

    def _render(self):
        p = getattr(self._pv_widget, "_plotter", None)
        if p is None:
            return
        name   = self.crystal_combo.currentText()
        faces  = [hkl for hkl, chk in self._face_checks.items() if chk.isChecked()]
        dirs   = [list(uvw) for uvw, chk in self._dir_checks.items() if chk.isChecked()]
        try:
            _build_crystal_scene(
                p, name,
                show_faces=faces,
                show_dirs=dirs,
                show_cell=self.cell_check.isChecked(),
                face_opacity=self.opacity_spin.value(),
                show_labels=self.label_check.isChecked(),
                bg_color=self._bg_color,
                cell_color=self._cell_color,
            )
        except Exception as exc:
            print(f"[CrystalTab] render error: {exc}")


# ---------------------------------------------------------------------------
# Lightweight function to compute perceptual lightness for button contrast
# ---------------------------------------------------------------------------

def _lightness(hex_color: str) -> float:
    h = hex_color.lstrip("#")
    if len(h) < 6:
        return 0.5
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return (0.299 * r + 0.587 * g + 0.114 * b) / 255.0


# ---------------------------------------------------------------------------
# PyVista widget factory
# ---------------------------------------------------------------------------

class _PvContainer(QWidget):
    """Thin wrapper that embeds a pyvistaqt interactor (or fallback label)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._plotter = None
        try:
            from pyvistaqt import BackgroundPlotter
            # Do NOT pass parent= — newer pyvistaqt passes it positionally
            # inside QtInteractor and a duplicate keyword causes TypeError.
            p = BackgroundPlotter(show=False)
            p.set_background("#1e1e2e")
            self._plotter = p
            layout.addWidget(p)
        except Exception as e:
            lbl = QLabel(
                f"PyVista / pyvistaqt not available.\n\n{e}\n\n"
                "Install with:\n  pip install pyvista pyvistaqt"
            )
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setWordWrap(True)
            layout.addWidget(lbl)


def _make_pv_widget(parent=None) -> QWidget:
    return _PvContainer(parent)
