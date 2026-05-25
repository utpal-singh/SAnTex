import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt5.QtCore import Qt


def _make_pyvista_widget(parent=None):
    """Return a pyvistaqt.BackgroundPlotter or a fallback placeholder."""
    try:
        from pyvistaqt import BackgroundPlotter
        plotter = BackgroundPlotter(show=False)
        plotter.set_background("white")
        return plotter, True
    except ImportError:
        placeholder = QLabel(
            "PyVista / pyvistaqt not installed.\nInstall with: pip install pyvista pyvistaqt",
            parent=parent,
        )
        placeholder.setAlignment(Qt.AlignCenter)
        return placeholder, False


class PyVistaWidget(QWidget):
    """Wrapper that embeds a pyvistaqt BackgroundPlotter (or a placeholder)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        self._view, self._has_pyvista = _make_pyvista_widget(self)
        self._layout.addWidget(self._view)

    @property
    def plotter(self):
        """Return the BackgroundPlotter, or None if pyvista is unavailable."""
        return self._view if self._has_pyvista else None

    def clear(self):
        if self._has_pyvista:
            self._view.clear()

    def plot_velocity_surface(self, surface_data: dict, scalar: str = "vp"):
        """
        Render a velocity surface from surface_data dict produced by
        AnisotropyBackend.compute_velocity_surface().

        Parameters
        ----------
        surface_data : dict with keys nx, ny, nz, vp, vs1, vs2, avs, vpvs1, vpvs2
        scalar       : which quantity to color ('vp', 'vs1', 'vs2', 'avs', 'vpvs1', 'vpvs2')
        """
        if not self._has_pyvista:
            return
        try:
            import pyvista as pv
        except ImportError:
            return

        nx = surface_data["nx"].ravel()
        ny = surface_data["ny"].ravel()
        nz = surface_data["nz"].ravel()
        v  = surface_data[scalar].ravel()

        # Scale unit-sphere directions by the velocity magnitude so the shape
        # reflects wave-speed variation (in km/s for readability)
        v_km = v / 1000.0
        pts = np.column_stack([nx * v_km, ny * v_km, nz * v_km])

        n_theta = surface_data["vp"].shape[0]
        n_phi   = surface_data["vp"].shape[1]

        # Build a structured surface grid from the (n_theta × n_phi) parametric grid
        grid = pv.StructuredGrid()
        grid.points = pts
        grid.dimensions = [n_phi, n_theta, 1]
        grid[scalar] = v_km

        self._view.clear()
        self._view.add_mesh(grid, scalars=scalar, cmap="seismic",
                            show_scalar_bar=True, smooth_shading=True)
        self._view.add_axes()
        self._view.reset_camera()

    def plot_3d_stereonet(self, surface_data: dict, scalar: str = "vp"):
        """Simple scatter version (fallback for non-structured grids)."""
        if not self._has_pyvista:
            return
        try:
            import pyvista as pv
        except ImportError:
            return

        nx = surface_data["nx"].ravel()
        ny = surface_data["ny"].ravel()
        nz = surface_data["nz"].ravel()
        v  = surface_data[scalar].ravel()
        v_km = v / 1000.0

        pts = np.column_stack([nx * v_km, ny * v_km, nz * v_km])
        cloud = pv.PolyData(pts)
        cloud[scalar] = v_km

        self._view.clear()
        self._view.add_mesh(cloud, scalars=scalar, cmap="seismic",
                            point_size=3, render_points_as_spheres=True,
                            show_scalar_bar=True)
        self._view.add_axes()
        self._view.reset_camera()
