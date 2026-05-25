import numpy as np
from vtk import (
    vtkStructuredGrid, vtkPoints, vtkFloatArray,
    vtkStructuredGridGeometryFilter, vtkPolyDataMapper,
    vtkActor, vtkRenderer, vtkRenderWindow,
    vtkRenderWindowInteractor
)
from ..tensor.tensor import Tensor
from .utils import christoffel_tensor, wave_property


def _compute_velocities(c, rho):
    tensor = Tensor()
    cijkl = tensor.voigt_to_tensor(c)
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    theta, phi = np.meshgrid(theta, phi)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    vp = np.zeros((100, 100))
    vs1 = np.zeros((100, 100))
    vs2 = np.zeros((100, 100))
    for i in range(100):
        for j in range(100):
            n = np.array([x[i, j], y[i, j], z[i, j]])
            tik = christoffel_tensor(cijkl, n)
            wave_moduli, _ = wave_property(tik)
            vp[i, j] = np.sqrt(wave_moduli[0] / rho)
            vs1[i, j] = np.sqrt(wave_moduli[1] / rho)
            vs2[i, j] = np.sqrt(wave_moduli[2] / rho)
    return x, y, z, vp, vs1, vs2


def _render(x, y, z, scalar_xyz, scalar_color):
    grid = vtkStructuredGrid()
    points = vtkPoints()
    values = vtkFloatArray()
    for i in range(100):
        for j in range(100):
            s = scalar_xyz[i, j]
            points.InsertNextPoint(x[i, j] * s, y[i, j] * s, z[i, j] * s)
            values.InsertNextValue(scalar_color[i, j])
    grid.SetDimensions(100, 100, 1)
    grid.SetPoints(points)
    grid.GetPointData().SetScalars(values)
    gf = vtkStructuredGridGeometryFilter()
    gf.SetInputData(grid)
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(gf.GetOutputPort())
    mapper.SetScalarRange(np.min(scalar_color), np.max(scalar_color))
    actor = vtkActor()
    actor.SetMapper(mapper)
    renderer = vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1)
    render_window = vtkRenderWindow()
    render_window.AddRenderer(renderer)
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    interactor.Initialize()
    interactor.Start()


class Plotter:
    """A class for plotting seismic wave velocities in 3D using VTK."""

    @classmethod
    def plot_vs_splitting(cls, c, rho):
        """Plot Vs1-Vs2 splitting in 3D using VTK."""
        x, y, z, vp, vs1, vs2 = _compute_velocities(c, rho)
        _render(x, y, z, scalar_xyz=vs1 - vs2, scalar_color=vp)

    @classmethod
    def plot_vp(cls, c, rho):
        """Plot P-wave velocity (Vp) in 3D using VTK."""
        x, y, z, vp, vs1, vs2 = _compute_velocities(c, rho)
        _render(x, y, z, scalar_xyz=vp, scalar_color=vp)

    @classmethod
    def plot_vs1(cls, c, rho):
        """Plot fast S-wave velocity (Vs1) in 3D using VTK."""
        x, y, z, vp, vs1, vs2 = _compute_velocities(c, rho)
        _render(x, y, z, scalar_xyz=vs1, scalar_color=vs1)

    @classmethod
    def plot_vs2(cls, c, rho):
        """Plot slow S-wave velocity (Vs2) in 3D using VTK."""
        x, y, z, vp, vs1, vs2 = _compute_velocities(c, rho)
        _render(x, y, z, scalar_xyz=vs2, scalar_color=vs2)
