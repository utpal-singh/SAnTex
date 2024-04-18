import numpy as np
from vtk import (
    vtkStructuredGrid, vtkPoints, vtkFloatArray,
    vtkStructuredGridGeometryFilter, vtkPolyDataMapper,
    vtkActor, vtkRenderer, vtkRenderWindow,
    vtkRenderWindowInteractor
)
# from .anisotropy import Anisotropy
# from cijkl_3dplots import christoffel_tensor, wave_property
# import tensor_conversion

from .utils import christoffel_tensor, wave_property, tensor_to_voigt, voigt_to_tensor

class Plotter:
    @classmethod
    def plot_vs_splitting(cls, c, rho):
        cijkl = voigt_to_tensor(c)

        theta = np.linspace(0, np.pi, 100)
        phi = np.linspace(0, 2 * np.pi, 100)
        theta, phi = np.meshgrid(theta, phi)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        vp = np.zeros((100, 100))
        vs1 = np.zeros((100, 100))
        vs2 = np.zeros((100, 100))

        for i in range(len(theta)):
            for j in range(len(phi)):
                n = np.array([np.sin(theta[i][j]) * np.cos(phi[i][j]),
                              np.sin(theta[i][j]) * np.sin(phi[i][j]),
                              np.cos(theta[i][j])])
                tik = christoffel_tensor(cijkl, n)
                wave_moduli, _ = wave_property(tik)
                vp[i][j] = np.sqrt(wave_moduli[0] / rho)
                vs1[i][j] = np.sqrt(wave_moduli[1] / rho)
                vs2[i][j] = np.sqrt(wave_moduli[2] / rho)

        grid = vtkStructuredGrid()
        points = vtkPoints()
        values = vtkFloatArray()

        for i in range(len(theta)):
            for j in range(len(phi)):
                points.InsertNextPoint(x[i][j] * (vs1[i][j] - vs2[i][j]),
                                       y[i][j] * (vs1[i][j] - vs2[i][j]),
                                       z[i][j] * (vs1[i][j] - vs2[i][j]))
                values.InsertNextValue(vp[i][j])

        grid.SetDimensions(100, 100, 1)
        grid.SetPoints(points)
        grid.GetPointData().SetScalars(values)

        geometry_filter = vtkStructuredGridGeometryFilter()
        geometry_filter.SetInputData(grid)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(geometry_filter.GetOutputPort())
        mapper.SetScalarRange(np.min(vp), np.max(vp))

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

    @classmethod
    def plot_vp(cls, c, rho):
        cijkl = voigt_to_tensor(c)

        theta = np.linspace(0, np.pi, 100)
        phi = np.linspace(0, 2 * np.pi, 100)
        theta, phi = np.meshgrid(theta, phi)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        vp = np.zeros((100, 100))
        vs1 = np.zeros((100, 100))
        vs2 = np.zeros((100, 100))

        for i in range(len(theta)):
            for j in range(len(phi)):
                n = np.array([np.sin(theta[i][j]) * np.cos(phi[i][j]),
                              np.sin(theta[i][j]) * np.sin(phi[i][j]),
                              np.cos(theta[i][j])])
                tik = christoffel_tensor(cijkl, n)
                wave_moduli, _ = wave_property(tik)
                vp[i][j] = np.sqrt(wave_moduli[0] / rho)
                vs1[i][j] = np.sqrt(wave_moduli[1] / rho)
                vs2[i][j] = np.sqrt(wave_moduli[2] / rho)

        grid = vtkStructuredGrid()
        points = vtkPoints()
        values = vtkFloatArray()

        for i in range(len(theta)):
            for j in range(len(phi)):
                points.InsertNextPoint(x[i][j] * (vp[i][j]),
                                       y[i][j] * (vp[i][j]),
                                       z[i][j] * (vp[i][j]))
                values.InsertNextValue(vp[i][j])

        grid.SetDimensions(100, 100, 1)
        grid.SetPoints(points)
        grid.GetPointData().SetScalars(values)

        geometry_filter = vtkStructuredGridGeometryFilter()
        geometry_filter.SetInputData(grid)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(geometry_filter.GetOutputPort())
        mapper.SetScalarRange(np.min(vp), np.max(vp))

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

    @classmethod
    def plot_vs1(cls, c, rho):
        cijkl = voigt_to_tensor(c)

        theta = np.linspace(0, np.pi, 100)
        phi = np.linspace(0, 2 * np.pi, 100)
        theta, phi = np.meshgrid(theta, phi)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        vp = np.zeros((100, 100))
        vs1 = np.zeros((100, 100))
        vs2 = np.zeros((100, 100))

        for i in range(len(theta)):
            for j in range(len(phi)):
                n = np.array([np.sin(theta[i][j]) * np.cos(phi[i][j]),
                              np.sin(theta[i][j]) * np.sin(phi[i][j]),
                              np.cos(theta[i][j])])
                tik = christoffel_tensor(cijkl, n)
                wave_moduli, _ = wave_property(tik)
                vp[i][j] = np.sqrt(wave_moduli[0] / rho)
                vs1[i][j] = np.sqrt(wave_moduli[1] / rho)
                vs2[i][j] = np.sqrt(wave_moduli[2] / rho)

        grid = vtkStructuredGrid()
        points = vtkPoints()
        values = vtkFloatArray()

        for i in range(len(theta)):
            for j in range(len(phi)):
                points.InsertNextPoint(x[i][j] * (vs1[i][j]),
                                       y[i][j] * (vs1[i][j]),
                                       z[i][j] * (vs1[i][j]))
                values.InsertNextValue(vp[i][j])

        grid.SetDimensions(100, 100, 1)
        grid.SetPoints(points)
        grid.GetPointData().SetScalars(values)

        geometry_filter = vtkStructuredGridGeometryFilter()
        geometry_filter.SetInputData(grid)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(geometry_filter.GetOutputPort())
        mapper.SetScalarRange(np.min(vp), np.max(vp))

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
    @classmethod
    def plot_vp(cls, c, rho):
        cijkl = voigt_to_tensor(c)

        theta = np.linspace(0, np.pi, 100)
        phi = np.linspace(0, 2 * np.pi, 100)
        theta, phi = np.meshgrid(theta, phi)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        vp = np.zeros((100, 100))
        vs1 = np.zeros((100, 100))
        vs2 = np.zeros((100, 100))

        for i in range(len(theta)):
            for j in range(len(phi)):
                n = np.array([np.sin(theta[i][j]) * np.cos(phi[i][j]),
                              np.sin(theta[i][j]) * np.sin(phi[i][j]),
                              np.cos(theta[i][j])])
                tik = christoffel_tensor(cijkl, n)
                wave_moduli, _ = wave_property(tik)
                vp[i][j] = np.sqrt(wave_moduli[0] / rho)
                vs1[i][j] = np.sqrt(wave_moduli[1] / rho)
                vs2[i][j] = np.sqrt(wave_moduli[2] / rho)

        grid = vtkStructuredGrid()
        points = vtkPoints()
        values = vtkFloatArray()

        for i in range(len(theta)):
            for j in range(len(phi)):
                points.InsertNextPoint(x[i][j] * (vs2[i][j]),
                                       y[i][j] * (vs2[i][j]),
                                       z[i][j] * (vs2[i][j]))
                values.InsertNextValue(vp[i][j])

        grid.SetDimensions(100, 100, 1)
        grid.SetPoints(points)
        grid.GetPointData().SetScalars(values)

        geometry_filter = vtkStructuredGridGeometryFilter()
        geometry_filter.SetInputData(grid)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(geometry_filter.GetOutputPort())
        mapper.SetScalarRange(np.min(vp), np.max(vp))

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



if __name__ == "__main__":
    density = 3310
    c = np.array([[323.70, 66.40, 71.60, 0.000, 0.000, 0.000],
                [66.40, 197.60, 75.60, 0.000, 0.000, 0.000],
                [71.60, 75.60, 235.10, 0.000, 0.000, 0.000],
                [0.000, 0.000, 0.000, 64.62, 0.000, 0.000],
                [0.000, 0.000, 0.000, 0.000, 78.05, 0.000],
                [0.000, 0.000, 0.000, 0.000, 0.000, 79.04]]) * 1e9
    # cijkl = tensor_conversion.voigt_to_tensor(c)

    Plotter.plot_wave_velocities(c, density)