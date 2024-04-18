import numpy as np
import math
from santex import Tensor

import matplotlib.pyplot as plt
import numpy as np
from .vtkplotter import Plotter
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .plot_vel_grid import plot_velocity_grid

from .. import Material
from .. import EBSD

from scipy.interpolate import griddata
import math

class Anisotropy:
    """
    Anisotropy class for representing anisotropic material properties.

    Attributes:
        cijkl (numpy.ndarray or None): stiffness tensor in voigt notation.
            It is converted from the input stiffness matrix using voigt_to_tensor method.
            If stiffness_matrix is None, cijkl remains None. Units can be SI or CGS but should remain consistent
            between stiffness matrix and density
        rho (float or None): Density of the material.
            If density is None, rho remains None.

    Methods:
        __init__(stiffness_matrix, density):
            Initializes the Anisotropy object with the given stiffness_matrix and density.
            If either stiffness_matrix or density is None, cijkl and rho remain None.

    Example Usage:
        stiffness_matrix = np.array([[198.96, 73.595, 68.185, 0., 9.735, 0.],
                                [73.595, 155.94, 62.23, 0., 6.295, 0.],
                                [68.185, 62.23, 225.99, 0., 33.85, 0.],
                                [0., 0., 0., 65.66, 0., 6.415],
                                [9.735, 6.295, 33.85, 0., 60.23, 0.],
                                [0., 0., 0., 6.415, 0., 65.18]]) * 10**9
        density = 3.5 
        anisotropy = Anisotropy(stiffness_matrix, density)
    """
    
    def __init__(self, stiffness_matrix, density):
        """
        Initializes an Anisotropy object with the given stiffness_matrix and density.

        Parameters:
            stiffness_matrix (list of lists or numpy.ndarray): The 6x6 stiffness matrix representing
                the material's anisotropic properties in Voigt notation.
            density (float): The density of the material.

        Raises:
            TypeError: If stiffness_matrix is not a 6x6 list or numpy.ndarray.

        """
        self.cijkl = None
        self.rho = None

        if stiffness_matrix is not None and density is not None:
            tensor_object = Tensor()
            self.cijkl = tensor_object.voigt_to_tensor(stiffness_matrix)
            self.rho = density


    def christoffel_tensor(self, n):
        """
        Calculates the Christoffel tensor given a direction vector n.

        Parameters:
            n (numpy.ndarray): The direction vector for which the Christoffel tensor is calculated.

        Returns:
            numpy.ndarray: The Christoffel tensor Tik for the given direction vector n.

        Raises:
            ValueError: If an error occurs during the calculation.

        Example direction:
            n = np.array([1, 0, 0])
        
        """
        try:
            tik = np.zeros((3, 3))

            for i in range(3):
                for k in range(3):
                    tik[i, k] = np.tensordot(self.cijkl[i, :, k, :], np.outer(n, n))

            return tik
        except Exception as e:
            raise ValueError("Error in calculating the Christoffel tensor:", e)

    def wave_property(self, tik):
        """
        Calculates the wave properties (wave moduli and polarization directions) given the Christoffel tensor Tik.

        Parameters:
            tik (numpy.ndarray): The Christoffel tensor Tik for which wave properties are calculated.

        Returns:
            tuple: A tuple containing the wave moduli and polarization directions.

        Raises:
            ValueError: If an error occurs during the calculation.
        """
        try:
            eigenvalues, eigenvectors = np.linalg.eig(tik)
            indices = np.argsort(eigenvalues)[::-1]

            wave_moduli = [eigenvalues[i] for i in indices]
            polarization_directions = [eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]) for i in indices]

            return tuple(wave_moduli), tuple(polarization_directions)
        except Exception as e:
            raise ValueError("Error in calculating wave properties:", e)

    def phase_velocity(self):
        """
        Calculates the phase velocities (P-wave velocity, S1-wave velocity, and S2-wave velocity) for different
        directions using the Christoffel tensor and wave properties.

        Returns:
            tuple: A tuple containing the phase velocities (vp for P-wave, vs1 for S1-wave, and vs2 for S2-wave).

        Raises:
            ValueError: If an error occurs during the calculation.    
        """
        try:
            vp = []
            vs1 = []
            vs2 = []

            theta_values = np.arange(0, math.pi + math.pi / 180, math.pi / 180)
            phi_values = np.arange(0, 2 * math.pi + math.pi / 180, math.pi / 180)

            theta, phi = np.meshgrid(theta_values, phi_values, indexing='ij')

            sin_theta = np.sin(theta)
            cos_theta = np.cos(theta)
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)

            sin_theta_sin_phi = sin_theta * cos_phi
            sin_theta_cos_phi = sin_theta * sin_phi

            sqrt_rho = np.sqrt(self.rho)

            for i in range(theta.shape[0]):
                for j in range(phi.shape[1]):
                    n = np.array([sin_theta_sin_phi[i, j], sin_theta_cos_phi[i, j], cos_theta[i, j]])
                    tik = self.christoffel_tensor(n)
                    if tik is not None:
                        wave_moduli, polarization_directions = self.wave_property(tik)
                        vp.append(np.sqrt(wave_moduli[0] / self.rho))
                        vs1.append(np.sqrt(wave_moduli[1] / self.rho))
                        vs2.append(np.sqrt(wave_moduli[2] / self.rho))

            return tuple(vp), tuple(vs1), tuple(vs2)
        except Exception as e:
            raise ValueError("Error in calculating phase velocity:", e)
        
    def velocities(self):
        vp, vs1, vs2 = self.phase_velocity()
        return vp, vs1, vs2

    def anisotropy_values(self, stiffness_matrix = None, density = None, method = None, return_values=None):
        """
        Calculates various anisotropy values based on the velocities calculated from the given stiffness matrix and density.

        Parameters:
            stiffness_matrix (list or None): The stiffness matrix representing the material's anisotropic properties.
            density (float or None): The density of the material.
            method (str or None): The method to use for calculating anisotropy values. Options: 'array'.
            return_values (str or None): The specific anisotropy value to return. Options: 'maxvp', 'minvp', 'maxvs1',
                'minvs1', 'maxvs2', 'minvs2', or None (default) to print all values.

        Returns:
            float or None: A dictionary containing the calculated anisotropy values if return_values is None,
                or a single anisotropy value if return_values is specified.

        Raises:
            ValueError: If an error occurs during the calculation.

        Notes:
            - If method='array', an array of Anisotropy objects can be provided for batch calculation.
        
        """

        if method == "array":
            anis_mat = []
            for i in stiffness_matrix:
                anis = Anisotropy(stiffness_matrix[i], density[i])
                anis_mat.append(self.anisotropy_values(anis))
                self.anisotropy_values(anis)
            return anis_mat

        vp, vs1, vs2 = self.velocities()
        maxvp = max(vp)
        minvp = min(vp)
        maxvs1 = max(vs1)
        minvs1 = min(vs1)
        maxvs2 = max(vs2)
        minvs2 = min(vs2)
        swaveAnisotropy_percent = 200 * (np.array(vs1) - np.array(vs2)) / (np.array(vs1) + np.array(vs2))
        max_vs_anisotropy_percent = max(swaveAnisotropy_percent)
        min_vs_anisotropy_percent = min(swaveAnisotropy_percent)
        p_wave_anisotropy_percent = 200 * (maxvp - minvp) / (maxvp + minvp)
        s1_wave_anisotropy_percent = 200 * (maxvs1 - minvs1) / (maxvs1 + minvs1)
        s2_wave_anisotropy_percent = 200 * (maxvs2 - minvs2) / (maxvs2 + minvs2)
        dvs = np.array(vs1) - np.array(vs2)
        maxdvs = max(dvs)
        vp_vs1 = np.array(vp) / np.array(vs1)
        AVpVs1 = 200 * (max(vp_vs1) - min(vp_vs1)) / (max(vp_vs1) + min(vp_vs1))

        if return_values == 'maxvp':
            return maxvp
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'maxvs1':
            return maxvs1
        
        else:
            print("Max Vp: ", maxvp)
            print("Min Vp: ", minvp)
            print("Max Vs1: ", maxvs1)
            print("Min Vs1: ", minvs1)
            print("Max Vs2: ", maxvs2)
            print("Min Vs2: ", minvs2)
            print("Max vs anisotropy percent: ", max_vs_anisotropy_percent)
            print("Min vs anisotropy percent: ", min_vs_anisotropy_percent)
            print("P wave anisotropy percent: ", p_wave_anisotropy_percent)
            print("S1 Wave anisotropy percent: ", s1_wave_anisotropy_percent)
            print("S2 Wave anisotropy percent: ", s2_wave_anisotropy_percent)
            print("Velocity difference: ", maxdvs)
            print("Vp/Vs1 ratio: ", AVpVs1)
            return {
                'maxvp': maxvp,
                'minvp': minvp,
                'maxvs1': maxvs1,
                'minvs1': minvs1,
                'maxvs2': maxvs2,
                'minvs2': minvs2,
                'max_vs_anisotropy_percent': max_vs_anisotropy_percent,
                'min_vs_anisotropy_percent': min_vs_anisotropy_percent,
                'p_wave_anisotropy_percent': p_wave_anisotropy_percent,
                's1_wave_anisotropy_percent': s1_wave_anisotropy_percent,
                's2_wave_anisotropy_percent': s2_wave_anisotropy_percent,
                'maxdvs': maxdvs,
                'AVpVs1': AVpVs1
            }
        
    def plot_velocities(self, pressure_range, temperature_range, return_type, is_ebsd = False, phase = None, grid = [5, 5], filename = None, *args):
        """
        Plots velocities based on specified ranges and return types.

            Parameters:
                pressure_range (tuple): The range of pressures for which velocities will be plotted.
                temperature_range (tuple): The range of temperatures for which velocities will be plotted.
                return_type (str): The type of velocity to plot. Options: 'maxvp', 'minvp', 'maxvs1', 'minvs1', 'maxvs2', 'minvs2'.
                is_ebsd (bool): Whether the data comes from electron backscatter diffraction (EBSD). Default is False.
                phase (str or None): The phase of the material. Only required if is_ebsd is True.
                grid (list): The grid dimensions for the plot. Default is [5, 5].
                filename (str or None): The filename to save the plot. Required if is_ebsd is True.
                *args: can be [01, 2, 3]

            Raises:
                ValueError: If required parameters are not provided or if an error occurs during plotting.

            Notes:
                - If is_ebsd is True, phase and filename must be provided.
                - The plot type and appearance can be customized using *args.

            Example Usage:
                anisotropy.plot_velocities((0, 100), (500, 1000), 'maxvp', is_ebsd=True, phase='phase1', filename='velocity_plot.png', 'ro-')
        """
        plot_velocity_grid(pressure_range, temperature_range, return_type, is_ebsd = False, phase = phase, grid = [5, 5], filename = None, *args)



    def plot(self, colormap="RdBu_r", step = 180, savefig = False, figname = None, dpi = 300):
        """
        Plots various anisotropic maps based on the Christoffel tensor.

        Parameters:
            colormap (str): The colormap to use for plotting. Default is "RdBu_r".
            step (int): The step size for theta and phi values. Default is 180.
            savefig (bool): Whether to save the plot as an image. Default is False.
            figname (str or None): The filename to save the plot. Required if savefig is True.
            dpi (int): The resolution of the saved image. Default is 300.

        Raises:
            ValueError: If an error occurs during the plotting process.

        Notes:
            - This method generates a 2x3 grid of subplots, each representing different anisotropic maps based on the Christoffel tensor.
            - The colormap, step size, and other parameters can be customized.

        Example Usage:
            anisotropy.plot(colormap="viridis", step=120, savefig=True, figname="anisotropy_plot.png", dpi=600)
        """
        try:
            fig, axs = plt.subplots(2, 3, figsize=(15, 10))
            step = math.pi / step

            # texts for each subplot
            texts = ['Ratio of VP to VS1', 'Velocity of P-waves (VP)', 'Velocity of S1-waves (VS1)', 
                    'Velocity of S2-waves (VS2)', 'Anisotropy measure for VP and VS1', 
                    'Anisotropy measure for VP and VS2']
            
            for i, ax in enumerate(axs.flat):
                x = []
                y = []
                c = []

                for theta in np.arange(0, math.pi / 2 + step, step):
                    for phi in np.arange(0, 2 * math.pi + step, step):
                        n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
                        tik = self.christoffel_tensor(n)
                        if tik is not None:
                            wave_moduli, _ = self.wave_property(tik)
                            vp = math.sqrt(wave_moduli[0] / self.rho)
                            vs1 = math.sqrt(wave_moduli[1] / self.rho)
                            vs2 = math.sqrt(wave_moduli[2] / self.rho)
                            vpvs1 = vp/vs1
                            x.append(n[0] / (1 + n[2]))
                            y.append(n[1] / (1 + n[2]))
                            if i == 0:
                                c.append(vpvs1)
                            elif i == 1:
                                c.append(vp)
                            elif i == 2:
                                c.append(vs1)
                            elif i == 3:
                                c.append(vs2)
                            elif i == 4:
                                c.append((vp - vs1) / (vp + vs1))
                            elif i == 5:
                                c.append((vp - vs2) / (vp + vs2))

                # Interpolate onto a regular grid
                xi = np.linspace(min(x), max(x), 100)
                yi = np.linspace(min(y), max(y), 100)
                xi, yi = np.meshgrid(xi, yi)
                zi = griddata((x, y), c, (xi, yi), method='linear')

                # Plotting contour lines for each subplot
                contours = ax.contour(xi, yi, zi, 5, colors='black')
                ax.clabel(contours, inline=True, fontsize=8)

                sc = ax.scatter(x, y, c=c, cmap=colormap, s=5)  # Reduce scatter dot size s for clarity
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect('equal', 'box')

                # Adding above text at the bottom of each plot
                ax.text(0.5, -0.15, texts[i], ha='center', transform=ax.transAxes)

            axs[0, 0].set_title('VP/VS1')
            axs[0, 1].set_title('VP')
            axs[0, 2].set_title('VS1')
            axs[1, 0].set_title('VS2')
            axs[1, 1].set_title('AVpVs1')
            axs[1, 2].set_title('AVpVs2')

            plt.tight_layout()
            plt.show()

            if savefig:
                plt.savefig(f"{figname}", dpi = dpi)
        except Exception as e:
            print(f"An error occurred: {e}")


    def plotly(self):
        try:
            fig = go.Figure()

            step = math.pi / 180

            fig = make_subplots(rows=2, cols=3, subplot_titles=("VP/VS1", "VP", "VS1", "VS2", "AVpVs1", "AVpVs2"))

            for i in range(6):
                x = []
                y = []
                c = []

                for theta in np.arange(0, math.pi / 2 + step, step):
                    for phi in np.arange(0, 2 * math.pi + step, step):
                        n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
                        tik = self.christoffel_tensor(n)
                        if tik is not None:
                            wave_moduli, _ = self.wave_property(tik)
                            vp = math.sqrt(wave_moduli[0] / self.rho)
                            vs1 = math.sqrt(wave_moduli[1] / self.rho)
                            vs2 = math.sqrt(wave_moduli[2] / self.rho)
                            vpvs1 = vp/vs1
                            x.append(n[0] / (1 + n[2]))
                            y.append(n[1] / (1 + n[2]))
                            if i == 0:
                                c.append(vpvs1)
                            elif i == 1:
                                c.append(vp)
                            elif i == 2:
                                c.append(vs1)
                            elif i == 3:
                                c.append(vs2)
                            elif i == 4:
                                c.append((vp - vs1) / (vp + vs1))
                            elif i == 5:
                                c.append((vp - vs2) / (vp + vs2))

                row = i // 3 + 1
                col = i % 3 + 1
                fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode='markers',
                    marker=dict(
                        size=5,
                        color=c,
                        colorscale='RdBu',
                        colorbar=dict(title='Colorbar Title'),
                    ),
                    name='',
                ), row=row, col=col)

            fig.update_layout(
                title="Plot Title",
                xaxis_title="x",
                yaxis_title="y"
            )

            fig.show()
        except Exception as e:
            raise ValueError("Error in plotting:", e)


            
    def plotter_vs_splitting(self, density, voigt_stiffness):

        Plotter.plot_vs_splitting(voigt_stiffness, density)

    def plotter_vp(self, density, voigt_stiffness):

        Plotter.plot_vp(voigt_stiffness, density)

    def plotter_vs1(self, density, voigt_stiffness):

        Plotter.plot_vp(voigt_stiffness, density)

    def plotter_vs2(self, density, voigt_stiffness):

        Plotter.plot_vp(voigt_stiffness, density)
