import numpy as np
import math
from ..tensor import Tensor
from matplotlib import pyplot as plt
import numpy as np
from .vtkplotter import Plotter
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .plot_vel_grid import plot_velocity_grid
from .utils import christoffel_tensor, wave_property

from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

import math

class Anisotropy:
    """

    A class to represent anisotropic material properties.

    Attributes:

        cijkl (Optional[numpy.ndarray]): Stiffness tensor in Voigt notation.
            It is converted from the input stiffness matrix using the `voigt_to_tensor` method.
            If `stiffness_matrix` is None, `cijkl` remains None. Units should be consistent between the 
            stiffness matrix and density, either in SI or CGS.
        rho (Optional[float]): Density of the material.
            If `density` is None, `rho` remains None.

    Parameters:

        stiffness_matrix (Optional[numpy.ndarray]): A 6x6 matrix representing the stiffness tensor in Voigt notation.
        density (Optional[float]): The density of the material.

    Methods:

        __init__(stiffness_matrix, density):
            Initializes the Anisotropy object with the given stiffness_matrix and density.
            If either `stiffness_matrix` or `density` is None, `cijkl` and `rho` remain None.

    Example:

        >>> import numpy as np
        >>> stiffness_matrix = np.array([[198.96, 73.595, 68.185, 0., 9.735, 0.],
        >>>                              [73.595, 155.94, 62.23, 0., 6.295, 0.],
        >>>                              [68.185, 62.23, 225.99, 0., 33.85, 0.],
        >>>                              [0., 0., 0., 65.66, 0., 6.415],
        >>>                              [9.735, 6.295, 33.85, 0., 60.23, 0.],
        >>>                              [0., 0., 0., 6.415, 0., 65.18]]) * 10**9
        >>> density = 3.5  # density in g/cm³ or consistent units
        >>> anisotropy = Anisotropy(stiffness_matrix, density)
        >>> print(anisotropy.cijkl)
        >>> print(anisotropy.rho)

    """

    
    def __init__(self, stiffness_matrix, density):
        """
        Initializes an Anisotropy object with the given stiffness matrix and density.

        Parameters:
            stiffness_matrix (list[list[float]] or numpy.ndarray): A 6x6 stiffness matrix representing
                the material's anisotropic properties in Voigt notation.
            density (float): The density of the material.

        Raises:
            TypeError: If `stiffness_matrix` is not a 6x6 list of lists or a numpy.ndarray.

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
            return christoffel_tensor(self.cijkl, n)
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
            return wave_property(tik)
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

    def anisotropy_values(self, stiffness_matrix=None, density=None, method=None, return_values=None):
        """
        Calculates various anisotropy values based on the velocities calculated from the given stiffness matrix and density.

        Parameters:
            stiffness_matrix (list or None): The stiffness matrix representing the material's anisotropic properties.
                If None, the current object's stiffness matrix is used.
            
            density (float or None): The density of the material. If None, the current object's density is used.
            
            method (str or None): The method to use for calculating anisotropy values. Options include 'array' for batch 
                processing of multiple anisotropy objects.
            
            return_values (str or None): The specific anisotropy value to return. Options: 'maxvp', 'minvp', 'maxvs1',
                'minvs1', 'maxvs2', 'minvs2', 'max_vs_anisotropy_percent', 'min_vs_anisotropy_percent',
                'p_wave_anisotropy_percent', 's1_wave_anisotropy_percent', 's2_wave_anisotropy_percent',
                'maxdvs', 'AVpVs1', or None (default) to print all values.

        Returns:
            float or dict: If `return_values` is specified, returns the corresponding anisotropy value. If None, 
            returns a dictionary containing all calculated anisotropy values.

        Raises:
            ValueError: If an error occurs during the calculation, such as invalid input values.

        Notes:
            - If `method='array'`, an array of Anisotropy objects can be provided for batch calculation.
            - The calculated anisotropy values include both P-wave and S-wave anisotropy percentages.

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
        meanvp = np.mean(vp)
        meanvs1 = np.mean(vs1)
        meanvs2 = np.mean(vs2)
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
        elif return_values == 'minvp':
            return minvp
        elif return_values == 'maxvs1':
            return maxvs1
        elif return_values == 'minvs1':
            return minvs1
        elif return_values == 'maxvs2':
            return maxvs2
        elif return_values == 'minvs2':
            return minvs2
        elif return_values == 'max_vs_anisotropy_percent':
            return max_vs_anisotropy_percent
        elif return_values == 'min_vs_anisotropy_percent':
            return min_vs_anisotropy_percent
        elif return_values == 'p_wave_anisotropy_percent':
            return p_wave_anisotropy_percent
        elif return_values == 's1_wave_anisotropy_percent':
            return s1_wave_anisotropy_percent
        elif return_values == 's2_wave_anisotropy_percent':
            return s2_wave_anisotropy_percent
        elif return_values == 'maxdvs':
            return maxdvs
        elif return_values == 'AVpVs1':
            return AVpVs1
        
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
            print("Velocity difference, maxdvs: ", maxdvs)
            print("Vp/Vs1 ratio: ", AVpVs1)
            print("Mean vp: ", meanvp)
            print("Mean Vs1: ", meanvs1)
            print("Mean Vs2: ", meanvs2)
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
        
    def plot_velocities(self, pressure_range, temperature_range, return_type, is_ebsd=False, phase=None, grid=[5, 5], filename=None, *args):
        """
        Plots velocities based on specified ranges and return types.

        Parameters:
            pressure_range (tuple): The range of pressures for which velocities will be plotted. Example: [2, 4].
            temperature_range (tuple): The range of temperatures for which velocities will be plotted. Example: [1000, 2000].
            return_type (str): The type of velocity to plot. Options: 'maxvp', 'minvp', 'maxvs1', 'minvs1', 'maxvs2', 'minvs2'.
            is_ebsd (bool): Whether the data comes from electron backscatter diffraction (EBSD). Default is False.
            phase (str or None): The phase of the material. Example: phase='Forsterite'.
            grid (list): The grid dimensions for the plot. Default is [5, 5].
            filename (str or None): The filename of the EBSD data. Required if `is_ebsd` is True.
            *args: Additional arguments for plot customization.

        Raises:
            ValueError: If required parameters are not provided or if an error occurs during plotting.

        Notes:
            - If `is_ebsd` is True, `phase` and `filename` must be provided.
            - The plot type and appearance can be customized using `*args`.

        Example:
            >>> anisotropy.plot_velocities((0, 100), (500, 1000), 'maxvp', is_ebsd=True, phase='phase1', filename='velocity_plot.png', 'ro-')
        """
        return plot_velocity_grid(pressure_range=pressure_range, temperature_range=temperature_range, return_type=return_type, is_ebsd=is_ebsd, phase=phase, grid=grid, filename=filename, *args)


    def plot(self, colormap="RdBu", step = 180, savefig = False, figname = None, dpi = 300, save_format = 'svg'):
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
                            vpvs2 = vp/vs2
                            avs = 200*(vs1-vs2)/(vs1+vs2)
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
                            # elif i == 4:
                            #     c.append((vp - vs1) / (vp + vs1))
                            elif i == 4:
                                c.append(vpvs2)
                            # elif i == 5:
                            #     c.append((vp - vs2) / (vp + vs2))
                            elif i == 5:
                                c.append(avs)

                # Interpolate onto a regular grid
                xi = np.linspace(min(x), max(x), 100)
                yi = np.linspace(min(y), max(y), 100)
                xi, yi = np.meshgrid(xi, yi)
                zi = griddata((x, y), c, (xi, yi), method='linear')

                # Plotting contour lines for each subplot
                contours = ax.contour(xi, yi, zi, 5, colors='black')
                ax.clabel(contours, inline=True, fontsize=8)

                sc = ax.scatter(x, y, c=c, cmap=colormap, s=5, rasterized=True)  # Reduce scatter dot size s for clarity
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect('equal', 'box')

                ax.text(0.5, -0.15, texts[i], ha='center', transform=ax.transAxes)
                # cbar = fig.colorbar(sc, ax=ax)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad = 0.1)
                cbar = fig.colorbar(sc, cax=cax, orientation='vertical')

            axs[0, 0].set_title('Vp/Vs1')
            axs[0, 1].set_title('Vp (m/s)')
            axs[0, 2].set_title('Vs1 (m/s)')
            axs[1, 0].set_title('Vs2')
            axs[1, 1].set_title('Vp/Vs2')
            axs[1, 2].set_title('S-wave Anisotropy (%)')

            plt.tight_layout()

            if savefig:
                plt.savefig(f"{figname}.{save_format}", dpi = dpi, format = save_format)
            
            plt.show()
        except Exception as e:
            print(f"An error occurred: {e}")


    def plotly(self):
        """
        Generates a Plotly figure to visualize material properties based on calculations using Christoffel and wave tensors.

        This method creates a Plotly figure with subplots for various properties such as VP/VS1, VP, VS1, VS2, AVpVs1, and AVpVs2.
        The figure is generated by iterating over specified ranges of angles and calculating properties like wave moduli, VP, VS1, VS2,
        VP/VS1 ratios, and AVpVs1/AVpVs2 ratios. The resulting data points are plotted on scatter plots with a color scale representing
        different property values.

        Raises:
            ValueError: If an error occurs during the plotting process.

        Notes:
            - Requires the existence of the `christoffel_tensor` and `wave_property` methods in the class instance to calculate material properties.

        Example Usage:
            plotter = anisotropy_instance.plotly()
        """
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
        """
        Plots the Vs splitting based on density and Voigt stiffness.

        Parameters:
            density (float): The density of the material.
            voigt_stiffness (float): The Voigt stiffness of the material.
        """
        Plotter.plot_vs_splitting(voigt_stiffness, density)

    def plotter_vp(self, density, voigt_stiffness):
        """
        Plots the Vp based on density and Voigt stiffness.

        Parameters:
            density (float): The density of the material.
            voigt_stiffness (float): The Voigt stiffness of the material.
        """

        Plotter.plot_vp(voigt_stiffness, density)

    def plotter_vs1(self, density, voigt_stiffness):
        """
        Plots the Vs1 based on density and Voigt stiffness.

        Parameters:
            density (float): The density of the material.
            voigt_stiffness (float): The Voigt stiffness of the material.
        """

        Plotter.plot_vp(voigt_stiffness, density)

    def plotter_vs2(self, density, voigt_stiffness):
        """
        Plots the Vs2 based on density and Voigt stiffness.

        Parameters:
            density (float): The density of the material.
            voigt_stiffness (float): The Voigt stiffness of the material.
        """

        Plotter.plot_vp(voigt_stiffness, density)

    def voigt_velocity(self):
        tensor = Tensor()
        C = tensor.tensor_to_voigt(self.cijkl)
        # Bulk modulus
        Kv = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[0,2])) / 9
        # Shear modulus
        Gv = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[0,2] + C[1,2])) / 15 + (C[3,3] + C[4,4] + C[5,5]) / 5
        # velocities
        vp = np.sqrt((Kv + (4/3)*Gv) / self.rho)
        vs = np.sqrt(Gv / self.rho)

        return vp, vs

    def reuss_velocity(self):
        tensor = Tensor()
        C = tensor.tensor_to_voigt(self.cijkl)
        S = np.linalg.inv(C)
        Kr = 1 / (S[0,0] + S[1,1] + S[2,2] + 2*(S[0,1] + S[0,2] + S[1,2]))
        Gr = 15 / (4*(S[0,0]+S[1,1]+S[2,2]-S[0,1]-S[0,2]-S[1,2]) + 3*(S[3,3]+S[4,4]+S[5,5]))
        vp = np.sqrt((Kr + (4/3)*Gr) / self.rho)
        vs = np.sqrt(Gr / self.rho)

        return vp, vs

    def hill_velocity(self):
        vp_voigt, vs_voigt = self.voigt_velocity()
        vp_reuss, vs_reuss = self.reuss_velocity()
        vp_hill = (vp_voigt + vp_reuss) / 2
        vs_hill = (vs_voigt + vs_reuss) / 2

        return vp_hill, vs_hill, vs_hill

    def get_compliance(self):
        tensor = Tensor()
        C_voigt = tensor.tensor_to_voigt(self.cijkl)
        compliance = np.linalg.inv(C_voigt)
        return compliance