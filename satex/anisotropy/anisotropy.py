import numpy as np
import math
from ..tensor import Tensor

import matplotlib.pyplot as plt
import numpy as np
from .vtkplotter import Plotter
import plotly.graph_objects as go
from plotly.subplots import make_subplots



class Anisotropy:
    def __init__(self, stiffness_matrix, density):
        try:
            tensor_object = Tensor()
            self.cijkl = tensor_object.voigt_to_tensor(stiffness_matrix)
            self.rho = density
        except Exception as e:
            raise ValueError("Error in initializing Anisotropy object: " + str(e))

    def christoffel_tensor(self, n):
        try:
            tik = np.zeros((3, 3))

            for i in range(3):
                for k in range(3):
                    tik[i, k] = np.tensordot(self.cijkl[i, :, k, :], np.outer(n, n))

            return tik
        except Exception as e:
            raise ValueError("Error in calculating the Christoffel tensor:", e)

    def wave_property(self, tik):
        try:
            eigenvalues, eigenvectors = np.linalg.eig(tik)
            indices = np.argsort(eigenvalues)[::-1]

            wave_moduli = [eigenvalues[i] for i in indices]
            polarization_directions = [eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]) for i in indices]

            return tuple(wave_moduli), tuple(polarization_directions)
        except Exception as e:
            raise ValueError("Error in calculating wave properties:", e)

    def phase_velocity(self):
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


    def plot(self):
        try:
            fig, axs = plt.subplots(2, 3, figsize=(15, 10))
            step = math.pi / 180

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

                sc = ax.scatter(x, y, c=c, cmap='RdBu')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect('equal', 'box')
                
                # Adding text at the bottom of each plot
                text = 'p'  # You can customize this text
                ax.text(0.5, -0.15, text, ha='center', transform=ax.transAxes)

            axs[0, 0].set_title('VP/VS1')
            axs[0, 1].set_title('VP')
            axs[0, 2].set_title('VS1')
            axs[1, 0].set_title('VS2')
            axs[1, 1].set_title('AVpVs1')
            axs[1, 2].set_title('AVpVs2')

            plt.tight_layout()
            plt.show()
        except Exception as e:
            raise ValueError("Error in plotting:", e)
        

    def plotvs1(self):
        try:
            fig, ax = plt.subplots(figsize=(8, 6))

            step = math.pi / 180

            x = []
            y = []
            c = []

            for theta in np.arange(0, math.pi / 2 + step, step):
                for phi in np.arange(0, 2 * math.pi + step, step):
                    n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
                    tik = self.christoffel_tensor(n)
                    if tik is not None:
                        wave_moduli, _ = self.wave_property(tik)
                        vs1 = math.sqrt(wave_moduli[1] / self.rho)
                        x.append(n[0] / (1 + n[2]))
                        y.append(n[1] / (1 + n[2]))
                        c.append(vs1)

            sc = ax.scatter(x, y, c=c, cmap='RdBu')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_aspect('equal', 'box')

            # Adding text at the bottom of the plot
            text = 'Some text'  # You can customize this text
            ax.text(0.5, -0.15, text, ha='center', transform=ax.transAxes)

            ax.set_title('VS1')

            plt.colorbar(sc, ax=ax, label='VS1')
            plt.tight_layout()
            plt.show()
        except Exception as e:
            raise ValueError("Error in plotting:", e)

    def plotly(self):
        try:
            fig = go.Figure()

            step = math.pi / 180

            # Create subplots
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

                # Add trace to subplot
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



if __name__ == "__main__":
    try:
        stiffness_matrix = np.array([[198.96, 73.595, 68.185, 0., 9.735, 0.],
                                    [73.595, 155.94, 62.23, 0., 6.295, 0.],
                                    [68.185, 62.23, 225.99, 0., 33.85, 0.],
                                    [0., 0., 0., 65.66, 0., 6.415],
                                    [9.735, 6.295, 33.85, 0., 60.23, 0.],
                                    [0., 0., 0., 6.415, 0., 65.18]]) * 10**9

        density = 3500

        anisotropy_instance = Anisotropy(stiffness_matrix, density)

        vp, vs1, vs2 = anisotropy_instance.phase_velocity()

        print("The values of vp are:")
        print(np.min(vp))
        print("The values of vs1 are:")
        print(np.min(vs1))
        print("The values of vs2 are:")
        print(np.min(vs2))
    except ValueError as e:
        print("Error:", e)
