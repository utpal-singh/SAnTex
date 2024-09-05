import numpy as np
import math
from matplotlib import pyplot as plt
from .utils import christoffel_tensor, wave_property

def phase_velocity(cijkl, rho):
    """
    Calculate phase velocities for various propagation directions.

    Parameters:
    - cijkl (array): The fourth-rank stiffness tensor for the material.
    - rho (float): Density of the material.

    Returns:
    - tuple: Tuples containing phase velocities for P-wave (vp), S-wave with higher velocity (vs1), and S-wave with lower velocity (vs2).
    """
    vp = []
    vs1 = []
    vs2 = []
    step = math.pi / 180

    for theta in np.arange(0, math.pi + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):
            #propagation direction vector n from theta and phi using spherical coordinates
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            # print(n)
            tik = christoffel_tensor(cijkl, n)
            wave_moduli, polarization_directions = wave_property(tik)
            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))
    return tuple(vp), tuple(vs1), tuple(vs2)

def plot_vp_2d(cijkl, rho, save_plot=False, filename=None, dpi=300, cmap='RdBu'):
    """
    Plot phase velocities in a 2D stereographic projection.

    Parameters:
    - cijkl (array): The fourth-rank stiffness tensor for the material.
    - rho (float): Density of the material.
    - save_plot (bool): Whether to save the plot or not. Default is False.
    - filename (str): Name of the file to save the plot. Required if save_plot is True.
    - dpi (int): Dots per inch for saving the plot. Default is 300.
    - cmap (str or Colormap): The colormap to be used for the scatter plot. Default is 'RdBu'.

    Returns:
    - None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    step = math.pi / 180
    x = []
    y = []
    c = []

    for theta in np.arange(0, math.pi / 2 + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            tik = christoffel_tensor(cijkl, n)
            wave_moduli, polarization_directions = wave_property(tik)
            vp = math.sqrt(wave_moduli[0] / rho)
            x.append(n[0] / (1 + n[2]))
            y.append(n[1] / (1 + n[2]))
            c.append(vp)

    sc = ax.scatter(x, y, c=c, cmap=cmap)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('vp in 2d stereographic projection plots with velocity represented by colour')
    cb = plt.colorbar(sc)
    cb.set_label('vp')
    cb.ax.set_yticklabels(['{:.1f}'.format(v) for v in cb.get_ticks()])

    if save_plot:
        if filename is None:
            raise ValueError("Filename must be provided when saving the plot.")
        plt.savefig(filename, dpi=dpi)
    else:
        plt.show()


