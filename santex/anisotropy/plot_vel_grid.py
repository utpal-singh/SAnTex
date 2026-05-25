import numpy as np
from matplotlib import pyplot as plt

from .anisotropy import Anisotropy
from ..material import Material


def plot_velocity_grid(pressure_range, temperature_range, return_type, is_ebsd=False,
                       phase=None, grid=[5, 5], filename=None, fig_name=None,
                       save_plot=False, dpi=300, *args):
    """
    Generates a grid of velocity values based on specified pressure and temperature ranges.

    Parameters:
        pressure_range (tuple): Min and max pressure values (GPa).
        temperature_range (tuple): Min and max temperature values (K).
        return_type (str): Velocity value to return — one of 'maxvp', 'minvp',
            'maxvs1', 'minvs1', 'maxvs2', 'minvs2'.
        is_ebsd (bool): Whether EBSD data is used. Default is False.
        phase (str): Phase name. Required when is_ebsd is False.
        grid (list): Grid size as [rows, columns]. Default is [5, 5].
        filename (str): EBSD data file path. Required when is_ebsd is True.
        fig_name (str): Output filename when save_plot is True.
        save_plot (bool): Save the figure instead of showing it. Default is False.
        dpi (int): Resolution for saved figure. Default is 300.
        *args: Phase indices when is_ebsd is True.

    Returns:
        np.ndarray: Grid of velocity values.
    """
    if not is_ebsd:
        pressure_grid = np.linspace(pressure_range[0], pressure_range[1], grid[0])
        temperature_grid = np.linspace(temperature_range[0], temperature_range[1], grid[1])
        vp_matrix = np.zeros((grid[0], grid[1]))

        for i, p in enumerate(pressure_grid):
            for j, t in enumerate(temperature_grid):
                material = Material()
                cij_highpt = material.voigt_high_PT(phase, p, t) * 10**9
                density_highpt = material.load_density(phase, p, t)
                anisotropy = Anisotropy(cij_highpt, density_highpt)
                vp_matrix[i, j] = anisotropy.anisotropy_values(return_values=return_type)

        pressure_mesh, temperature_mesh = np.meshgrid(pressure_grid, temperature_grid)
        plt.figure(figsize=(8, 6))
        plt.pcolormesh(pressure_mesh, temperature_mesh, vp_matrix, shading='auto', cmap='viridis')
        plt.colorbar(label=f'{return_type}')
        plt.xlabel('Pressure')
        plt.ylabel('Temperature')
        plt.title(f'{return_type} vs. Pressure and Temperature')
        plt.tight_layout()
        if save_plot:
            if fig_name is None:
                raise ValueError("fig_name must be provided when save_plot is True.")
            plt.savefig(fig_name, dpi=dpi)
        else:
            plt.show()

        return vp_matrix

    else:
        # EBSD path — import here to avoid a cross-package dependency at module load time
        from ..ebsd import EBSD

        if filename is None:
            raise ValueError("filename must be provided when is_ebsd is True.")
        if phase is None:
            raise ValueError("phase must be provided when is_ebsd is True.")

        ebsd = EBSD(filename)
        df = ebsd.get_ebsd_data()
        euler_angles = [ebsd.get_euler_angles(phase=i, data=df) for i in args]

        average_tensor, average_density = ebsd.get_anisotropy_for_ebsd(
            cij=None, euler_angles=euler_angles, density=None
        )
        anis = Anisotropy(average_tensor * 10**9, average_density)
        return anis.anisotropy_values(return_values=return_type)
