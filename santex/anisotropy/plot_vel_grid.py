import numpy as np
from matplotlib import pyplot as plt

def plot_velocity_grid(pressure_range, temperature_range, return_type, is_ebsd = False, phase = None, grid = [5, 5], filename = None,fig_name=None, save_plot=False, dpi=300, *args):
    from santex.anisotropy import Anisotropy
    from santex.ebsd import EBSD
    from santex.material import Material
    """
    Generates a grid of velocity values based on specified pressure and temperature ranges.

    Parameters:
        pressure_range (tuple): A tuple containing the minimum and maximum pressure values (in GPa).
        temperature_range (tuple): A tuple containing the minimum and maximum temperature values (in Kelvin).
        return_type (str): The type of velocity values to return. Allowed values are 'maxvp', 'minv', 'maxvs1', 'minvs1', 'maxvs2', 'minvs2'.
        is_ebsd (bool, optional): Indicates whether EBSD data is used. Default is False.
        phase (str, optional): The phase for which velocity values are calculated. Required if is_ebsd is False.
        grid (list, optional): A list specifying the grid size as [rows, columns]. Default is [5, 5].
        filename (str, optional): The filename of the EBSD data file. Required if is_ebsd is True.
        *args: Additional phase indices if EBSD data is used.

    Returns:
        np.ndarray: A grid of velocity values corresponding to the specified pressure and temperature ranges.

    Notes:
        - If is_ebsd is False, phase must be provided to calculate velocity values for the specified phase.
        - If is_ebsd is True, filename must be provided along with additional phase indices (*args) to calculate velocity values from EBSD data.
        - The return_type parameter determines the type of velocity values returned, such as maximum velocity (maxvp), minimum velocity (minv), etc.
    
    Usage:
    args can be [0, 1, 2, 3]
    give filename is is_ebsd is True
    """

    if is_ebsd is False:
        cij = []
        density = []

        pressure_grid = np.linspace(pressure_range[0], pressure_range[1], grid[0])
        temperature_grid = np.linspace(temperature_range[0], temperature_range[1], grid[1])
        vp_matrix = np.zeros((grid[0], grid[1]))
        
        for i, p in enumerate(pressure_grid):
            for j, t in enumerate(temperature_grid):
                material = Material()
                cij_highpt = material.voigt_high_PT(phase, p, t)*10**9
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
        # plt.grid(True)
        plt.tight_layout()
        if save_plot:
            if fig_name is None:
                raise ValueError("Filename must be provided when saving the plot.")
            plt.savefig(fig_name, dpi=dpi)
        else:
            plt.show()

        return vp_matrix

    if is_ebsd is True:
        euler_angles = []
        df = ebsd.get_ebsd_data()
        for i in args:
            ebsd = EBSD(filename)
            euler_angles.append(ebsd.get_euler_angles(phase=i, data = df))

        average_tensor, average_density = ebsd.getAnisotropyForEBSD(cij, euler_angles, density)
        anis = Anisotropy(average_tensor*10**9, average_density)
    



