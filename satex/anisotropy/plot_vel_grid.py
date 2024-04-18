import numpy as np
from satex import Material
import matplotlib.pyplot as plt

def plot_velocity_grid(pressure_range, temperature_range, return_type, is_ebsd = False, phase = None, grid = [5, 5], filename = None, *args):
    from satex import Anisotropy
    """
    Return values can ve maxvp, minv, maxvs1, minvs1, maxvs2, minvs2 
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
                cij_highpt = material.voigthighPT(phase, p, t)*10**9
                density_highpt = material.load_density(phase, p, t)
                anisotropy = Anisotropy(cij_highpt, density_highpt)
                vp_matrix[i, j] = anisotropy.anisotropy_values(return_values=return_type)

        return vp_matrix

    if is_ebsd is True:
        euler_angles = []
        df = ebsd.get_ebsd_data()
        for i in args:
            ebsd = EBSD(filename)
            euler_angles.append(ebsd.get_euler_angles(phase=i, data = df))

        average_tensor, average_density = ebsd.getAnisotropyForEBSD(cij, euler_angles, density)
        anis = Anisotropy(average_tensor*10**9, average_density)


