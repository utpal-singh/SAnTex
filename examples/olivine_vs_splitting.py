import numpy as np

from sage import Anisotropy

density = 3310
c = np.array([[323.70, 66.40, 71.60, 0.000, 0.000, 0.000],
            [66.40, 197.60, 75.60, 0.000, 0.000, 0.000],
            [71.60, 75.60, 235.10, 0.000, 0.000, 0.000],
            [0.000, 0.000, 0.000, 64.62, 0.000, 0.000],
            [0.000, 0.000, 0.000, 0.000, 78.05, 0.000],
            [0.000, 0.000, 0.000, 0.000, 0.000, 79.04]]) * 1e9
# cijkl = tensor_conversion.voigt_to_tensor(c)

# Plotter.plot_wave_velocities(c, density)
anisotropy_instance = Anisotropy(c, density)
anisotropy_instance.plotter_vs_splitting(density, c)