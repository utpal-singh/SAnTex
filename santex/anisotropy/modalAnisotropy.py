from ..material import Material
from .anisotropy import Anisotropy

import numpy as np

def modal_anisotropy(material, fraction, pressure, temperature):
    voigtMatrixTotal = []
    density = []
    for i in range(len(material)):
        material_instance = Material()
        voigtMatrix = material_instance.voigt_high_PT(material[i], PRESSURE = pressure, TEMP = temperature)
        voigtMatrixTotal.append(voigtMatrix)
        density.append(material_instance.load_density(material[i]))

    fractions = np.array(fraction)
    aggregateStiffnessTensor = np.sum(
        np.array(voigtMatrixTotal) * fractions[:, np.newaxis, np.newaxis], axis=0)
    density_aggregate = np.dot(np.array(density), fractions)

    anisotropy_instance = Anisotropy(stiffness_matrix=aggregateStiffnessTensor, density=density_aggregate)



