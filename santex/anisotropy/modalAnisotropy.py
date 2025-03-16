from ..material import Material
from ..anisotropy import Anisotropy

import numpy as np

def modal_anisotropy(material, fraction, pressure, temperature):
    voigtMatrixTotal = []
    density = []
    for i in range(len(material)):
        material_instance = Material()
        phase = 'Diopside'
        voigtMatrix = material_instance.voigt_high_PT(material[i], PRESSURE = pressure, TEMP = temperature)
        voigtMatrixTotal.append(voigtMatrix)
        density.append(material_instance.load_density(material[i]))

    aggregateStiffnessTensor = np.prod(np.array(voigtMatrixTotal), np.array(fraction))
    density_aggregate = np.prod(np.array(density), np.array(fraction))
    
    anisotropy_instance = Anisotropy(stiffness_matrix=aggregateStiffnessTensor, density=density_aggregate)



