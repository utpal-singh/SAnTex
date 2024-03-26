import numpy as np
from .material import Material

def modalRock(rock, fraction, pressure, temperature, melt=0):
    # Rock is an array of minerals and fraction is another array corresponding to the particular phase. eg. rock = ["Forsterite", "Diopside", "Enstatite"]
    fraction = [0.2, 0.5, 0.3]
    cij = []
    rho = []
    for item in rock:
        material = Material()
        cij.append(material.voigthighPT(item, pressure, temperature))
        rho.append(material.load_density(item, pressure, temperature))
    rho_average = np.sum(np.multiply(fraction, rho))
    cij_average = np.zeros((6, 6))

    for i in range(len(fraction)):
        cij_average += fraction[i] * cij[i]
    return cij_average, rho_average