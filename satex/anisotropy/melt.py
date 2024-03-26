import numpy as np
from satex import Material

def calc_melt_density(weight_per, pressure, temperature):
    """
    oxide = ['sio2', 'tio2', 'al2o3', 'fe2o3', 'feo', 'mgo', 'cao', 'na2o', 'k2o', 'h2o']
    weight_per = [40, 1.84, 13.7, 2.7, 9.57, 6.67, 11.5, 2.68, 0.25, 10]
    """
    oxide = ['sio2', 'tio2', 'al2o3', 'fe2o3', 'feo', 'mgo', 'cao', 'na2o', 'k2o', 'h2o']
    total = np.sum(weight_per)
    normalised_wt_per = np.array(weight_per)/total*100
    mol_weight = [60.0855, 79.88, 101.96, 159.69, 71.85, 40.3, 56.08, 61.98, 94.2, 18]
    moles = normalised_wt_per/mol_weight
    total_moles = np.sum(moles)
    mole_fraction = moles/total_moles
    molar_volume = np.array([26.86, 28.32, 37.42, 41.5, 12.68, 12.02, 16.9, 29.65, 47.28, 22.9])
    tref = np.array([1773, 1773, 1773, 1723, 1723, 1773, 1773, 1773, 1773, 1273])
    dVdT = np.array([0, 0.00724, 0.00262, 0, 0.00369, 0.00327, 0.00374, 0.00768, 0.01208, 0.0095])
    dVdP = np.array([-0.000189, -0.000231, -0.000226, -0.000253, -0.000045, 0.000027, 0.000034, -0.00024, -0.000675, -0.00032])
    temp_k = np.array([temperature]*len(oxide))
    press_b = np.array([pressure]*len(oxide))
    component_density = (mole_fraction*mol_weight)/(molar_volume+(dVdT*(temp_k-tref))+(dVdP*(press_b-1)))
    div_vliq = (molar_volume+(dVdT*(temp_k-tref))+(dVdP*(press_b-1)))*mole_fraction
    vliq = np.sum(div_vliq)
    div_x_mol_wt = mole_fraction*mol_weight
    sum_div_x_mol_wt = np.sum(div_x_mol_wt)
    density = sum_div_x_mol_wt/vliq
    return density, component_density

def modalRock(rock, fraction, pressure, temperature, melt=0, weight_per = []):
    # Rock is an array of minerals and fraction is another array corresponding to the particular phase. eg. rock = ["Forsterite", "Diopside", "Enstatite"], fraction = [0.2, 0.5, 0.3]
    if melt:
        density, comp_density = calc_melt_density(weight_per, pressure*10000, temperature)
        voigt_melt = np.array([[16.1,   15.967,  15.967,   0.,      0,   0.   ],
                                [ 15.967, 16.1,   15.967,    0.,      0,   0.   ],
                                [ 15.967,  15.967,  16.1,    0.,     0,    0.   ],
                                [  0.,      0.,      0.,     0.01,    0.,      0],
                                [  0,   0,  0,    0.,     0.01,    0.   ],
                                [  0.,      0.,     0.,      0,   0.,     0.01 ]])
        cij = []
        rho = []
        material = Material()
        for item in rock:
            cij.append(material.voigthighPT(item, pressure, temperature))
            rho.append(material.load_density(item, pressure, temperature))
        rho_average = np.sum(np.multiply(fraction, rho))
        cij_average = np.zeros((6, 6))

        for i in range(len(fraction)):
            cij_average += fraction[i] * cij[i]

        cij_average_melt = ((cij_average*(1-melt))+(melt*(voigt_melt)))
        density_average_melt = rho_average*(1-melt) + density*melt
        return cij_average_melt, density_average_melt
    cij = []
    rho = []
    material = Material()
    for item in rock:
        cij.append(material.voigthighPT(item, pressure, temperature))
        rho.append(material.load_density(item, pressure, temperature))
    rho_average = np.sum(np.multiply(fraction, rho))
    cij_average = np.zeros((6, 6))

    for i in range(len(fraction)):
        cij_average += fraction[i] * cij[i]
    return cij_average, rho_average

if __name__ == "__main__":
    weight_per = [40, 1.84, 13.7, 2.7, 9.57, 6.67, 11.5, 2.68, 0.25, 10]
    dens, component_density = calc_melt_density(weight_per, 500, 1273)
    print(dens)
    print(component_density)
    rock = ["Forsterite", "Diopside", "Enstatite"]
    fraction = [0.2, 0.5, 0.3]
    cij_average, rho_average = modalRock(rock, fraction, 2, 1000)
    print(cij_average)
    print(rho_average)

    cij_average, rho_average = modalRock(rock, fraction, 2, 1000, 0.05, weight_per)
    print(cij_average)
    print(rho_average)





