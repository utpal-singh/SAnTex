import numpy as np
from ..material import Material
import pandas as pd

def calc_melt_density(weight_per, pressure, temperature):
    """
    Calculates the density and component density of a melt based on its composition.

    Parameters:
        weight_per (list or dict): List or dictionary representing the weight percentage of various oxides in the melt. 
                                   If a dictionary is provided, it should have oxide names as keys and corresponding weight percentages as values.
        pressure (float): Pressure of the melt in bars.
        temperature (float): Temperature of the melt in Kelvin.

    Returns:
        density (float): Density of the melt.
        component_density (array): Component-wise density of the melt.

    Example Usage: 
    oxide = ['sio2', 'tio2', 'al2o3', 'fe2o3', 'feo', 'mgo', 'cao', 'na2o', 'k2o', 'h2o']
    weight_per = [40, 1.84, 13.7, 2.7, 9.57, 6.67, 11.5, 2.68, 0.25, 10] or weight_per can be a dictionary like {'sio2': 40, 'tio2': 1.84, 'feo': 9.57, 'mgo': 6.67, 'cao': 11.5, 'na2o': 2.68, 'k2o': 0.25, 'h2o': 10}
    """
    if type(weight_per) == list:
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
    else:
        oxide_dict = weight_per
        oxides = ['sio2', 'tio2', 'al2o3', 'fe2o3', 'feo', 'mgo', 'cao', 'na2o', 'k2o', 'h2o']
        mol_weight = np.array([60.0855, 79.88, 101.96, 159.69, 71.85, 40.3, 56.08, 61.98, 94.2, 18])
        
        weight_per = {oxide: 0 for oxide in oxides}
        weight_per.update(oxide_dict)
        
        weight_per_list = [weight_per[oxide] for oxide in oxides]
        
        total = sum(weight_per_list)
        normalised_wt_per = np.array(weight_per_list) / total * 100
        moles = normalised_wt_per / mol_weight
        total_moles = np.sum(moles)
        mole_fraction = moles / total_moles
        molar_volume = np.array([26.86, 28.32, 37.42, 41.5, 12.68, 12.02, 16.9, 29.65, 47.28, 22.9])
        tref = np.array([1773, 1773, 1773, 1723, 1723, 1773, 1773, 1773, 1773, 1273])
        dVdT = np.array([0, 0.00724, 0.00262, 0, 0.00369, 0.00327, 0.00374, 0.00768, 0.01208, 0.0095])
        dVdP = np.array([-0.000189, -0.000231, -0.000226, -0.000253, -0.000045, 0.000027, 0.000034, -0.00024, -0.000675, -0.00032])
        temp_k = np.array([temperature] * len(oxides))
        press_b = np.array([pressure] * len(oxides))
        component_density = (mole_fraction * mol_weight) / (molar_volume + (dVdT * (temp_k - tref)) + (dVdP * (press_b - 1)))
        div_vliq = (molar_volume + (dVdT * (temp_k - tref)) + (dVdP * (press_b - 1))) * mole_fraction
        vliq = np.sum(div_vliq)
        div_x_mol_wt = mole_fraction * mol_weight
        sum_div_x_mol_wt = np.sum(div_x_mol_wt)
        density = sum_div_x_mol_wt / vliq
        
        return density, component_density

def calc_melt_density_from_excel(input_file, output_file):
    """
    Calculates melt densities for samples listed in an input Excel file and saves the results to an output Excel file.

    Parameters:
        input_file (str): Path to the input Excel file containing sample data.
        output_file (str): Path to save the output Excel file with calculated melt densities.

    Sample excel file:
        Sample Name	Pressure (bars)	Temperature (K)	    SiO2	TiO2	Al2O3	Fe2O3	FeO	    MgO	    CaO	    Na2O	K2O	    H2O
        sample_wir_01	500	            1273	        45.5	2.1	    15.8	3.4	    8.2	    4.9	    10.2	3.3	    0.6	    5.5
        sample_wir_02	500	            1273	        49.7	1.9	    14.3	2.8	7   .1	    5.6	    9.8	    4.1	    0.7	    4.0

    """
    df = pd.read_excel(input_file)
    
    results_df = pd.DataFrame(columns=['Sample Name', 'Density', 'Component Density'])

    for index, row in df.iterrows():
        sample_name = row['Sample Name']
        pressure = row['Pressure (bars)']
        temperature = row['Temperature (K)']
        
        oxide_dict = row[3:].to_dict()
        
        density, component_density = calc_melt_density(oxide_dict, pressure, temperature)
        
        results_df = results_df.append({
            'Sample Name': sample_name, 
            'Density': density, 
            'Component Density': component_density
        }, ignore_index=True)

    results_df.to_excel(output_file, index=False)



def modal_rock(rock, fraction, pressure, temperature, melt=0, weight_per = []):

    """
    Calculates the average elastic constants and density of a rock mixture.

    Parameters:
        rock (list): An array of mineral names constituting the rock mixture.
        fraction (list): An array corresponding to the fractional composition of each mineral in the rock.
        pressure (float): Pressure of the rock mixture in GPa.
        temperature (float): Temperature of the rock mixture in Kelvin.
        melt (float, optional): Fraction of melt in the rock mixture. Set to a non-zero value if there is a melt component. Default is 0.
        weight_per (list, optional): Weight percentage of oxides in the melt component, required if `melt` is non-zero. Default is [].

    Returns:
        cij_average (array): Average elastic constants of the rock mixture.
        rho_average (float): Average density of the rock mixture.

    Example Usage:
        Rock is an array of minerals and fraction is another array corresponding to the particular phase. 
        eg. rock = ["Forsterite", "Diopside", "Enstatite"], fraction = [0.2, 0.5, 0.3]
    """
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
            cij.append(material.voigt_high_PT(item, pressure, temperature))
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
        cij.append(material.voigt_high_PT(item, pressure, temperature))
        rho.append(material.load_density(item, pressure, temperature))
    rho_average = np.sum(np.multiply(fraction, rho))
    cij_average = np.zeros((6, 6))

    for i in range(len(fraction)):
        cij_average += fraction[i] * cij[i]
    return cij_average, rho_average






