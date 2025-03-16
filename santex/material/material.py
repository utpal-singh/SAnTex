import os
import json
import numpy as np
import pandas as pd
from tabulate import tabulate
from ..isotropy import Isotropy

class Material:
    def __init__(self, database_path=os.path.join(os.path.dirname(__file__), 'data/materials_database.json'), database_path2=os.path.join(os.path.dirname(__file__), 'data/derivatives_P.json'), database_path3=os.path.join(os.path.dirname(__file__),'data/derivatives_T.json')):
        """
        Initialize a Material object with data loaded from JSON files.

        Parameters:
        - database_path (str): The file path to the main materials database JSON file. Default is 'data/materials_database.json'.
        - database_path2 (str): The file path to the pressure derivatives JSON file. Default is 'data/derivatives_P.json'.
        - database_path3 (str): The file path to the temperature derivatives JSON file. Default is 'data/derivatives_T.json'.

        Returns:
        - None
        """
        self.materials_data = self.load_materials_data(database_path)
        self.pressure_deriv = self.load_pressure_deriv(database_path2)
        self.temperature_deriv = self.load_temperature_deriv(database_path3)
        # self.database_path_mat = os.path.join(os.path.dirname(__file__), 'data/materials_database.json')


    def load_materials_data(self, database_path):

        """
        Load materials data from a JSON file.

        Parameters:
        - database_path (str): The file path to the materials database JSON file.

        Returns:
        - dict: A dictionary containing the loaded materials data.
        """
        try:
            with open(database_path, 'r') as file:
                materials_data = json.load(file)
            return materials_data
        except FileNotFoundError:
            print("Database file not found.")
            return None
        except json.JSONDecodeError:
            print("Error decoding JSON data.")
            return None


    def load_pressure_deriv(self, database_path2):
        """
        Load pressure derivatives from a JSON file.

        Parameters:
        - database_path (str): The file path to the pressure derivatives JSON file.

        Returns:
        - dict: A dictionary containing the loaded pressure derivatives data.
        """
        try:
            with open(database_path2, 'r') as file:
                pressure_data = json.load(file)
            return pressure_data
        except FileNotFoundError:
            print("Database file not found.")
            return None
        except json.JSONDecodeError:
            print("Error decoding JSON data.")
            return None


    def load_temperature_deriv(self, database_path3):

        """
        Load temperature derivatives from a JSON file.

        Parameters:
        - database_path (str): The file path to the temperature derivatives JSON file.

        Returns:
        - dict: A dictionary containing the loaded temperature derivatives data.
        """
        try:
            with open(database_path3, 'r') as file:
                temperature_data = json.load(file)
            return temperature_data
        except FileNotFoundError:
            print("Database file not found.")
            return None
        except json.JSONDecodeError:
            print("Error decoding JSON data.")
            return None

    
    def load_density(self, phase, pressure = None, temperature = None):

        """
        Load the density of a material based on the given phase.
        
        Parameters:
        - phase (str): The phase of the material.
        - pressure (float or None): The pressure in MPa (optional).
        - temperature (float or None): The temperature in Kelvin (optional).
        
        Returns:
        - float or None: The density in g/cm³ or None if data not available.
        """
        if pressure is None or temperature is None:

            try:
                for material in self.materials_data:
                    try:
                        if material['Phase'] == phase:
                            return material["Density(g/cm3)"] * 1000
                    except KeyError:
                        print("Error: Density data not found or incorrect format.")
                        return None
                print(f"Error: Phase '{phase}' not found in materials data.")
                return None
            except AttributeError:
                print("Error: Materials data not available.")
                return None
            
        else:
            isotropy = Isotropy()
            density, aks, amu = isotropy.calculate_seismic_properties(phase, temperature=temperature, pressure=pressure, return_vp_vs_vbulk=False, return_aktout=False)
            return density

    

    def get_properties_by_phase(self, phase):

        """
        Get the properties of a material based on the given phase.
        
        Parameters:
        - phase (str): The phase of the material.
        
        Returns:
        - dict or None: The properties of the material or None if not found.
        """
        for material in self.materials_data:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found
    
    def get_pressure_properties_by_phase(self, phase):
        """
        Get the pressure properties of a material based on the given phase.
        
        Parameters:
        - phase (str): The phase of the material.
        
        Returns:
        - dict or None: The pressure properties or None if not found.
        """
        for material in self.pressure_deriv:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found
    
    def get_temperature_properties_by_phase(self, phase):
        """
        Get the temperature properties of a material based on the given phase.
        
        Parameters:
        - phase (str): The phase of the material.
        
        Returns:
        - dict or None: The temperature properties or None if not found.
        """

        for material in self.temperature_deriv:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found

    def get_voigt_matrix(self, phase):
        """
        Calculate and return the Voigt matrix of a material based on the given phase.
        
        Parameters:
        - phase (str): The phase of the material.
        
        Returns:
        - numpy array or None: The Voigt matrix or None if data not available.
        """
        material_properties = self.get_properties_by_phase(phase)

        if not material_properties:
            return None
        
        c11 = material_properties.get('C11', 0)
        c22 = material_properties.get('C22', 0)
        c33 = material_properties.get('C33', 0)
        c44 = material_properties.get('C44', 0)
        c55 = material_properties.get('C55', 0)
        c66 = material_properties.get('C66', 0)
        c12 = material_properties.get('C12', 0)
        c13 = material_properties.get('C13', 0)
        c23 = material_properties.get('C23', 0)
        c14 = material_properties.get('C14', 0)
        c15 = material_properties.get('C15', 0)
        c16 = material_properties.get('C16', 0)
        c24 = material_properties.get('C24', 0)
        c25 = material_properties.get('C25', 0)
        c26 = material_properties.get('C26', 0)
        c34 = material_properties.get('C34', 0)
        c35 = material_properties.get('C35', 0)
        c36 = material_properties.get('C36', 0)
        c45 = material_properties.get('C45', 0)
        c46 = material_properties.get('C46', 0)
        c56 = material_properties.get('C56', 0)

        voigt_matrix = np.array([
            [c11, c12, c13, c14, c15, c16],
            [c12, c22, c23, c24, c25, c26],
            [c13, c23, c33, c34, c35, c36],
            [c14, c24, c34, c44, c45, c46],
            [c15, c25, c35, c45, c55, c56],
            [c16, c26, c36, c46, c56, c66]
        ])

        # Make the matrix symmetric by copying the upper triangular part to the lower triangular part
        voigt_matrix = (voigt_matrix + voigt_matrix.T)/2

        return voigt_matrix
    
    def get_voigt_matrix_pressure(self, phase):
        material_properties = self.get_pressure_properties_by_phase(phase)

        if not material_properties:
            return None

        c11 = material_properties.get('C11', 0)
        c22 = material_properties.get('C22', 0)
        c33 = material_properties.get('C33', 0)
        c44 = material_properties.get('C44', 0)
        c55 = material_properties.get('C55', 0)
        c66 = material_properties.get('C66', 0)
        c12 = material_properties.get('C12', 0)
        c13 = material_properties.get('C13', 0)
        c23 = material_properties.get('C23', 0)
        c14 = material_properties.get('C14', 0)
        c15 = material_properties.get('C15', 0)
        c16 = material_properties.get('C16', 0)
        c24 = material_properties.get('C24', 0)
        c25 = material_properties.get('C25', 0)
        c26 = material_properties.get('C26', 0)
        c34 = material_properties.get('C34', 0)
        c35 = material_properties.get('C35', 0)
        c36 = material_properties.get('C36', 0)
        c45 = material_properties.get('C45', 0)
        c46 = material_properties.get('C46', 0)
        c56 = material_properties.get('C56', 0)

        voigt_matrix = np.array([
            [c11, c12, c13, c14, c15, c16],
            [c12, c22, c23, c24, c25, c26],
            [c13, c23, c33, c34, c35, c36],
            [c14, c24, c34, c44, c45, c46],
            [c15, c25, c35, c45, c55, c56],
            [c16, c26, c36, c46, c56, c66]
        ])

        # Make the matrix symmetric by copying the upper triangular part to the lower triangular part
        voigt_matrix = (voigt_matrix + voigt_matrix.T)/2

        return voigt_matrix
    
    def get_voigt_matrix_temperature(self, phase):
        material_properties = self.get_temperature_properties_by_phase(phase)

        if not material_properties:
            return None

        c11 = material_properties.get('C11', 0)
        c22 = material_properties.get('C22', 0)
        c33 = material_properties.get('C33', 0)
        c44 = material_properties.get('C44', 0)
        c55 = material_properties.get('C55', 0)
        c66 = material_properties.get('C66', 0)
        c12 = material_properties.get('C12', 0)
        c13 = material_properties.get('C13', 0)
        c23 = material_properties.get('C23', 0)
        c14 = material_properties.get('C14', 0)
        c15 = material_properties.get('C15', 0)
        c16 = material_properties.get('C16', 0)
        c24 = material_properties.get('C24', 0)
        c25 = material_properties.get('C25', 0)
        c26 = material_properties.get('C26', 0)
        c34 = material_properties.get('C34', 0)
        c35 = material_properties.get('C35', 0)
        c36 = material_properties.get('C36', 0)
        c45 = material_properties.get('C45', 0)
        c46 = material_properties.get('C46', 0)
        c56 = material_properties.get('C56', 0)

        voigt_matrix = np.array([
            [c11, c12, c13, c14, c15, c16],
            [c12, c22, c23, c24, c25, c26],
            [c13, c23, c33, c34, c35, c36],
            [c14, c24, c34, c44, c45, c46],
            [c15, c25, c35, c45, c55, c56],
            [c16, c26, c36, c46, c56, c66]
        ])

        # Make the matrix symmetric by copying the upper triangular part to the lower triangular part
        voigt_matrix = (voigt_matrix + voigt_matrix.T)/2

        return voigt_matrix
    

    
    def available_phases(self):
        """
        Generate a formatted table listing available phases along with their crystal systems and primary phases.
        
        Returns:
        - str: A formatted table listing available phases.
        """
        headers = ["Phase", "Crystal System", "Primary Phase"]
        data = [(material['Phase'], material['Crystal System'], material['Primary Phase']) for material in self.materials_data]
        return tabulate(data, headers=headers, tablefmt="pretty")

    def available_primary_phases(self):
        """
        Generate a formatted table listing available primary phases along with their crystal systems and phases.
        
        Returns:
        - str: A formatted table listing available primary phases.
        """
        headers = ["Primary Phase", "Crystal Systems", "Phases"]
        info = {}
        for material in self.materials_data:
            primary_phase = material['Primary Phase']
            crystal_system = material['Crystal System']
            phase = material['Phase']
            if primary_phase not in info:
                info[primary_phase] = {'Crystal Systems': [], 'Phases': []}
            info[primary_phase]['Crystal Systems'].append(crystal_system)
            info[primary_phase]['Phases'].append(phase)
        data = [(key, ", ".join(info[key]['Crystal Systems']), ", ".join(info[key]['Phases'])) for key in info]
        return tabulate(data, headers=headers, tablefmt="pretty")

    def available_crystal_systems(self):
        """
        Generate a formatted table listing available crystal systems along with their associated phases.
        
        Returns:
        - str: A formatted table listing available crystal systems.
        """
        headers = ["Crystal System", "Phases"]
        crystal_systems = set(material['Crystal System'] for material in self.materials_data)
        info = {system: [] for system in crystal_systems}
        for material in self.materials_data:
            system = material['Crystal System']
            phase = material['Phase']
            info[system].append(phase)
        data = [(key, ", ".join(info[key])) for key in info]
        return tabulate(data, headers=headers, tablefmt="pretty")
    

    def voigt_high_PT(self, phase, PRESSURE=0, TEMP=300):
        """
        Perform calculations to obtain the Voigt matrix at high pressure and temperature conditions.
        
        Parameters:
        - phase (str): The phase of the material.
        - PRESSURE (float): The pressure in MPa (default 0).
        - TEMP (float): The temperature in Kelvin (default 300).
        
        Returns:
        - numpy array or None: The calculated Voigt matrix or None if data not available.
        """
        dCdT = self.get_voigt_matrix_temperature(phase)
        dCdP = self.get_voigt_matrix_pressure(phase)
        
        tensor_calc = self.get_voigt_matrix(phase) + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
        return tensor_calc
    
    def modal_rock(self, rock, fraction, pressure, temperature, melt=0):
        """
        Calculate the average Voigt matrix and density of a rock composed of multiple minerals.
        
        Parameters:
        - rock (list): An array of minerals.
        - fraction (list): Fractions corresponding to each mineral in the rock.
        - pressure (float): The pressure in MPa.
        - temperature (float): The temperature in Kelvin.
        - melt (float): Amount of melting (default 0).
        
        Returns:
        - tuple: Average Voigt matrix as a numpy array and average density as a float.

        # Rock is an array of minerals and fraction is another array corresponding to the particular phase. eg. rock = ["Forsterite", "Diopside", "Enstatite"], fraction = [0.2, 0.5, 0.3]

        """
        cij = []
        rho = []
        for item in rock:
            cij.append(self.voigt_high_PT(item, pressure, temperature))
            rho.append(self.load_density(item, pressure, temperature))
        rho_average = np.sum(np.multiply(fraction, rho))
        cij_average = np.zeros((6, 6))

        for i in range(len(fraction)):
            cij_average += fraction[i] * cij[i]
        return cij_average, rho_average
    

    def modal_rock_from_excel(self, excel_path, pressure, temperature, melt=0):
        """
        Calculate the average Voigt matrix and density of rocks from an Excel worksheet.
        
        Parameters:
        - excel_path (str): Path to the Excel worksheet.
        - pressure (float): The pressure in MPa.
        - temperature (float): The temperature in Kelvin.
        - melt (float): Amount of melting (default 0).
        
        Returns:
        - tuple: Average Voigt matrix as a numpy array and average density as a float.
        """

        df = pd.read_excel(excel_path)
        rock_data = df.values 

        cij_average = np.zeros((6, 6))
        rho_total = 0

        for rock, fraction in rock_data:
            cij = self.voigthighPT(rock, pressure, temperature)
            rho = self.load_density(rock, pressure, temperature)
            rho_total += fraction * rho
            cij_average += fraction * cij

        return cij_average, rho_total

            