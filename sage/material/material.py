import os
import json
import numpy as np
from tabulate import tabulate

# database_path = os.path.join(os.path.dirname(__file__), "materials_database.csv")

class Material:
    def __init__(self, database_path=os.path.join(os.path.dirname(__file__), 'data/materials_database.json'), database_path2=os.path.join(os.path.dirname(__file__), 'data/derivatives_P.json'), database_path3=os.path.join(os.path.dirname(__file__),'data/derivatives_T.json')):
        self.materials_data = self.load_materials_data(database_path)
        self.pressure_deriv = self.load_pressure_deriv(database_path2)
        self.temperature_deriv = self.load_temperature_deriv(database_path3)


    def load_materials_data(self, database_path):
        with open(database_path, 'r') as file:
            materials_data = json.load(file)
        return materials_data
    

    def load_pressure_deriv(self, database_path2):
        with open(database_path2, 'r') as file:
            pressure_data = json.load(file)
        return pressure_data
    

    def load_temperature_deriv(self, database_path3):
        with open(database_path3, 'r') as file:
            temperature_data = json.load(file)
        return temperature_data
    

    def get_properties_by_phase(self, phase):
        for material in self.materials_data:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found
    
    def get_pressure_properties_by_phase(self, phase):
        for material in self.pressure_deriv:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found
    
    def get_temperature_properties_by_phase(self, phase):
        for material in self.temperature_deriv:
            if material['Phase'] == phase:
                return material
        return None  # Return None if the phase is not found

    def get_voigt_matrix(self, phase):
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
        # voigt_matrix = (voigt_matrix + voigt_matrix.T)/2 - np.diag(voigt_matrix.diagonal())
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
    

    
    def availablePhases(self):
        headers = ["Phase", "Crystal System", "Primary Phase"]
        data = [(material['Phase'], material['Crystal System'], material['Primary Phase']) for material in self.materials_data]
        return tabulate(data, headers=headers, tablefmt="pretty")

    def availablePrimaryPhases(self):
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

    def availableCrystalSystems(self):
        headers = ["Crystal System", "Phases"]
        crystal_systems = set(material['Crystal System'] for material in self.materials_data)
        info = {system: [] for system in crystal_systems}
        for material in self.materials_data:
            system = material['Crystal System']
            phase = material['Phase']
            info[system].append(phase)
        data = [(key, ", ".join(info[key])) for key in info]
        return tabulate(data, headers=headers, tablefmt="pretty")
    



    def voigthighPT(self, phase, PRESSURE=0, TEMP=300):
            dCdT = self.get_voigt_matrix_temperature(phase)
            dCdP = self.get_voigt_matrix_pressure(phase)
            print(density)
            
            tensor_calc = self.get_voigt_matrix(phase) + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
            return tensor_calc

            

if __name__ == "__main__":
    mat = Material("olivine")
    print(calculate_tensor_deriv("olivine", 2, 1000))
