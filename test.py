# import sage

# material_instance = sage.Material(name="Example Material", properties={"density": 2.5})
# ebsd_instance = sage.EBSD(data="some_data")
# tensor_instance = sage.Tensor(values=[1, 2, 3, 4])
# anisotropy_instance = sage.Anisotropy(material=material_instance, tensor=tensor_instance)

# material_instance.display_properties()
# ebsd_instance.analyze_data()
# print(f"Tensor Magnitude: {tensor_instance.calculate_magnitude()}")
# anisotropy_instance.calculate_effect()




# from sage.material import load_material_data, get_material_properties

# # Load material data
# materials = load_material_data("path/to/your/materials_database.json")

# # Get properties of Lawsonite
# lawsonite_properties = get_material_properties(materials, phase="Lawsonite")

# # Print or use the properties as needed
# print(lawsonite_properties)

import numpy as np
from tabulate import tabulate
from sage import Material

material_instance = Material()

# Get properties for 'Almandine-pyrope'
almandine_properties = material_instance.get_properties_by_phase('Almandine-pyrope')
print("Material Properties for Almandine-pyrope:")
print(tabulate(almandine_properties.items(), headers=["Property", "Value"], tablefmt="fancy_grid"))
print("\n")

# Get voigt matrix for 'Amphobole #3 Taramite-Tschermakite1'
non_isotropic_phase = 'Amphobole #3 Taramite-Tschermakite1'
non_isotropic_voigt_matrix = material_instance.get_voigt_matrix(non_isotropic_phase)
print(f'Voigt matrix for {non_isotropic_phase}:')
print(tabulate(non_isotropic_voigt_matrix, tablefmt="pretty"))


# almandine_properties = material_instance.get_properties_by_phase('Almandine-pyrope')
# if almandine_properties:
#     density = almandine_properties['Density(g/cm3)']
#     c11 = almandine_properties['C11']
#     # Access other properties as needed
#     print(f'Density: {density}, C11: {c11}')
# else:
#     print('Phase not found in the materials database.')


# Create an instance of Material
material_instance = Material()

# availablePhases
phases_info = material_instance.availablePhases()
print("Available Phases:")
print(phases_info)
print()

# availablePrimaryPhases
primary_phases_info = material_instance.availablePrimaryPhases()
print("Available Primary Phases:")
print(primary_phases_info)
print()

# availableCrystalSystems
crystal_systems_info = material_instance.availableCrystalSystems()
print("Available Crystal Systems:")
print(crystal_systems_info)
print()


from sage import Tensor
voigt_matrix = np.array([[198.96, 73.595, 68.185, 0., 9.735, 0.],
                                [73.595, 155.94, 62.23, 0., 6.295, 0.],
                                [68.185, 62.23, 225.99, 0., 33.85, 0.],
                                [0., 0., 0., 65.66, 0., 6.415],
                                [9.735, 6.295, 33.85, 0., 60.23, 0.],
                                [0., 0., 0., 6.415, 0., 65.18]]) * 10**9
# Example usage
# tensor_data = ...  # Replace with your actual tensor data
# tensor_instance = Tensor(tensor_data)

# Convert to Voigt form
voigt_matrix = Tensor.tensor_to_voigt(tensor_data)

# Convert back to tensor
tensor_instance_back = Tensor.voigt_to_tensor(voigt_matrix)
