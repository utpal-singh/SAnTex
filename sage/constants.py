import pandas as pd
import numpy as np
import math

df = pd.read_csv("single_crystal_library.csv", header=0, index_col=None)

df.set_index('Phase', inplace=True)


import numpy as np

# Elastic stiffness constants (Cij) in the order C11, C22, C33, C44, C55, C66, C12, C13, C23, C15, C25, C35, C46, C14, C16, C24, C26, C34, C36, C45, C56

def constants_help(help_info):
    if help_info == "phase":
        print("""Almandine-pyrope
Grossular
Majorite
Pyrope
a_quartz_1
a_quartz_2
a_Quartz_3
a_quartz_4
a_quartz_696C
a_quartz_700C
Calcite_1
Calcite_2
Forsterite (San Carlos)
Fayalite (0.3GPa)
Lawsonite
Orthoenstatite (MgSiO3)
Orthoenstatite (MgSiO3)
Enstatite
Bronzite (Mg0.8Fe0.2SiO3)
Ferrosilite (FeSiO3)
Biotite
Muscovite
Phlogopite
Illite-smectite
Dickite
Augite
Diopside (Di72He9Jd3Cr3Ts12)
Chrome-diopside
Jadeite
Omphacite
Coesite
Amphobole #1 Richterite1
Amphobole #2 Kataphorite1
Amphobole #3 Taramite-Tschermakite1
Amphobole #4 Hornblende-Tschermakite1
Amphobole #5 Tremolite1
Amphobole #6 Edenite1
Amphobole #7 Edenite1
Amphobole #8 Pargasite1
Amphobole #9 Pargasite1
Hornblende (#1)
Hornblende (#2)
Glaucophane
Sanidine (Or83Ab15)
Sanidine (Or89Ab11)
Orthoclase (Or93Ab7)
Albite (Or0Ab100)
An0 (Albite)
An25 (Oligoclase)
An37 (Andesine)
An48 (Andesine)
An60 (Labradorite)
An67 (Labradorite)
An78 (Bytownite)
An96 (Anorthite)
Kaolinite
Nacrite""")

constants_help("phase")
# Function to build symmetric Voigt matrix
def build_symmetric_voigt_matrix(cij, system):
    # Define Voigt matrix indices for the elastic constants
    indices = {
        "Triclinic": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), 
                      (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)],
        "Monoclinic": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), 
                       (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)],
        "Cubic": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 4)],
        "Hexagonal/ Trigonal": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 4)],
        "Tetragonal": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 4)],
        "Orthorhombic": [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (0, 4)]
    }

    # Create an empty 6x6 Voigt matrix
    voigt_matrix = np.zeros((6, 6))

    # Fill the Voigt matrix with the elastic constants using indices for the specified crystal system
    for i, j in indices[system]:
        voigt_matrix[i, j] = cij.pop(0)
        voigt_matrix[j, i] = voigt_matrix[i, j]  # Ensure symmetry

    return voigt_matrix

# Crystal systems: "triclinic", "monoclinic", "cubic", "hexagonal", "tetragonal", "orthorhombic"
# crystal_systems = ["triclinic", "monoclinic", "cubic", "hexagonal", "tetragonal", "orthorhombic"]

# # Build Voigt matrices for all crystal systems
# for system in crystal_systems:
#     cij_copy = cij.copy()  # Make a copy of the original list
#     voigt_matrix = build_symmetric_voigt_matrix(cij_copy, system)
#     print(f"Voigt matrix for {system} crystal system:")
#     print(voigt_matrix)
#     print("\n")

# Define a function to build a symmetric voigt matrix from a list of elastic stiffness
def build_voigt_matrix(cij):
  # Check the length of the list
  if len(cij) != 21:
    print("The list must have 21 elements")
    return None
  # Initialize an empty 6x6 matrix
  matrix = [[0 for i in range(6)] for j in range(6)]
  # Fill the matrix according to the crystal system
  # Cubic: C11, C12, C44
  if cij[3] != 0 and cij[4] != 0 and cij[5] != 0 and all(cij[6:]) == 0:
    matrix[0][0] = matrix[1][1] = matrix[2][2] = cij[0]
    matrix[0][1] = matrix[0][2] = matrix[1][0] = matrix[1][2] = matrix[2][0] = matrix[2][1] = cij[6]
    matrix[3][3] = matrix[4][4] = matrix[5][5] = cij[3]
  # Hexagonal: C11, C12, C13, C33, C44
  elif cij[4] != 0 and cij[5] != 0 and all(cij[7:10]) == 0 and all(cij[11:]) == 0:
    matrix[0][0] = matrix[1][1] = cij[0]
    matrix[0][1] = matrix[1][0] = cij[6]
    matrix[0][2] = matrix[1][2] = matrix[2][0] = matrix[2][1] = cij[7]
    matrix[2][2] = cij[2]
    matrix[3][3] = matrix[4][4] = cij[3]
    matrix[5][5] = (cij[0] - cij[6]) / 2
  # Tetragonal: C11, C12, C13, C16, C33, C44, C66
  elif cij[5] != 0 and all(cij[8:10]) == 0 and all(cij[12:14]) == 0 and all(cij[17:]) == 0:
    matrix[0][0] = matrix[1][1] = cij[0]
    matrix[0][1] = matrix[1][0] = cij[6]
    matrix[0][2] = matrix[1][2] = matrix[2][0] = matrix[2][1] = cij[7]
    matrix[0][5] = matrix[1][4] = matrix[4][1] = matrix[5][0] = cij[15]
    matrix[2][2] = cij[2]
    matrix[3][3] = matrix[4][4] = cij[3]
    matrix[5][5] = cij[5]
  # Orthorhombic: C11, C12, C13, C22, C23, C33, C44, C55, C66
  elif all(cij[9:15]) == 0 and all(cij[18:]) == 0:
    matrix[0][0] = cij[0]
    matrix[0][1] = matrix[1][0] = cij[6]
    matrix[0][2] = matrix[2][0] = cij[7]
    matrix[1][1] = cij[1]
    matrix[1][2] = matrix[2][1] = cij[8]
    matrix[2][2] = cij[2]
    matrix[3][3] = cij[3]
    matrix[4][4] = cij[4]
    matrix[5][5] = cij[5]
  # Monoclinic: C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66
  elif all(cij[16:18]) == 0 and all(cij[19:]) == 0:
    matrix[0][0] = cij[0]
    matrix[0][1] = matrix[1][0] = cij[6]
    matrix[0][2] = matrix[2][0] = cij[7]
    matrix[0][4] = matrix[4][0] = cij[9]
    matrix[1][1] = cij[1]
    matrix[1][2] = matrix[2][1] = cij[8]
    matrix[1][4] = matrix[4][1] = cij[10]
    matrix[2][2] = cij[2]
    matrix[2][4] = matrix[4][2] = cij[11]
    matrix[3][3] = cij[3]
    matrix[3][5] = matrix[5][3] = cij[12]
    matrix[4][4] = cij[4]
    matrix[5][5] = cij[5]
  # Triclinic: C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66
  else:
    matrix[0][0] = cij[0]
    matrix[0][1] = matrix[1][0] = cij[6]
    matrix[0][2] = matrix[2][0] = cij[7]
    matrix[0][3] = matrix[3][0] = cij[13]
    matrix[0][4] = matrix[4][0] = cij[9]
    matrix[0][5] = matrix[5][0] = cij[15]
    matrix[1][1] = cij[1]
    matrix[1][2] = matrix[2][1] = cij[8]
    matrix[1][3] = matrix[3][1] = cij[16]
    matrix[1][4] = matrix[4][1] = cij[10]
    matrix[1][5] = matrix[5][1] = cij[17]
    matrix[2][2] = cij[2]
    matrix[2][3] = matrix[3][2] = cij[18]
    matrix[2][4] = matrix[4][2] = cij[11]
    matrix[2][5] = matrix[5][2] = cij[19]
    matrix[3][3] = cij[3]
    matrix[3][4] = matrix[4][3] = cij[14]
    matrix[3][5] = matrix[5][3] = cij[12]
    matrix[4][4] = cij[4]
    matrix[4][5] = matrix[5][4] = cij[20]
    matrix[5][5] = cij[5]
  # Return the matrix
  return matrix


# print(df.loc["Almandine-pyrope", 'C11'])
# print(df.loc["Almandine-pyrope", 'Crystal System'])


# def help():
#     print(df.columns)
#     print(df.loc["Almandine-pyrope", 'C11'])
#     print(df.loc["Almandine-pyrope", 'Crystal System'])

# def voigt_cubic(c11, c12, c44):
#     voigt = np.array([
#         [c11, c12, c12, 0, 0, 0],
#         [c12, c11, c12, 0, 0, 0],
#         [c12, c12, c11, 0, 0, 0],
#         [0, 0, 0, c44, 0, 0],
#         [0, 0, 0, 0, c44, 0],
#         [0, 0, 0, 0, 0, c44]
#     ])
#     return voigt

# def voigt_hexagonal(c11, c12, c13, c33, c44):
#     voigt = np.array([
#         [c11, c12, c13, 0, 0, 0],
#         [c12, c11, c13, 0, 0, 0],
#         [c13, c13, c33, 0, 0, 0],
#         [0, 0, 0, c44, 0, 0],
#         [0, 0, 0, 0, c44, 0],
#         [0, 0, 0, 0, 0, (c11-c12)/2]
#     ])
#     return voigt

# def voigt_tetragonal(c11, c12, c13, c33, c44, c66):
#     voigt = np.array([
#         [c11, c12, c13, 0, 0, 0],
#         [c12, c11, c13, 0, 0, 0],
#         [c13, c13, c33, 0, 0, 0],
#         [0, 0, 0, c44, 0, 0],
#         [0, 0, 0, 0, c44, 0],
#         [0, 0, 0, 0, 0, c66]
#     ])
#     return voigt


# def triclinic_voigt(cij):
#     # if len(cij) != 21:
#     #     raise ValueError("The input list cij should contain 21 elastic stiffness coefficients for a triclinic crystal.")
    
#     C_voigt = np.zeros((6, 6))
    
#     # Fill the Voigt matrix
#     C_voigt[0, 0] = cij[0]
#     C_voigt[1, 1] = cij[1]
#     C_voigt[2, 2] = cij[2]
#     C_voigt[3, 3] = cij[3]
#     C_voigt[4, 4] = cij[4]
#     C_voigt[5, 5] = cij[5]
#     C_voigt[0, 1] = cij[6]
#     C_voigt[1, 0] = cij[6]
#     C_voigt[0, 2] = cij[7]
#     C_voigt[2, 0] = cij[7]
#     C_voigt[1, 2] = cij[8]
#     C_voigt[2, 1] = cij[8]
#     C_voigt[0, 3] = cij[9]
#     C_voigt[3, 0] = cij[9]
#     C_voigt[1, 4] = cij[10]
#     C_voigt[4, 1] = cij[10]
#     C_voigt[2, 5] = cij[11]
#     C_voigt[5, 2] = cij[11]
#     C_voigt[3, 4] = cij[12]
#     C_voigt[4, 3] = cij[12]
#     C_voigt[3, 5] = cij[13]
#     C_voigt[5, 3] = cij[13]
#     C_voigt[4, 5] = cij[14]
#     C_voigt[5, 4] = cij[14]
    
#     return C_voigt

# def monoclinic_voigt(cij):
#     # if len(cij) != 13:
#     #     raise ValueError("The input list cij should contain 13 elastic stiffness coefficients for a monoclinic crystal.")
    
#     C_voigt = np.zeros((6, 6))
    
#     # Fill the Voigt matrix
#     C_voigt[0, 0] = cij[0]
#     C_voigt[1, 1] = cij[1]
#     C_voigt[2, 2] = cij[2]
#     C_voigt[3, 3] = cij[3]
#     C_voigt[4, 4] = cij[4]
#     C_voigt[5, 5] = cij[5]
#     C_voigt[0, 1] = cij[6]
#     C_voigt[1, 0] = cij[6]
#     C_voigt[0, 2] = cij[7]
#     C_voigt[2, 0] = cij[7]
#     C_voigt[0, 3] = cij[8]
#     C_voigt[3, 0] = cij[8]
#     C_voigt[1, 5] = cij[9]
#     C_voigt[5, 1] = cij[9]
#     C_voigt[2, 4] = cij[10]
#     C_voigt[4, 2] = cij[10]
#     C_voigt[3, 4] = cij[11]
#     C_voigt[4, 3] = cij[11]
#     C_voigt[2, 5] = cij[12]
#     C_voigt[5, 2] = cij[12]
    
#     return C_voigt

# def orthorhombic_voigt(cij):
#     # if len(cij) != 9:
#     #     raise ValueError("The input list cij should contain 9 elastic stiffness coefficients for an orthorhombic crystal.")
    
#     C_voigt = np.zeros((6, 6))
    
#     # Fill the Voigt matrix
#     C_voigt[0, 0] = cij[0]
#     C_voigt[1, 1] = cij[1]
#     C_voigt[2, 2] = cij[2]
#     C_voigt[3, 3] = cij[3]
#     C_voigt[4, 4] = cij[4]
#     C_voigt[5, 5] = cij[5]
    
#     return C_voigt

# # Example usage:
# # triclinic_cij = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
# # monoclinic_cij = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
# # orthorhombic_cij = [1, 2, 3, 4, 5, 6, 7, 8, 9]

# # C_triclinic = triclinic_voigt(triclinic_cij)
# # C_monoclinic = monoclinic_voigt(monoclinic_cij)
# # C_orthorhombic = orthorhombic_voigt(orthorhombic_cij)

# # print("Triclinic Voigt matrix:")
# # print(C_triclinic)
# # print("\nMonoclinic Voigt matrix:")
# # print(C_monoclinic)
# # print("\nOrthorhombic Voigt matrix:")
# # print(C_orthorhombic)



def get_voigt_phase(phase, *temperature):
    col_data = df.loc[phase]
    system = df.loc[phase]["Crystal System"]
    print(phase)
    k= 3
    cij = list()
    while k<=23:
        cij_elem = col_data.iloc[k]
        cij.append(cij_elem)
        k+=1

    cij = [0 if isinstance(x, float) and math.isnan(x) else x for x in cij]
    return build_voigt_matrix(cij)


if __name__ == "__main__":
    # pass
    # get_voigt_phase("fefr")
    # help()
    # phase = "Almandine-pyrope"
    # # col_data = df.loc[phase]
    # system = df.loc[phase]["Crystal System"]
    # # print(system)
    # # k= 3
    # # cij = list()
    # # while k<=23:
    # #     cij_elem = col_data.iloc[k]
    # #     cij.append(cij_elem)
    # #     k+=1
    # # print(cij)
    # # # print(system)
    # # # cij = [1, 2, 3, 0.5, 0.5, 0.6, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.05, 0.08, 0.07, 0.06, 0.09, 0.15, 0.12, 0.11, 0.17, 0.16]

    # cij_ = df.loc[phase]

    # # # print(cij_["C11"])
    # # print(build_voigt_matrix(cij, system))
    # print(np.array(get_voigt_phase("Chrome-diopside")))
    pass