import numpy as np

def build_voigt_matrix(cij, crystal_system):
  """Builds a symmetric Voigt matrix from a list of elastic stiffness coefficients.

  Args:
    cij: A list of 21 elastic stiffness coefficients, in the order
      [C11, C22, C33, C44, C55, C66, C12, C13, C23, C15, C25, C35, C46, C14, C16,
       C24, C26, C34, C36, C45, C56].
    crystal_system: The crystal system of the material, either 'Triclinic',
      'Monoclinic', 'Cubic', 'Hexagonal', 'Tetragonal', or 'Orthorhombic'.

  Returns:
    A numpy array representing the symmetric Voigt matrix.
  """

  # Create a zero-filled Voigt matrix.
  voigt_matrix = np.zeros([6, 6])

  # Fill in the Voigt matrix according to the crystal system.
  if crystal_system == 'Triclinic':
    # All 21 elastic stiffness coefficients are present.
    voigt_matrix[:] = cij
  elif crystal_system == 'Monoclinic':
    # The elastic stiffness coefficients C12, C13, and C23 are equal.
    voigt_matrix[0, 1] = voigt_matrix[1, 0] = cij[6]
    voigt_matrix[0, 2] = voigt_matrix[2, 0] = cij[7]
    voigt_matrix[1, 2] = voigt_matrix[2, 1] = cij[8]
    voigt_matrix[:] = cij
  elif crystal_system == 'Cubic':
    # The elastic stiffness coefficients C11, C12, and C44 are equal.
    voigt_matrix[0, 1] = voigt_matrix[1, 0] = cij[6]
    voigt_matrix[4, 4] = cij[3]
    voigt_matrix[:] = cij
  elif crystal_system == 'Hexagonal':
    # The elastic stiffness coefficients C11, C12, and C33 are equal, and
    # the elastic stiffness coefficients C13 and C44 are equal.
    voigt_matrix[0, 1] = voigt_matrix[1, 0] = cij[6]
    voigt_matrix[0, 2] = voigt_matrix[2, 0] = cij[10]
    voigt_matrix[4, 4] = cij[3]
    voigt_matrix[:] = cij
  elif crystal_system == 'Tetragonal':
    # The elastic stiffness coefficients C11, C12, and C66 are equal, and
    # the elastic stiffness coefficients C13 and C44 are equal.
    voigt_matrix[0, 1] = voigt_matrix[1, 0] = cij[6]
    voigt_matrix[4, 4] = cij[3]
    voigt_matrix[:] = cij
  elif crystal_system == 'Orthorhombic':
    # The elastic stiffness coefficients C11, C22, and C33 are distinct.
    voigt_matrix[:] = cij

  # Make the Voigt matrix symmetric.
  voigt_matrix = np.triu(voigt_matrix) + np.tril(voigt_matrix.T)

  return voigt_matrix


# Example usage:

cij = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]

# Build the Voigt matrix for a triclinic crystal system.
voigt_matrix_triclinic = build_voigt_matrix(cij, crystal_system='Triclinic')

# Build the Voigt matrix for a monoclinic crystal system.
voigt_matrix_monoclinic = build_voigt_matrix(cij, crystal_system)
