# Import numpy library for array operations
import numpy as np

# Define a function to convert a 3*3*3*3 cijkl stiffness tensor to 6*6 Voigt matrix
def cijkl_to_voigt(cijkl):
    # Initialize an empty 6*6 Voigt matrix
    voigt = np.zeros((6, 6))
    # Define a mapping from ij and kl pairs to I and J indices
    mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}
    # Loop through all possible values of i, j, k and l from 1 to 3
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):
                    # Get the corresponding values of I and J from the mapping
                    I = mapping[(i, j)]
                    J = mapping[(k, l)]
                    # Assign the value of cijkl to the Voigt matrix element at row I and column J
                    voigt[I-1][J-1] = cijkl[i-1][j-1][k-1][l-1]
    # Return the Voigt matrix
    return voigt

# Define a function to convert a 6*6 Voigt matrix to a 3*3*3*3 cijkl stiffness tensor
def voigt_to_cijkl(voigt):
    # Initialize an empty 3*3*3*3 cijkl stiffness tensor
    cijkl = np.zeros((3, 3, 3, 3))
    # Define a mapping from I and J indices to ij and kl pairs
    mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}
    # Loop through all possible values of I and J from 1 to 6
    for I in range(1,7):
        for J in range(1,7):
            # Get the corresponding values of ij and kl from the mapping
            ij = mapping[I]
            kl = mapping[J]
            # Assign the value of the Voigt matrix element at row I and column J to cijkl
            cijkl[ij[0]-1][ij[1]-1][kl[0]-1][kl[1]-1] = voigt[I-1][J-1]
            # If ij is not equal to ji or kl is not equal to lk, assign the same value to the symmetric elements of cijkl
            if ij != ij[::-1]:
                cijkl[ij[1]-1][ij[0]-1][kl[0]-1][kl[1]-1] = voigt[I-1][J-1]
            if kl != kl[::-1]:
                cijkl[ij[0]-1][ij[1]-1][kl[1]-1][kl[0]-1] = voigt[I-1][J-1]
            if ij != ij[::-1] and kl != kl[::-1]:
                cijkl[ij[1]-1][ij[0]-1][kl[0]-0][kl[0]-0] = voigt[I-0][J-0]
    # Return the cijkl stiffness tensor
    return cijkl

# Test the functions with an example cijkl stiffness tensor
cijkl = np.array([[[[100.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]]],


       [[[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   50.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   50.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ]]],


       [[[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ],
         [ -5.0 ,   50.0 ,    1.2 ]],

        [[ -5.0 ,   50.0 ,    1.2 ],
         [ 50.0 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ]],

        [[ 50.0 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ]]]])

# Convert the cijkl stiffness tensor to a Voigt matrix
voigt = cijkl_to_voigt(cijkl)
# Print the Voigt matrix
print("The Voigt matrix is:")
print(voigt)

# Convert the Voigt matrix back to a cijkl stiffness tensor
cijkl = voigt_to_cijkl(voigt