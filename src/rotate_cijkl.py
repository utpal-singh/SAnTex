# Import numpy library for array operations and math library for trigonometric functions
import numpy as np
import math

# Define a function to rotate a 3*3*3*3 stiffness tensor using the bunge convention phi1, phi and phi2 in the ZXZ format
def rotate_cijkl(cijkl, phi1, phi, phi2):
    # Convert the angles from degrees to radians
    phi1 = math.radians(phi1)
    phi = math.radians(phi)
    phi2 = math.radians(phi2)
    # Define the rotation matrices for each angle
    R1 = np.array([[math.cos(phi1), math.sin(phi1), 0],
                   [-math.sin(phi1), math.cos(phi1), 0],
                   [0, 0, 1]])
    R2 = np.array([[1, 0, 0],
                   [0, math.cos(phi), math.sin(phi)],
                   [0, -math.sin(phi), math.cos(phi)]])
    R3 = np.array([[math.cos(phi2), math.sin(phi2), 0],
                   [-math.sin(phi2), math.cos(phi2), 0],
                   [0, 0, 1]])
    # Compute the total rotation matrix by multiplying the three matrices
    R = np.matmul(np.matmul(R3, R2), R1)
    # Initialize an empty rotated stiffness tensor
    cijkl_rotated = np.zeros((3, 3, 3, 3))
    # Loop through all possible values of i, j, k and l from 0 to 2
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # Loop through all possible values of p, q, r and s from 0 to 2
                    for p in range(3):
                        for q in range(3):
                            for r in range(3):
                                for s in range(3):
                                    # Apply the rotation formula to the stiffness tensor elements
                                    cijkl_rotated[i][j][k][l] += R[i][p] * R[j][q] * R[k][r] * R[l][s] * cijkl[p][q][r][s]
    # Return the rotated stiffness tensor
    return cijkl_rotated

# Test the function with an example cijkl stiffness tensor and some angles
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
         [ -5.0 ,   -5.0 ,   50.0 ]],

        [[ -5.0 ,   -5.0 ,   50.0 ],
         [ 50.0 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ]]],


       [[[ -5.0 ,   -5.0 ,   50.0 ],
         [ 50.0 ,    1.2 ,    1.2 ],
         [ 50.00,    1.,     . ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]]]])

phi1 = 30
phi = 45
phi2 = 60

# Rotate the cijkl stiffness tensor using the bunge convention
cijkl_rotated = rotate_cijkl(cijkl, phi1, phi, phi2)
# Print the rotated stiffness tensor
print("The rotated stiffness tensor is:")
print(cijkl_rotated)
