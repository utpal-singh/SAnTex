import numpy as np
import math

def rotate_cijkl(cijkl, phi1, phi, phi2):
    phi1 = math.radians(phi1)
    phi = math.radians(phi)
    phi2 = math.radians(phi2)
    R1 = np.array([[math.cos(phi1), math.sin(phi1), 0],
                   [-math.sin(phi1), math.cos(phi1), 0],
                   [0, 0, 1]])
    R2 = np.array([[1, 0, 0],
                   [0, math.cos(phi), math.sin(phi)],
                   [0, -math.sin(phi), math.cos(phi)]])
    R3 = np.array([[math.cos(phi2), math.sin(phi2), 0],
                   [-math.sin(phi2), math.cos(phi2), 0],
                   [0, 0, 1]])
    R = np.matmul(np.matmul(R3, R2), R1)
    cijkl_rotated = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for p in range(3):
                        for q in range(3):
                            for r in range(3):
                                for s in range(3):
                                    # Apply the rotation formula to the stiffness tensor elements
                                    cijkl_rotated[i][j][k][l] += R[i][p] * R[j][q] * R[k][r] * R[l][s] * cijkl[p][q][r][s]
    return cijkl_rotated

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

cijkl_rotated = rotate_cijkl(cijkl, phi1, phi, phi2)
print("The rotated stiffness tensor is:")
print(cijkl_rotated)
