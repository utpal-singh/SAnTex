import numpy as np

def euler_to_rotation(phi1, phi, phi2):

    phi1 = np.radians(phi1)
    phi = np.radians(phi)
    phi2 = np.radians(phi2)

    Rz1 = np.array([[np.cos(phi1), -np.sin(phi1), 0], [np.sin(phi1), np.cos(phi1), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    Rz2 = np.array([[np.cos(phi2), -np.sin(phi2), 0], [np.sin(phi2), np.cos(phi2), 0], [0, 0, 1]])

    R = Rz2 @ Rx @ Rz1

    return R

def rotate_tensor(tensor, phi1, phi, phi2):

    R = euler_to_rotation(phi1, phi, phi2)

    output = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):

                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    # Apply the rotation matrix to each pair of indices
                                    output[i][j][k][l] += R[i][m] * R[j][n] * R[k][o] * R[l][p] * tensor[m][n][o][p]

    return output

if __name__ == "__main__":
    import numpy as np
    import tensor_conversion
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
    
    tensor = tensor_conversion.voigt_to_tensor(M)
    print(np.array(tensor_conversion.tensor_to_voigt(tensor)))

    alpha = 60 # Rotation angle around x axis in degrees
    beta = 40 # Rotation angle around y axis in degrees
    gamma = 30    # Rotation angle around z axis in degrees
    output = rotate_tensor(tensor, alpha, beta, gamma)
    print("Output tensor:")

    print(np.array(tensor_conversion.tensor_to_voigt(output)))