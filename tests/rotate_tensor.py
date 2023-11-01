# Import numpy library for matrix operations
import numpy as np
import tensor_conversion

# Define a function to generate a rotation matrix from Euler angles
def euler_to_rotation(alpha, beta, gamma):
    # Convert angles from degrees to radians
    alpha = np.radians(alpha)
    beta = np.radians(beta)
    gamma = np.radians(gamma)
    # Define the rotation matrices around x, y and z axes
    Rx = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
    # Compute the total rotation matrix by multiplying the matrices in order
    R = Rz @ Ry @ Rx
    # Return the rotation matrix
    return R

# Define a function to rotate a 3x3x3x3 stiffness tensor with Euler angles
def rotate_tensor(tensor, alpha, beta, gamma):
    # Get the rotation matrix from the Euler angles
    R = euler_to_rotation(alpha, beta, gamma)
    # Initialize an empty 3x3x3x3 tensor for the output
    output = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    # Loop through the output tensor indices and assign the rotated values
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # Loop through the input tensor indices and sum over the products with the rotation matrix elements
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    # Apply the rotation matrix to each pair of indices
                                    output[i][j][k][l] += R[i][m] * R[j][n] * R[k][o] * R[l][p] * tensor[m][n][o][p]
    # Return the output tensor
    return output

if __name__ == "__main__":
    import numpy as np
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
    
    tensor = tensor_conversion.voigt_to_tensor(M)
    print(np.array(tensor_conversion.tensor_to_voigt(tensor)))

    alpha = 30 # Rotation angle around x axis in degrees
    beta = 30 # Rotation angle around y axis in degrees
    gamma = 30 # Rotation angle around z axis in degrees
    output = rotate_tensor(tensor, alpha, beta, gamma)
    print("Output tensor:")

    print(np.array(tensor_conversion.tensor_to_voigt(output)))
