# Define a function to convert 3x3x3x3 stiffness tensor to Voigt matrix notation
def tensor_to_voigt(tensor):
    # Initialize an empty 6x6 matrix
    voigt = [[0 for i in range(6)] for j in range(6)]
    # Define a mapping from tensor indices to Voigt indices
    mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}
    # Loop through the tensor indices and assign the corresponding values to the Voigt matrix
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):
                    # Get the Voigt indices from the mapping
                    I = mapping[(i, j)]
                    J = mapping[(k, l)]
                    # Assign the tensor value to the Voigt matrix
                    voigt[I-1][J-1] = tensor[i-1][j-1][k-1][l-1]
    # Return the Voigt matrix
    return voigt

# Define a function to convert Voigt matrix notation to 3x3x3x3 stiffness tensor
def voigt_to_tensor(voigt):
    # Initialize an empty 3x3x3x3 tensor
    tensor = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    # Define a mapping from Voigt indices to tensor indices
    mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}
    # Loop through the Voigt indices and assign the corresponding values to the tensor
    for I in range(1, 7):
        for J in range(1, 7):
            # Get the tensor indices from the mapping
            i, j = mapping[I]
            k, l = mapping[J]
            # Assign the Voigt value to the tensor
            tensor[i-1][j-1][k-1][l-1] = voigt[I-1][J-1]
            # Apply symmetry conditions
            tensor[j-1][i-1][k-1][l-1] = voigt[I-1][J-1]
            tensor[i-1][j-1][l-1][k-1] = voigt[I-1][J-1]
            tensor[j-1][i-1][l-1][k-1] = voigt[I-1][J-1]
    # Return the tensor
    return tensor

if __name__ == "__main__":
    import numpy as np
    # Test functions with an example stiffness tensor
    tensor = [[[[1000.0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    voigt = tensor_to_voigt(tensor)
    tensor_2 = voigt_to_tensor(voigt)

    voigt_2 = tensor_to_voigt(tensor_2)

    print(tensor == tensor_2)
    print("Voigt matrix:")
    for row in voigt:
        print(row)
    tensor = voigt_to_tensor(voigt)
    print("Stiffness tensor:")
    for i in range(3):
        print(tensor[i])

    print(tensor == tensor_2)
    print(voigt == voigt_2)

    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
    
    tensor_M = voigt_to_tensor(M)
    print(tensor_M)
    print(M)
    print(np.array(tensor_to_voigt(tensor_M)))