import numpy as np

mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}
inverse_mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}

def tensor_to_voigt(tensor):
    voigt = np.zeros((6, 6))
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):
                    I = mapping[(i, j)]
                    J = mapping[(k, l)]
                    voigt[I-1][J-1] = tensor[i-1][j-1][k-1][l-1]
    return voigt

def voigt_to_tensor(voigt):
    tensor = np.zeros((3, 3, 3, 3))
    for I in range(1, 7):
        for J in range(1, 7):
            i, j = inverse_mapping[I]
            k, l = inverse_mapping[J]
            tensor[i-1][j-1][k-1][l-1] = voigt[I-1][J-1]
            tensor[j-1][i-1][k-1][l-1] = voigt[I-1][J-1]
            tensor[i-1][j-1][l-1][k-1] = voigt[I-1][J-1]
            tensor[j-1][i-1][l-1][k-1] = voigt[I-1][J-1]
    return tensor

if __name__ == "__main__":
    tensor = np.full((3, 3, 3, 3), 1000.0)
    voigt = tensor_to_voigt(tensor)
    tensor_2 = voigt_to_tensor(voigt)
    voigt_2 = tensor_to_voigt(tensor_2)

    print(np.array_equal(tensor, tensor_2))
    print("Voigt matrix:")
    print(voigt)
    tensor = voigt_to_tensor(voigt)
    print("Stiffness tensor:")
    print(tensor)

    print(np.array_equal(tensor, tensor_2))
    print(np.array_equal(voigt, voigt_2))

    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
    
    tensor_M = voigt_to_tensor(M)
    print(tensor_M)
    print(M)
    print(tensor_to_voigt(tensor_M))
