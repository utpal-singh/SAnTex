import numpy as np
class Tensor:
    def __init__(self, data=None):
        self.data = np.array(data) if data is not None else None

    def tensor_to_voigt(self, tensor):
        voigt = [[0 for i in range(6)] for j in range(6)]
        mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}

        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    for l in range(1, 4):

                        I = mapping[(i, j)]
                        J = mapping[(k, l)]
                        voigt[I-1][J-1] = tensor[i-1][j-1][k-1][l-1]

        return np.array(voigt)


    def voigt_to_tensor(self, voigt):
        tensor = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}
        
        for I in range(1, 7):
            for J in range(1, 7):
                i, j = mapping[I]
                k, l = mapping[J]
                tensor[i-1][j-1][k-1][l-1] = voigt[I-1][J-1]
                tensor[j-1][i-1][k-1][l-1] = voigt[I-1][J-1]
                tensor[i-1][j-1][l-1][k-1] = voigt[I-1][J-1]
                tensor[j-1][i-1][l-1][k-1] = voigt[I-1][J-1]
        return np.array(tensor)
    
# import numpy as np

    def euler_to_rotation(self, phi1, phi, phi2):
        phi1 = np.radians(phi1)
        phi = np.radians(phi)
        phi2 = np.radians(phi2)
        Rz1 = np.array([[np.cos(phi1), -np.sin(phi1), 0], [np.sin(phi1), np.cos(phi1), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
        Rz2 = np.array([[np.cos(phi2), -np.sin(phi2), 0], [np.sin(phi2), np.cos(phi2), 0], [0, 0, 1]])
        R = Rz1 @ Rx @ Rz2 

        return R

    def rotate_tensor(self, tensor, phi1, phi, phi2):
        R = self.euler_to_rotation(phi1, phi, phi2)
        tensor = np.array(tensor)
        R_tensor = np.einsum('im,jn,ko,lp,mnop->ijkl', R, R, R, R, tensor)

        return R_tensor

if __name__ == "__main__":
    x = Tensor()
    tensor = [[[[1000.0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    # print(tensor)
    voigt = x.tensor_to_voigt(tensor)
    print(voigt)
    tensor_2 = x.voigt_to_tensor(voigt)

    voigt_2 = x.tensor_to_voigt(tensor_2)

    print(tensor == tensor_2)
    print("Voigt matrix:")
    for row in voigt:
        print(row)
    tensor = x.voigt_to_tensor(voigt)
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
    
    tensor_M = x.voigt_to_tensor(M)
    print(tensor_M)
    print(M)
    print((x.tensor_to_voigt(tensor_M)))
