import numpy as np

def christoffel_tensor(cijkl, n):
    tik = np.zeros((3, 3))

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    tik[i][k] += cijkl[i][j][k][l] * n[j] * n[l]
                    
    return tik

def wave_property(tik):
    eigenvalues, eigenvectors = np.linalg.eig(tik)
    indices = np.argsort(eigenvalues)[::-1]
    wave_moduli = []
    polarization_directions = []

    for i in indices:
        wave_moduli.append(eigenvalues[i])
        polarization_directions.append(eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]))
    return tuple(wave_moduli), tuple(polarization_directions)

def tensor_to_voigt(tensor):
    voigt = [[0 for i in range(6)] for j in range(6)]
    mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}

    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):

                    I = mapping[(i, j)]
                    J = mapping[(k, l)]
                    voigt[I-1][J-1] = tensor[i-1][j-1][k-1][l-1]

    return voigt


def voigt_to_tensor(voigt):
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
    return tensor