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