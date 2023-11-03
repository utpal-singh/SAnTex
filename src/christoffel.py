import numpy as np


def christoffel_tensor(stiffness, direction):

    # Normalize the direction vector to unit length
    direction = direction / np.linalg.norm(direction)
    print("Direction: ", direction)

    christoffel = np.zeros((3, 3))

    for i in range(3):
        for k in range(3):
            # Sum over  products of the stiffness tensor and the direction vector elements
            for j in range(3):
                for l in range(3):
                    christoffel[i][k] += stiffness[i][j][k][l] * direction[j] * direction[l]


    return christoffel

def wave_properties(christoffel, density):
    eigenvalues, eigenvectors = np.linalg.eig(christoffel)
    indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[indices]
    eigenvectors = eigenvectors[:, indices]

    wave_moduli = []
    velocities = []

    for i in range(3):
        # Wave modulus equal to eigenvalue
        wave_modulus = eigenvalues[i]
        wave_moduli.append(wave_modulus)

        velocity = np.sqrt(wave_modulus / density)
        velocities.append(velocity)

    return wave_moduli, velocities, eigenvectors


if __name__=="__main__":
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                    [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                    [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                    [  0.,      0.,      0.,     65.66,    0.,      6.415],
                    [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                    [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
    import tensor_conversion
    tensor_M = tensor_conversion.voigt_to_tensor(M)
    stiffness = tensor_M
    direction = np.array([1, 0, 0])
    density = 3500.0
    christoffel = christoffel_tensor(stiffness, direction)
    wave_moduli, velocities, polarization_directions = wave_properties(christoffel, density)
    print("Christoffel tensor:")
    print(christoffel)
    print("Wave moduli:")
    print(wave_moduli)
    print("Velocities:")
    print(velocities)
    print("Polarization directions:")
    print(polarization_directions)