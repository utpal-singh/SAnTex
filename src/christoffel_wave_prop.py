import numpy as np
import math

# n - propagation direction
def christoffel_tensor(cijkl, n):
    tik = np.zeros((3, 3))

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    # Eistein summation
                    tik[i][k] += cijkl[i][j][k][l] * n[j] * n[l]

    return tik

# wave moduli and polarization directions from a christoffel tensor
def wave_property(tik):
    eigenvalues, eigenvectors = np.linalg.eig(tik)

    # Sort eigenvalues in descending order and get the corresponding indices
    indices = np.argsort(eigenvalues)[::-1]

    wave_moduli = []
    polarization_directions = []

    for i in indices:
        wave_moduli.append(eigenvalues[i])

        # Append the normalized eigenvector to polarization directions list
        polarization_directions.append(eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]))

    return tuple(wave_moduli), tuple(polarization_directions)

# calculate vp, vs1 and vs2 over all the directions possible given a stiffness tensor and a density
def phase_velocity(cijkl, rho):

    vp = []
    vs1 = []
    vs2 = []

    step = math.pi / 180

    for theta in np.arange(0, math.pi + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):

            # propagation direction vector n from theta and phi using spherical coordinates
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            
            # christoffel tensor for this direction
            tik = christoffel_tensor(cijkl, n)
            
            # wave moduli and polarization directions for this direction
            wave_moduli, polarization_directions = wave_property(tik)
            
            # vp, vs1 and vs2 from the wave moduli and density using square root formula and append them to the lists
            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))
    return tuple(vp), tuple(vs1), tuple(vs2)

if __name__ == "__main__":
    import numpy as np
    import tensor_conversion
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                    [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                    [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                    [  0.,      0.,      0.,     65.66,    0.,      6.415],
                    [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                    [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9

    cijkl = tensor_conversion.voigt_to_tensor(M)
    rho = 3500


    vp, vs1, vs2 = phase_velocity(cijkl, rho)

    print("The values of vp are:")
    print(min(vp))
    print("The values of vs1 are:")
    print(min(vs1))
    print("The values of vs2 are:")
    print(min(vs2))
