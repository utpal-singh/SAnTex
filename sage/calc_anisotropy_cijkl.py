import numpy as np
import math

# Calculate christoffel tensor, given stiffness tensor and direction vector
def christoffel_tensor(cijkl, n):
    tik = np.zeros((3, 3))

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    tik[i][k] += cijkl[i][j][k][l] * n[j] * n[l]

    return tik

# Calculate wave moduli and polarization directions from christofel tensor
def wave_property(tik):
    eigenvalues, eigenvectors = np.linalg.eig(tik)

    # Sort eigenvalues in descending order and get corresponding indices
    indices = np.argsort(eigenvalues)[::-1]
    wave_moduli = []
    polarization_directions = []

    for i in indices:
        wave_moduli.append(eigenvalues[i])
        # Append normalized eigenvector
        polarization_directions.append(eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]))

    return tuple(wave_moduli), tuple(polarization_directions)

# Calculate vp, vs1 and vs2 over all the directions possible given a stiffness tensor and a density
def phase_velocity(cijkl, rho):
    vp = []
    vs1 = []
    vs2 = []
    #step to sample step for radians
    step = math.pi / 180

    for theta in np.arange(0, math.pi + step, step):

        for phi in np.arange(0, 2 * math.pi + step, step):
            # Calculate propagation direction vector n from theta and phi using spherical coordinates

            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])

            # Christoffel tensor for this direction
            tik = christoffel_tensor(cijkl, n)

            # Wave moduli and polarization directions for this direction
            wave_moduli, polarization_directions = wave_property(tik)

            # Calculate vp, vs1 and vs2 from the wave moduli and density using square root formula and append them to the lists
            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))

    return tuple(vp), tuple(vs1), tuple(vs2)

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
         [ 50.00,    1.,     0. ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]]]])

rho = 3500


vp, vs1, vs2 = phase_velocity(cijkl, rho)

print("The values of vp are:")
print(vp)
print("The values of vs1 are:")
print(vs1)
print("The values of vs2 are:")
print(vs2)
