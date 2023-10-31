# Import numpy library for array operations and math library for trigonometric functions
import numpy as np
import math

# Define a function to calculate the christoffel tensor given a stiffness tensor and a propagation direction
def christoffel_tensor(cijkl, n):
    # Initialize an empty 3*3 christoffel tensor
    tik = np.zeros((3, 3))
    # Loop through all possible values of i, k, j and l from 0 to 2
    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    # Apply the formula to the christoffel tensor elements
                    tik[i][k] += cijkl[i][j][k][l] * n[j] * n[l]
    # Return the christoffel tensor
    return tik

# Define a function to calculate the wave moduli and polarization directions given a christoffel tensor
def wave_property(tik):
    # Compute the eigenvalues and eigenvectors of the christoffel tensor
    eigenvalues, eigenvectors = np.linalg.eig(tik)
    # Sort the eigenvalues in descending order and get the corresponding indices
    indices = np.argsort(eigenvalues)[::-1]
    # Initialize empty lists for wave moduli and polarization directions
    wave_moduli = []
    polarization_directions = []
    # Loop through the indices
    for i in indices:
        # Append the eigenvalue to the wave moduli list
        wave_moduli.append(eigenvalues[i])
        # Append the normalized eigenvector to the polarization directions list
        polarization_directions.append(eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]))
    # Return the wave moduli and polarization directions as tuples
    return tuple(wave_moduli), tuple(polarization_directions)

# Define a function to calculate vp, vs1 and vs2 over all the directions possible given a stiffness tensor and a density
def phase_velocity(cijkl, rho):
    # Initialize empty lists for vp, vs1 and vs2
    vp = []
    vs1 = []
    vs2 = []
    # Define a step size for sampling the angles in radians
    step = math.pi / 180
    # Loop through all possible values of theta from 0 to pi with the step size
    for theta in np.arange(0, math.pi + step, step):
        # Loop through all possible values of phi from 0 to 2*pi with the step size
        for phi in np.arange(0, 2 * math.pi + step, step):
            # Calculate the propagation direction vector n from theta and phi using spherical coordinates
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            # Calculate the christoffel tensor for this direction
            tik = christoffel_tensor(cijkl, n)
            # Calculate the wave moduli and polarization directions for this direction
            wave_moduli, polarization_directions = wave_property(tik)
            # Calculate vp, vs1 and vs2 from the wave moduli and density using square root formula and append them to the lists
            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))
    # Return vp, vs1 and vs2 as tuples
    return tuple(vp), tuple(vs1), tuple(vs2)

# Test the functions with an example cijkl stiffness tensor and a density value
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
         [ 50.00,    1.,     . ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]],

        [[ .00,     .00,     .00 ],
         [.00,     .00,     .00 ],
         [.00,     .00,     .00 ]]]])

rho = 10

# Calculate vp, vs1 and vs2 over all the directions possible
vp, vs1, vs2 = phase_velocity(cijkl, rho)
# Print vp, vs1 and vs2
print("The values of vp are:")
print(vp)
print("The values of vs1 are:")
print(vs1)
print("The values of vs2 are:")
print(vs2)
