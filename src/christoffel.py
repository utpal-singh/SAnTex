
import numpy as np

# Define a function to calculate the Christoffel tensor from the stiffness tensor and the propagation direction
def christoffel_tensor(stiffness, direction):
    # Normalize the direction vector to unit length
    direction = direction / np.linalg.norm(direction)
    print("Direction: ", direction)
    # Initialize an empty 3x3 matrix for the Christoffel tensor
    christoffel = np.zeros((3, 3))
    # Loop through the Christoffel tensor indices and assign the values
    for i in range(3):
        for k in range(3):
            # Sum over the products of the stiffness tensor and the direction vector elements
            for j in range(3):
                for l in range(3):
                    christoffel[i][k] += stiffness[i][j][k][l] * direction[j] * direction[l]
    # Return the Christoffel tensor
    return christoffel

# Define a function to calculate the wave moduli and polarization directions from the Christoffel tensor and the density
def wave_properties(christoffel, density):
    # Compute the eigenvalues and eigenvectors of the Christoffel tensor
    eigenvalues, eigenvectors = np.linalg.eig(christoffel)
    # Sort the eigenvalues and eigenvectors in descending order
    indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[indices]
    eigenvectors = eigenvectors[:, indices]
    # Initialize empty lists for the wave moduli and velocities
    wave_moduli = []
    velocities = []
    # Loop through the eigenvalues and compute the wave moduli and velocities
    for i in range(3):
        # Wave modulus is equal to eigenvalue
        wave_modulus = eigenvalues[i]
        wave_moduli.append(wave_modulus)
        # Velocity is equal to square root of wave modulus divided by density
        velocity = np.sqrt(wave_modulus / density)
        velocities.append(velocity)
    # Return the wave moduli, velocities and polarization directions
    return wave_moduli, velocities, eigenvectors

# Test the functions with an example stiffness tensor, propagation direction and density

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