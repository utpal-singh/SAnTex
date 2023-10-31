# Import numpy for matrix operations and linear algebra
import numpy as np

# Define a function to calculate vp, vs1 and vs2 given a stiffness tensor and a propagation direction
def calculate_velocities(tensor, n):
    # Define the Christoffel tensor as T = C * n * n
    T = np.zeros((3, 3)) # Initialize an empty 3*3 matrix
    # Loop through the tensor indices
    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    # Sum over the repeated indices j and l
                    T[i][k] += tensor[i][j][k][l] * n[j] * n[l]
    # Find the eigenvalues of the Christoffel tensor
    eigenvalues = np.linalg.eigvals(T)
    print(T)
    # Sort the eigenvalues in ascending order
    eigenvalues = np.sort(eigenvalues)
    # Divide the eigenvalues by the density to get the wave moduli
    rho = 1000 # Assume a constant density of 1000 kg/m^3
    M = eigenvalues / rho
    # Take the square root of the wave moduli to get the velocities
    velocities = np.sqrt(M)
    # Assign the velocities to vp, vs1 and vs2
    vp = velocities[2] # The largest eigenvalue corresponds to vp
    vs1 = velocities[1] # The middle eigenvalue corresponds to vs1
    vs2 = velocities[0] # The smallest eigenvalue corresponds to vs2
    # Return the velocities as a tuple
    return (vp, vs1, vs2)

# Example usage
# Define a sample stiffness tensor
M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
import tensor_conversion
tensor_M = tensor_conversion.voigt_to_tensor(M)
stiffness = tensor_M# Define a sample propagation direction
n = [0, 1, 0] # A unit vector along (1, 1, sqrt(2))
# Calculate the velocities using the function
vp, vs1, vs2 = calculate_velocities(stiffness, n)
# Print the velocities
print("vp =", vp)
print("vs1 =", vs1)
print("vs2 =", vs2)
