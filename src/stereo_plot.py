# Import numpy for matrix operations and linear algebra
import numpy as np
# Import matplotlib for plotting
import matplotlib.pyplot as plt

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
    # Sort the eigenvalues in ascending order
    eigenvalues = np.sort(eigenvalues)
    # Divide the eigenvalues by the density to get the wave moduli
    rho = 3500 # Assume a constant density of 1000 kg/m^3
    M = eigenvalues / rho
    # Take the square root of the wave moduli to get the velocities
    velocities = np.sqrt(M)
    # Assign the velocities to vp, vs1 and vs2
    vp = velocities[2] # The largest eigenvalue corresponds to vp
    vs1 = velocities[1] # The middle eigenvalue corresponds to vs1
    vs2 = velocities[0] # The smallest eigenvalue corresponds to vs2
    # Return the velocities as a tuple
    return (vp, vs1, vs2)

# Define a sample stiffness tensor
M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9
import tensor_conversion
tensor_M = tensor_conversion.voigt_to_tensor(M)
stiffness = tensor_M
# Define a grid of propagation directions on a 2d circle using polar coordinates
theta = np.linspace(0, 2*np.pi, 50) # Angle from x-axis
n_grid = np.array([np.cos(theta), np.sin(theta), np.zeros_like(theta)]) # Convert to cartesian coordinates
print(n_grid)

# Calculate the velocities for each direction on the grid using the function
vp_grid = np.zeros_like(theta) # Initialize an empty grid for vp
vs1_grid = np.zeros_like(theta) # Initialize an empty grid for vs1
vs2_grid = np.zeros_like(theta) # Initialize an empty grid for vs2

# Loop through the grid indices
for i in range(len(theta)):
    # Get the direction vector at each grid point
    n = n_grid[:, i]
    # Calculate the velocities using the function
    vp, vs1, vs2 = calculate_velocities(tensor_M, n)
    # Assign the velocities to the grids
    vp_grid[i] = vp
    vs1_grid[i] = vs1
    vs2_grid[i] = vs2

# Plot the velocities on a stereo plot of 2d using matplotlib

# Create a figure and an axis object
fig = plt.figure()
ax = fig.add_subplot(111)
print("wrjbvre")
print(max(vp_grid))
print(min(vp_grid))
print(vp_grid.shape)
print("vereve")

# Plot the vp grid as a scatter plot with a red color map and a black edge color
ax.scatter(n_grid[0], n_grid[1], c=plt.cm.Reds(vp_grid), edgecolors='k')
# Plot the vs1 grid as a scatter plot with a green color map and a black edge color
# ax.scatter(n_grid[0], n_grid[1], c=plt.cm.Greens(vs1_grid), edgecolors='k')
# # Plot the vs2 grid as a scatter plot with a blue color map and a black edge color
# ax.scatter(n_grid[0], n_grid[1], c=plt.cm.Blues(vs2_grid), edgecolors='k')

# Set the axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')

# Set the title
ax.set_title('Velocities on a 2d circle')

# Show the plot
plt.show()
