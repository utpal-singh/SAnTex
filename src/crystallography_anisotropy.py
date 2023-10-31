# Import numpy for matrix operations and linear algebra
import numpy as np
# Import matplotlib for plotting
import matplotlib.pyplot as plt
# Import mpl_toolkits for 3d projection
from mpl_toolkits.mplot3d import Axes3D

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

# Define a sample stiffness tensor
tensor = [[[[100, -10, -10], [-10, -10, -10], [-10, -10, -10]], [[-10, -10, -10], [-10, -10, -10], [-10, -10, -10]], [[-10, -10, -10], [-10, -10, -10], [-10, -10, -10]]], [[[-10, -10, -10], [-10, -20, -20], [-10, -20, -20]], [[-20, -20, -20], [200, -20, -20], [-20,-20,-20]], [[-20,-20,-20], [-20,-20,-20], [200,-20,-20]]], [[[-10,-10,-10], [-20,-20,-20], [-20,-20,-30]], [[-20,-20,-20], [-20,-30,-30], [-30,-30,-30]], [[-30,-30,-30], [-30,-30,-40], [300,-30,-30]]]]

# Define a grid of propagation directions on a 3d sphere using spherical coordinates
theta = np.linspace(0, np.pi, 50) # Angle from z-axis
phi = np.linspace(0, 2*np.pi, 50) # Angle from x-axis
theta_grid, phi_grid = np.meshgrid(theta, phi) # Create a grid of angles

# Convert the spherical coordinates to crystallographic coordinates using [1]
a1 = np.array([1/2., 0., 0.]) # The first basis vector of the crystallographic coordinate system
a2 = np.array([0., 1/2., 0.]) # The second basis vector of the crystallographic coordinate system
a3 = np.array([0., 0., 1/2.]) # The third basis vector of the crystallographic coordinate system

n_grid = np.zeros((3, len(theta), len(phi))) # Initialize an empty grid for n

# Loop through the grid indices
for i in range(len(theta)):
    for j in range(len(phi)):
        # Get the spherical coordinates at each grid point
        r = 1 # The radius of the sphere is 1
        t = theta_grid[i][j] # The angle from z-axis
        p = phi_grid[i][j] # The angle from x-axis
        
        # Convert the spherical coordinates to cartesian coordinates using [2]
        x = r * np.sin(t) * np.cos(p)
        y = r * np.sin(t) * np.sin(p)
        z = r * np.cos(t)
        
        # Convert the cartesian coordinates to crystallographic coordinates using [3]
        n = x * a1 + y * a2 + z * a3
        
        # Normalize the crystallographic coordinates to get a unit vector
        n = n / np.linalg.norm(n)
        
        # Assign the unit vector to the grid
        n_grid[:, i, j] = n

# Calculate the velocities for each direction on the grid using the function
vp_grid = np.zeros_like(theta_grid) # Initialize an empty grid for vp
vs1_grid = np.zeros_like(theta_grid) # Initialize an empty grid for vs1
vs2_grid = np.zeros_like(theta_grid) # Initialize an empty grid for vs2

# Loop through the grid indices
for i in range(len(theta)):
    for j in range(len(phi)):
        # Get the direction vector at each grid point
        n = n_grid[:, i, j]
        # Calculate the velocities using the function
        vp, vs1, vs2 = calculate_velocities(tensor, n)
        # Assign the velocities to the grids
        vp_grid[i][j] = vp
        vs1_grid[i][j] = vs1
        vs2_grid[i][j] = vs2

# Plot the velocities on a spherical plot using matplotlib

# Create a figure and a 3d axis object
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the vp grid as a surface plot with a red color map
ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Reds(vp_grid), alpha=0.5)
# Plot the vs1 grid as a surface plot with a green color map
ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Greens(vs1_grid), alpha=0.5)
# Plot the vs2 grid as a surface plot with a blue color map
ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Blues(vs2_grid), alpha=0.5)

# Set the axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Set the title
ax.set_title('Velocities on a 3d sphere')

# Show the plot
plt.show()
print((vp_grid))
