# Import numpy library for array operations and math library for trigonometric functions
import numpy as np
import math

# Import matplotlib library for plotting and mplot3d toolkit for 3d plots
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

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


def plot_vp_2d(cijkl, rho):
    # Create a figure and a 2d axes object
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Define a step size for sampling the angles in radians
    step = math.pi / 180
    # Initialize empty lists for x, y coordinates and colours
    x = []
    y = []
    c = []
    # Loop through all possible values of theta from 0 to pi/2 with the step size
    for theta in np.arange(0, math.pi / 2 + step, step):
        # Loop through all possible values of phi from 0 to 2*pi with the step size
        for phi in np.arange(0, 2 * math.pi + step, step):
            # Calculate the propagation direction vector n from theta and phi using spherical coordinates
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            # Calculate the christoffel tensor for this direction
            tik = christoffel_tensor(cijkl, n)
            # Calculate the wave moduli and polarization directions for this direction
            wave_moduli, polarization_directions = wave_property(tik)
            # Calculate vp from the first wave modulus and density using square root formula
            vp = math.sqrt(wave_moduli[0] / rho)
            # Calculate the stereographic projection of the propagation direction vector using tan formula
            x.append(n[0] / (1 + n[2]))
            y.append(n[1] / (1 + n[2]))
            # Append vp as the colour value to the list
            c.append(vp)
    # Plot the scatter plot of x, y coordinates with colours c using red-blue colour map
    sc = ax.scatter(x, y, c=c, cmap='RdBu')
    # Set the axis labels and title
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('vp in 2d stereographic projection plots with velocity represented by colour')
    # Add a colour bar with label and format
    cb = plt.colorbar(sc)
    cb.set_label('vp')
    cb.ax.set_yticklabels(['{:.1f}'.format(v) for v in cb.get_ticks()])
    # Show the plot
    plt.show()

if __name__=="__main__":
    import tensor_conversion
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                    [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                    [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                    [  0.,      0.,      0.,     65.66,    0.,      6.415],
                    [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                    [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9

    cijkl = tensor_conversion.voigt_to_tensor(M)

    rho = 3500


    # Define a function to plot the vp in 2d stereographic projection plots with velocity represented by colour


    # Plot the vp in 2d stereographic projection plots with velocity represented by colour
    plot_vp_2d(cijkl, rho)
