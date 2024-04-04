import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def calculate_odf(euler_angles):
    """
    Calculate and plot an ODF from Bunge Euler angles.
    
    Parameters:
    - euler_angles: A Nx3 numpy array of Euler angles (phi1, Phi, phi2) in degrees.
    
    """
    euler_angles_rad = np.radians(euler_angles)
    
    # using a Gaussian KDE to estimate the distribution
    kde = gaussian_kde(euler_angles_rad.T)
    
    phi1, Phi, phi2 = np.mgrid[0:np.pi:100j, 0:np.pi/2:100j, 0:2*np.pi:100j]
    positions = np.vstack([phi1.ravel(), Phi.ravel(), phi2.ravel()])
    
    density = np.reshape(kde(positions).T, phi1.shape)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(phi1, Phi, phi2, c=density, cmap='viridis')
    ax.set_xlabel('Phi1')
    ax.set_ylabel('Phi')
    ax.set_zlabel('Phi2')
    ax.set_title('ODF Representation (Simplified)')
    plt.show()

if __name__ == "__main__":
    euler_angles = np.array([
        [10, 20, 30],
        [40, 50, 60],
        [20, 25, 46],
        [15, 45, 34],
        [23, 45, 67],
        [45, 56, 34]
    ])
    
    calculate_odf(euler_angles)
