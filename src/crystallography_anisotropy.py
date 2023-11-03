import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_velocities(tensor, n):
    T = np.zeros((3, 3)) # Initialize an empty 3*3 matrix

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    T[i][k] += tensor[i][j][k][l] * n[j] * n[l]

    eigenvalues = np.linalg.eigvals(T)

    eigenvalues = np.sort(eigenvalues)

    rho = 3500 
    M = eigenvalues / rho
    velocities = np.sqrt(M)

    vp = velocities[2] 
    vs1 = velocities[1] 
    vs2 = velocities[0] 

    return (vp, vs1, vs2)

tensor = [[[[100, -10, -10], [-10, -10, -10], [-10, -10, -10]], [[-10, -10, -10], [-10, -10, -10], [-10, -10, -10]], [[-10, -10, -10], [-10, -10, -10], [-10, -10, -10]]], [[[-10, -10, -10], [-10, -20, -20], [-10, -20, -20]], [[-20, -20, -20], [200, -20, -20], [-20,-20,-20]], [[-20,-20,-20], [-20,-20,-20], [200,-20,-20]]], [[[-10,-10,-10], [-20,-20,-20], [-20,-20,-30]], [[-20,-20,-20], [-20,-30,-30], [-30,-30,-30]], [[-30,-30,-30], [-30,-30,-40], [300,-30,-30]]]]

theta = np.linspace(0, np.pi, 50) 
phi = np.linspace(0, 2*np.pi, 50) 
theta_grid, phi_grid = np.meshgrid(theta, phi) 

#basis vectors in crystallographic cordinate
a1 = np.array([1/2., 0., 0.]) 
a2 = np.array([0., 1/2., 0.]) 
a3 = np.array([0., 0., 1/2.]) 

n_grid = np.zeros((3, len(theta), len(phi))) 

for i in range(len(theta)):
    for j in range(len(phi)):
        r = 1 #radius of sphere
        t = theta_grid[i][j] 
        p = phi_grid[i][j]
        
        x = r * np.sin(t) * np.cos(p)
        y = r * np.sin(t) * np.sin(p)
        z = r * np.cos(t)
        
        # Convert the cartesian coordinates to crystallographic coordinates 
        n = x * a1 + y * a2 + z * a3
        
        # Normalize the crystallographic coordinates to get a unit vector
        n = n / np.linalg.norm(n)
        
        n_grid[:, i, j] = n

vp_grid = np.zeros_like(theta_grid) 
vs1_grid = np.zeros_like(theta_grid)
vs2_grid = np.zeros_like(theta_grid)

for i in range(len(theta)):
    for j in range(len(phi)):
        n = n_grid[:, i, j]
        vp, vs1, vs2 = calculate_velocities(tensor, n)
        vp_grid[i][j] = vp
        vs1_grid[i][j] = vs1
        vs2_grid[i][j] = vs2


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Reds(vp_grid), alpha=0.5)
ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Greens(vs1_grid), alpha=0.5)
ax.plot_surface(n_grid[0], n_grid[1], n_grid[2], facecolors=plt.cm.Blues(vs2_grid), alpha=0.5)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_title('Velocities on a 3d sphere')

plt.show()
print((vp_grid))
