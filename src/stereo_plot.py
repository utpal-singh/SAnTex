import numpy as np
import matplotlib.pyplot as plt

def calculate_velocities(tensor, n):
    T = np.zeros((3, 3)) 

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

M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9
import tensor_conversion
tensor_M = tensor_conversion.voigt_to_tensor(M)
stiffness = tensor_M

theta = np.linspace(0, 2*np.pi, 50) # Angle from x-axis
n_grid = np.array([np.cos(theta), np.sin(theta), np.zeros_like(theta)]) # Convert to cartesian coordinates
print(n_grid)

vp_grid = np.zeros_like(theta) 
vs1_grid = np.zeros_like(theta) 
vs2_grid = np.zeros_like(theta)

for i in range(len(theta)):
    n = n_grid[:, i]
    vp, vs1, vs2 = calculate_velocities(tensor_M, n)

    vp_grid[i] = vp
    vs1_grid[i] = vs1
    vs2_grid[i] = vs2



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

ax.set_xlabel('x')
ax.set_ylabel('y')

ax.set_title('Velocities on a 2d circle')

plt.show()
