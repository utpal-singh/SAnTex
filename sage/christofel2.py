import numpy as np

# Calculate vp, vs1 and vs2 from stiffness tensor and a propagation direction
def calculate_velocities(tensor, n):

    T = np.zeros((3, 3))
    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    # Sum over repeated indices j and l
                    T[i][k] += tensor[i][j][k][l] * n[j] * n[l]

    eigenvalues = np.linalg.eigvals(T)
    print(T)
    # Sort eigenvalues in ascending
    eigenvalues = np.sort(eigenvalues)


    rho = 3500
    M = eigenvalues / rho #wave moduli - get it by dividing eigenvalues with rho

    # square root of the wave moduli is velocity
    velocities = np.sqrt(M)

    vp = velocities[2] # largest eigenvalue
    vs1 = velocities[1] # middle eigenvalue
    vs2 = velocities[0] # smallest eigenvalue

    return (vp, vs1, vs2)


M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                 [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                 [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                 [  0.,      0.,      0.,     65.66,    0.,      6.415],
                 [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                  [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])
import tensor_conversion
tensor_M = tensor_conversion.voigt_to_tensor(M)
stiffness = tensor_M
n = [0, 1, 0]


vp, vs1, vs2 = calculate_velocities(stiffness, n)

print("vp =", vp)
print("vs1 =", vs1)
print("vs2 =", vs2)
