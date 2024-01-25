import numpy as np
import matplotlib.pyplot as plt

# Define the domain
x_values = np.linspace(0, 17, 21)
y_values = np.linspace(1423, 2500, 21)

# Generate the 2D matrix
matrix = np.zeros((len(y_values), len(x_values)))

# Fill the matrix with values
for i, y in enumerate(y_values):
    for j, x in enumerate(x_values):
        # Given equations
        y1 = -4309 * x**2 / 714 + 115253 * x / 714 + 1423
        y2 = -2959 * x**2 / 714 + 73403 * x / 714 + 1873

        # Interpolate values between y1 and y2
        if y < y1:
            matrix[i, j] = 0
        elif y > y2:
            matrix[i, j] = 1
        else:
            matrix[i, j] = (y - y1) / (y2 - y1)

# Display the 2D matrix
print("2D Matrix:")
print(matrix)

# Plot the 2D matrix
plt.imshow(matrix, cmap='viridis', extent=[0, 17, 1423, 2500], origin='lower', aspect='auto')
plt.colorbar(label='Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.title('2D Matrix Visualization')
plt.show()
