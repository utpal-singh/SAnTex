import numpy as np
import matplotlib.pyplot as plt

# Define the domain
x_values = np.linspace(0, 13, 21)
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
            matrix[i, j] = 0.05
        elif y > y2:
            matrix[i, j] = 0  # Cap the value at 0.05
        else:
            interpolated_value = ((y - y1) / (y2 - y1))
            matrix[i, j] = min(interpolated_value, 0.05)  # Cap the interpolated value at 0.05

# Multiply the entire matrix by 2


# Display the 2D matrix
print("2D Matrix:")
transmat = (matrix).T
print(transmat)

plt.imshow(transmat, cmap='viridis', extent=[1423, 2500, 0, 17], origin='lower', aspect='auto')
plt.colorbar(label='Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.title('2D Matrix Visualization')
plt.show()

# Define pressure and temperature
pressure = np.linspace(0, 13, 21)
temperature = np.linspace(1423, 2500, 21)
