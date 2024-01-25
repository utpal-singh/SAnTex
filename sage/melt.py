import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, solve

def inverted_quadratic_interpolation(x0, y0, x1, y1, x2, y2):
    x, a, b, c = symbols('x a b c')
    
    # Quadratic equation: y = ax^2 + bx + c
    eq1 = a*x0**2 + b*x0 + c - y0
    eq2 = a*x1**2 + b*x1 + c - y1
    eq3 = a*x2**2 + b*x2 + c - y2
    
    # Solving the system of equations
    solution = solve((eq1, eq2, eq3), (a, b, c))
    
    # Forming the quadratic function
    a_val, b_val, c_val = solution[a], solution[b], solution[c]
    quadratic_function = a_val*x**2 + b_val*x + c_val
    
    return quadratic_function

# Given points for the first set
x0_1, y0_1 = 0, 1423
x1_1, y1_1 = 17, 2423
x2_1, y2_1 = 14, 2500

# Get the inverted quadratic function for the first set
quadratic_function_1 = inverted_quadratic_interpolation(x0_1, y0_1, x1_1, y1_1, x2_1, y2_1)

# Given points for the second set
x0_2, y0_2 = 0, 1873
x1_2, y1_2 = 17, 2423
x2_2, y2_2 = 14, 2500

# Get the inverted quadratic function for the second set
quadratic_function_2 = inverted_quadratic_interpolation(x0_2, y0_2, x1_2, y1_2, x2_2, y2_2)

# Generate x values for plotting
x_values = np.linspace(0, 17, 1000)

# Calculate y values for both functions
y_values_1 = [quadratic_function_1.subs('x', val) for val in x_values]
y_values_2 = [quadratic_function_2.subs('x', val) for val in x_values]

# Plot both quadratic functions
plt.plot(x_values, y_values_1, label='Set 1 - Inverted Quadratic Function')
plt.plot(x_values, y_values_2, label='Set 2 - Inverted Quadratic Function')

# Mark the given points for both sets
plt.scatter([x0_1, x1_1, x2_1], [y0_1, y1_1, y2_1], color='red', label='Set 1 - Given Points')
plt.scatter([x0_2, x1_2, x2_2], [y0_2, y1_2, y2_2], color='blue', label='Set 2 - Given Points')

# Mark the maximum point for both sets
plt.scatter(14, 2500, color='green', label='Maximum Point (x=14, y=2500)')

# Add labels and legend
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solidus and liquidus')
plt.legend()

# Show the plot
plt.show()
