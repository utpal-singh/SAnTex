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

# Given points
x0, y0 = 0, 1873
x1, y1 = 17, 2423
x2, y2 = 14, 2500

# inverted quadratic function
quadratic_function = inverted_quadratic_interpolation(x0, y0, x1, y1, x2, y2)

# quadratic function
print("Inverted Quadratic Function:")
print(quadratic_function)

# x values for plotting
x_values = np.linspace(0, 17, 1000)
y_values = [quadratic_function.subs('x', val) for val in x_values]

# Plot quadratic function
plt.plot(x_values, y_values, label='Inverted Quadratic Function')

# Mark given points
plt.scatter([x0, x1, x2], [y0, y1, y2], color='red', label='Given Points')

# Mark maximum point
plt.scatter(14, 2500, color='green', label='Maximum Point (x=14, y=2500)')

# Add labels and legend
plt.xlabel('x')
plt.ylabel('y')
plt.title('Inverted Quadratic Interpolation')
plt.legend()

# Show plot
plt.show()
