import sympy as sp

# Define symbols
σ, T, G, e, V, N = sp.symbols('σ T G e V N')

# Define the equations
G = sp.Function('G')(σ, T)
e = -1/V * sp.diff(G, σ)
N = -sp.diff(G, T)
V = sp.diff(G, sp.symbols('P'))

# Display the equations
print("Molar Gibbs energy (G):")
sp.pprint(G)
print("\nFinite strain measure (e):")
sp.pprint(e)
print("\nMolar isostress entropy (N):")
sp.pprint(N)
print("\nMolar volume (V):")
sp.pprint(V)
