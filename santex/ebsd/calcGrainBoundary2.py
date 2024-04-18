import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from joblib import Parallel, delayed

def euler_to_quaternion(phi, theta, psi):
    """
    Convert Euler angles (in radians) to a quaternion.
    """
    c1 = np.cos(phi / 2)
    s1 = np.sin(phi / 2)
    c2 = np.cos(theta / 2)
    s2 = np.sin(theta / 2)
    c3 = np.cos(psi / 2)
    s3 = np.sin(psi / 2)

    q0 = c1 * c2 * c3 + s1 * s2 * s3
    q1 = s1 * c2 * c3 - c1 * s2 * s3
    q2 = c1 * s2 * c3 + s1 * c2 * s3
    q3 = c1 * c2 * s3 - s1 * s2 * c3

    return np.array([q0, q1, q2, q3])

def quaternion_to_angle(q1, q2):
    """
    Calculate the misorientation angle between two quaternions.
    """
    dot_product = np.dot(q1, q2)
    angle_radians = 2 * np.arccos(np.abs(dot_product))
    return np.degrees(angle_radians)

def misorientation_angle(euler1, euler2):
    """
    Calculate the misorientation angle between two sets of Euler angles.
    """
    q1 = euler_to_quaternion(*np.radians(euler1))
    q2 = euler_to_quaternion(*np.radians(euler2))
    return quaternion_to_angle(q1, q2)

def assign_to_grains_parallel(df, threshold):
    grains = []
    grain_indices = {}
    total_rows = len(df)
    processed_rows = 0

    def process_row(row):
        nonlocal grains, grain_indices, processed_rows

        euler_angles = row[['Euler1', 'Euler2', 'Euler3']].values
        assigned = False

        for grain_idx, grain_euler in enumerate(grains):
            misorientation = misorientation_angle(euler_angles, grain_euler)
            if misorientation <= threshold:
                grain_indices[row.name] = grain_idx
                grains[grain_idx] = np.mean([grain_euler, euler_angles], axis=0)
                assigned = True
                break

        if not assigned:
            grain_indices[row.name] = len(grains)
            grains.append(euler_angles)

        processed_rows += 1
        progress = (processed_rows / total_rows) * 100
        print(f"Progress: {progress:.2f}%\r", end='')

    Parallel(n_jobs=-1, prefer="threads")(delayed(process_row)(row) for _, row in df.iterrows())

    return grain_indices

# Sample phase names
phase_names = ['nothn', 'Forsterite', 'Enstatite Opx AV77', 'Diopside CaMgSi2O6',
               'Anorthite', 'Hornblende', 'Chromite', 'Pyrope']

# Assign rows to grains in parallel
threshold = 20  # Misorientation threshold
grain_indices = assign_to_grains_parallel(df[['Euler1', 'Euler2', 'Euler3']], threshold)

# Add grain indices to the dataframe
df['Grain'] = pd.Series(grain_indices)

# Plot grains and grain boundaries
plt.figure(figsize=(10, 8))

unique_grains = df['Grain'].unique()
cmap = cm.get_cmap('viridis', len(unique_grains))
colors = [cmap(i) for i in range(len(unique_grains))]

for grain, color in zip(unique_grains, colors):
    grain_pixels = df[df['Grain'] == grain]
    plt.scatter(grain_pixels['X'], grain_pixels['Y'], color=color, s=10, label=f'Grain {grain}')

grain_boundaries = df[df['Grain'].diff() != 0]
plt.scatter(grain_boundaries['X'], grain_boundaries['Y'], color='black', s=10, label='Grain Boundary')

plt.legend(fontsize=8)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Grains and Grain Boundaries')
plt.show()

# Display grains with their respective phase names
for grain_idx, group in df.groupby('Grain'):
    phase_counts = group['Phase'].value_counts()
    dominant_phase = phase_counts.idxmax()
    print(f"Grain {grain_idx}: Dominant Phase - {phase_names[dominant_phase]}, Size - {len(group)}")