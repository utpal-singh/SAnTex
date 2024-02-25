import numpy as np

def calculate_misorientation(quat1, quat2):
    # Compute misorientation between two quaternions
    dot_product = np.dot(quat1, quat2)
    misorientation = 2 * np.arccos(np.abs(dot_product)) * 180 / np.pi
    return misorientation

def calculate_grain_boundaries(quaternion_array, threshold=5):
    grain_boundaries = np.zeros_like(quaternion_array[:, :, 0])
    
    for i in range(1, quaternion_array.shape[0] - 1):
        for j in range(1, quaternion_array.shape[1] - 1):
            # Get quaternions of neighboring pixels
            center_quat = quaternion_array[i, j]
            top_quat = quaternion_array[i - 1, j]
            bottom_quat = quaternion_array[i + 1, j]
            left_quat = quaternion_array[i, j - 1]
            right_quat = quaternion_array[i, j + 1]
            
            # Calculate misorientation with neighboring pixels
            top_misorientation = calculate_misorientation(center_quat, top_quat)
            bottom_misorientation = calculate_misorientation(center_quat, bottom_quat)
            left_misorientation = calculate_misorientation(center_quat, left_quat)
            right_misorientation = calculate_misorientation(center_quat, right_quat)
            
            # If any misorientation exceeds threshold, mark as grain boundary
            if (top_misorientation > threshold or
                bottom_misorientation > threshold or
                left_misorientation > threshold or
                right_misorientation > threshold):
                grain_boundaries[i, j] = 1
    
    return grain_boundaries


if __name__ == "__main__":
    quaternion_array = ebsd_data.apply(lambda row: euler_to_quaternion(row['Euler1'], row['Euler2'], row['Euler3']), axis=1).values
    # grain_boundaries = calculate_grain_boundaries(ebsd_data['Quaternion'].values)
    grain_boundaries = calculate_grain_boundaries(quaternion_array)

