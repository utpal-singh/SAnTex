import numpy as np



def calculate_tensor_deriv(mineral, PRESSURE=0, TEMP=300):
    """
    pressure in gpa
    temperature in kelvin
    """

    if mineral == "diopside":
        tensor = np.zeros((6, 6))
        dCdP = np.zeros((6, 6))
        dCdT = np.zeros((6, 6))
        # dCdT_m2 = np.zeros((6, 6))
        tensor_calc = np.zeros((6, 6))
        tensor[0, 0] = 224.7
        tensor[1, 1] = 177.9
        tensor[2, 2] = 240.5
        tensor[3, 3] = 75.2
        tensor[4, 4] = 67.0
        tensor[5, 5] = 79.1
        tensor[0, 1] = 83.2
        tensor[0, 2] = 72.4
        tensor[1, 2] = 60.2
        tensor[0, 4] = 7.6
        tensor[1, 4] = 4.2
        tensor[2, 4] = 37.7
        tensor[3, 5] = 3.9
        tensor[4, 0] = tensor[0, 4]
        tensor[4, 1] = tensor[1, 4]
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]
        tensor[4, 2] = tensor[2, 4]
        tensor[5, 3] = tensor[3, 5]
        
        dCdP[0, 0] = 6.72
        dCdP[1, 1] = 5.68
        dCdP[2, 2] = 6.78
        dCdP[3, 3] = 1.52
        dCdP[4, 4] = 1.86
        dCdP[5, 5] = 2.56
        dCdP[0, 1] = 4.59
        dCdP[0, 2] = 4.37
        dCdP[1, 2] = 4.06
        dCdP[0, 4] = -0.73
        dCdP[1, 4] = -0.21
        dCdP[2, 4] = 1.50
        dCdP[3, 5] = -0.17
        dCdP[4, 0] = dCdP[0, 4]
        dCdP[4, 1] = dCdP[1, 4]
        dCdP[1, 0] = dCdP[0, 1]
        dCdP[2, 0] = dCdP[0, 2]
        dCdP[2, 1] = dCdP[1, 2]
        dCdP[4, 2] = dCdP[2, 4]
        dCdP[5, 3] = dCdP[3, 5]
        
        dCdT[0, 0] = -26.8
        dCdT[1, 1] = -20.7
        dCdT[2, 2] = -25.4
        dCdT[3, 3] = -8.3
        dCdT[4, 4] = -10.3
        dCdT[5, 5] = -12.1
        dCdT[0, 1] = -11.5
        dCdT[0, 2] = -10.6
        dCdT[1, 2] = -9.9
        dCdT[0, 4] = 1.2
        dCdT[1, 4] = 1.0
        dCdT[2, 4] = 2.1
        dCdT[3, 5] = 4.2
        dCdT[4, 0] = dCdT[0, 4]
        dCdT[4, 1] = dCdT[1, 4]
        dCdT[1, 0] = dCdT[0, 1]
        dCdT[2, 0] = dCdT[0, 2]
        dCdT[2, 1] = dCdT[1, 2]
        dCdT[4, 2] = dCdT[2, 4]
        dCdT[5, 3] = dCdT[3, 5]
        tensor_calc = tensor + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
        return tensor_calc
    
    if mineral == "forsterite":
        tensor = np.zeros((6, 6))
        dCdP = np.zeros((6, 6))
        dCdT = np.zeros((6, 6))
        # dCdT_m2 = np.zeros((6, 6))
        tensor_calc = np.zeros((6, 6))
        tensor[0, 0] = 328.4
        tensor[1, 1] = 199.8
        tensor[2, 2] = 235.3
        tensor[3, 3] = 65.89
        tensor[4, 4] = 81.2
        tensor[5, 5] = 80.88
        tensor[0, 1] = 63.9
        tensor[0, 2] = 68.8
        tensor[1, 2] = 73.8
        # tensor[0, 4] = 7.6
        # tensor[1, 4] = 4.2
        # tensor[2, 4] = 37.7
        # tensor[3, 5] = 3.9
        # tensor[4, 0] = tensor[0, 4]
        # tensor[4, 1] = tensor[1, 4]
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]
        # tensor[4, 2] = tensor[2, 4]
        # tensor[5, 3] = tensor[3, 5]
        
        dCdP[0, 0] = 8.47
        dCdP[1, 1] = 6.56
        dCdP[2, 2] = 6.37
        dCdP[3, 3] = 2.12
        dCdP[4, 4] = 1.66
        dCdP[5, 5] = 2.37
        dCdP[0, 1] = 4.67
        dCdP[0, 2] = 4.84
        dCdP[1, 2] = 4.11
        # dCdP[0, 4] = -0.73
        # dCdP[1, 4] = -0.21
        # dCdP[2, 4] = 1.50
        # dCdP[3, 5] = -0.17
        # dCdP[4, 0] = dCdP[0, 4]
        # dCdP[4, 1] = dCdP[1, 4]
        dCdP[1, 0] = dCdP[0, 1]
        dCdP[2, 0] = dCdP[0, 2]
        dCdP[2, 1] = dCdP[1, 2]
        # dCdP[4, 2] = dCdP[2, 4]
        # dCdP[5, 3] = dCdP[3, 5]
        
        dCdT[0, 0] = -33.1
        dCdT[1, 1] = -28.1
        dCdT[2, 2] = -28.6
        dCdT[3, 3] = -12.8
        dCdT[4, 4] = -13
        dCdT[5, 5] = -15.1
        dCdT[0, 1] = -10.4
        dCdT[0, 2] = -8.2
        dCdT[1, 2] = -4.6
        # dCdT[0, 4] = 1.2
        # dCdT[1, 4] = 1.0
        # dCdT[2, 4] = 2.1
        # dCdT[3, 5] = 4.2
        # dCdT[4, 0] = dCdT[0, 4]
        # dCdT[4, 1] = dCdT[1, 4]
        dCdT[1, 0] = dCdT[0, 1]
        dCdT[2, 0] = dCdT[0, 2]
        dCdT[2, 1] = dCdT[1, 2]
        # dCdT[4, 2] = dCdT[2, 4]
        # dCdT[5, 3] = dCdT[3, 5]
        tensor_calc = tensor + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
        return tensor_calc

    if mineral == "olivine":
        tensor = np.zeros((6, 6))
        dCdP = np.zeros((6, 6))
        dCdT = np.zeros((6, 6))
        # dCdT_m2 = np.zeros((6, 6))
        tensor_calc = np.zeros((6, 6))
        tensor[0, 0] = 323.7
        tensor[1, 1] = 197.6
        tensor[2, 2] = 235.1
        tensor[3, 3] = 64.62
        tensor[4, 4] = 78.05
        tensor[5, 5] = 79.04
        tensor[0, 1] = 66.4
        tensor[0, 2] = 71.6
        tensor[1, 2] = 75.6
        # tensor[0, 4] = 7.6
        # tensor[1, 4] = 4.2
        # tensor[2, 4] = 37.7
        # tensor[3, 5] = 3.9
        # tensor[4, 0] = tensor[0, 4]
        # tensor[4, 1] = tensor[1, 4]
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]
        # tensor[4, 2] = tensor[2, 4]
        # tensor[5, 3] = tensor[3, 5]
        
        dCdP[0, 0] = 7.98
        dCdP[1, 1] = 6.37
        dCdP[2, 2] = 6.38
        dCdP[3, 3] = 2.17
        dCdP[4, 4] = 1.64
        dCdP[5, 5] = 2.31
        dCdP[0, 1] = 4.74
        dCdP[0, 2] = 4.48
        dCdP[1, 2] = 3.76
        # dCdP[0, 4] = -0.73
        # dCdP[1, 4] = -0.21
        # dCdP[2, 4] = 1.50
        # dCdP[3, 5] = -0.17
        # dCdP[4, 0] = dCdP[0, 4]
        # dCdP[4, 1] = dCdP[1, 4]
        dCdP[1, 0] = dCdP[0, 1]
        dCdP[2, 0] = dCdP[0, 2]
        dCdP[2, 1] = dCdP[1, 2]
        # dCdP[4, 2] = dCdP[2, 4]
        # dCdP[5, 3] = dCdP[3, 5]
        
        dCdT[0, 0] = -34
        dCdT[1, 1] = -28.5
        dCdT[2, 2] = -28.6
        dCdT[3, 3] = -12.8
        dCdT[4, 4] = -13
        dCdT[5, 5] = -15.7
        dCdT[0, 1] = -10.5
        dCdT[0, 2] = -9.4
        dCdT[1, 2] = -5.1
        # dCdT[0, 4] = 1.2
        # dCdT[1, 4] = 1.0
        # dCdT[2, 4] = 2.1
        # dCdT[3, 5] = 4.2
        # dCdT[4, 0] = dCdT[0, 4]
        # dCdT[4, 1] = dCdT[1, 4]
        dCdT[1, 0] = dCdT[0, 1]
        dCdT[2, 0] = dCdT[0, 2]
        dCdT[2, 1] = dCdT[1, 2]
        # dCdT[4, 2] = dCdT[2, 4]
        # dCdT[5, 3] = dCdT[3, 5]
        tensor_calc = tensor + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
        return tensor_calc

    if mineral == "enstatite":
        tensor = np.zeros((6, 6))
        dCdP = np.zeros((6, 6))
        dCdT = np.zeros((6, 6))
        # dCdT_m2 = np.zeros((6, 6))
        tensor_calc = np.zeros((6, 6))
        tensor[0, 0] = 228.6
        tensor[1, 1] = 160.5
        tensor[2, 2] = 210.4
        tensor[3, 3] = 81.75
        tensor[4, 4] = 75.48
        tensor[5, 5] = 77.66
        tensor[0, 1] = 78.74
        tensor[0, 2] = 79.6
        tensor[1, 2] = 76.64
        # tensor[0, 4] = 7.6
        # tensor[1, 4] = 4.2
        # tensor[2, 4] = 37.7
        # tensor[3, 5] = 3.9
        # tensor[4, 0] = tensor[0, 4]
        # tensor[4, 1] = tensor[1, 4]
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]
        # tensor[4, 2] = tensor[2, 4]
        # tensor[5, 3] = tensor[3, 5]
        
        dCdP[0, 0] = 11.04
        dCdP[1, 1] = 9.19
        dCdP[2, 2] = 16.42
        dCdP[3, 3] = 2.38
        dCdP[4, 4] = 2.92
        dCdP[5, 5] = 2.75
        dCdP[0, 1] = 6.97
        dCdP[0, 2] = 9.09
        dCdP[1, 2] = 8.73
        # dCdP[0, 4] = -0.73
        # dCdP[1, 4] = -0.21
        # dCdP[2, 4] = 1.50
        # dCdP[3, 5] = -0.17
        # dCdP[4, 0] = dCdP[0, 4]
        # dCdP[4, 1] = dCdP[1, 4]
        dCdP[1, 0] = dCdP[0, 1]
        dCdP[2, 0] = dCdP[0, 2]
        dCdP[2, 1] = dCdP[1, 2]
        # dCdP[4, 2] = dCdP[2, 4]
        # dCdP[5, 3] = dCdP[3, 5]
        
        dCdT[0, 0] = -35.2
        dCdT[1, 1] = -32.8
        dCdT[2, 2] = -51.6
        dCdT[3, 3] = -13.1
        dCdT[4, 4] = -13.8
        dCdT[5, 5] = -14.5
        dCdT[0, 1] = -21.2
        dCdT[0, 2] = -31.8
        dCdT[1, 2] = -10.7
        # dCdT[0, 4] = 1.2
        # dCdT[1, 4] = 1.0
        # dCdT[2, 4] = 2.1
        # dCdT[3, 5] = 4.2
        # dCdT[4, 0] = dCdT[0, 4]
        # dCdT[4, 1] = dCdT[1, 4]
        dCdT[1, 0] = dCdT[0, 1]
        dCdT[2, 0] = dCdT[0, 2]
        dCdT[2, 1] = dCdT[1, 2]
        # dCdT[4, 2] = dCdT[2, 4]
        # dCdT[5, 3] = dCdT[3, 5]
        tensor_calc = tensor + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
        return tensor_calc

        

if __name__ == "__main__":
    print(calculate_tensor_deriv("olivine", 2, 1000))
    # tensor = np.zeros((6, 6))
    # dCdP = np.zeros((6, 6))
    # dCdT = np.zeros((6, 6))
    # dCdT_m2 = np.zeros((6, 6))
    # tensor[0, 0] = 224.7
    # tensor[1, 1] = 177.9
    # tensor[2, 2] = 240.5
    # tensor[3, 3] = 75.2
    # tensor[4, 4] = 67.0
    # tensor[5, 5] = 79.1
    # tensor[0, 1] = 83.2
    # tensor[0, 2] = 72.4
    # tensor[1, 2] = 60.2
    # tensor[0, 4] = 7.6
    # tensor[1, 4] = 4.2
    # tensor[2, 4] = 37.7
    # tensor[3, 5] = 3.9
    # tensor[4, 0] = tensor[0, 4]
    # tensor[4, 1] = tensor[1, 4]
    # tensor[1, 0] = tensor[0, 1]
    # tensor[2, 0] = tensor[0, 2]
    # tensor[2, 1] = tensor[1, 2]
    # tensor[4, 2] = tensor[2, 4]
    # tensor[5, 3] = tensor[3, 5]
    
    # dCdP[0, 0] = 6.72
    # dCdP[1, 1] = 5.68
    # dCdP[2, 2] = 6.78
    # dCdP[3, 3] = 1.52
    # dCdP[4, 4] = 1.86
    # dCdP[5, 5] = 2.56
    # dCdP[0, 1] = 4.59
    # dCdP[0, 2] = 4.37
    # dCdP[1, 2] = 4.06
    # dCdP[0, 4] = -0.73
    # dCdP[1, 4] = -0.21
    # dCdP[2, 4] = 1.50
    # dCdP[3, 5] = -0.17
    # dCdP[4, 0] = dCdP[0, 4]
    # dCdP[4, 1] = dCdP[1, 4]
    # dCdP[1, 0] = dCdP[0, 1]
    # dCdP[2, 0] = dCdP[0, 2]
    # dCdP[2, 1] = dCdP[1, 2]
    # dCdP[4, 2] = dCdP[2, 4]
    # dCdP[5, 3] = dCdP[3, 5]
    
    # dCdT[0, 0] = -26.8
    # dCdT[1, 1] = -20.7
    # dCdT[2, 2] = -25.4
    # dCdT[3, 3] = -8.3
    # dCdT[4, 4] = -10.3
    # dCdT[5, 5] = -12.1
    # dCdT[0, 1] = -11.5
    # dCdT[0, 2] = -10.6
    # dCdT[1, 2] = -9.9
    # dCdT[0, 4] = 1.2
    # dCdT[1, 4] = 1.0
    # dCdT[2, 4] = 2.1
    # dCdT[3, 5] = 4.2
    # dCdT[4, 0] = dCdT[0, 4]
    # dCdT[4, 1] = dCdT[1, 4]
    # dCdT[1, 0] = dCdT[0, 1]
    # dCdT[2, 0] = dCdT[0, 2]
    # dCdT[2, 1] = dCdT[1, 2]
    # dCdT[4, 2] = dCdT[2, 4]
    # dCdT[5, 3] = dCdT[3, 5]
    
    # dCdT_m2[0, 0] = -29.1
    # dCdT_m2[1, 1] = -24.8
    # dCdT_m2[2, 2] = -17.9
    # dCdT_m2[3, 3] = -10.3
    # dCdT_m2[4, 4] = -7.7
    # dCdT_m2[5, 5] = -15.2
    # dCdT_m2[0, 1] = -11.9
    # dCdT_m2[0, 2] = -6.4
    # dCdT_m2[1, 2] = 0
    # dCdT_m2[0, 4] = 2.5
    # dCdT_m2[1, 4] = 2.2
    # dCdT_m2[2, 4] = -4.6
    # dCdT_m2[3, 5] = 2.6
    # dCdT_m2[4, 0] = dCdT_m2[0, 4]
    # dCdT_m2[4, 1] = dCdT_m2[1, 4]
    # dCdT_m2[1, 0] = dCdT_m2[0, 1]
    # dCdT_m2[2, 0] = dCdT_m2[0, 2]
    # dCdT_m2[2, 1] = dCdT_m2[1, 2]
    # dCdT_m2[4, 2] = dCdT_m2[2, 4]
    # dCdT_m2[5, 3] = dCdT_m2[3, 5]
    
    # TEMP = float(input("Enter Temperature: "))
    # PRESSURE = float(input("Enter pressure: "))
    # STUDY = input("Enter dft or exp or both: ")
    
    # tensor_calc = np.zeros((6, 6))
    
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print("Elastic stiffness tensors using ambient conditions at T=300K and P=0GPa")
    # print(tensor)
    
    # if STUDY == "dft":
    #     tensor_calc = tensor + (dCdT * (TEMP - 300) * 0.001) + (dCdP * PRESSURE)
    #     print("------------------------------------------------------------------------------------------------")
    #     print("Elastic stiffness tensors using Quantum Espresso DFT code")
    #     print(tensor_calc)
    # elif STUDY == "exp":
    #     tensor_calc = tensor + (dCdT_m2 * (TEMP - 300) * 10**-3) + (dCdP * PRESSURE)
    #     print("------------------------------------------------------------------------------------------------")
    #     print("Elastic stiffness tensors using previous experimental results")
    #     print(tensor_calc)
    # else:
    #     t_dft = tensor + (dCdT * (TEMP - 300) * 10**-3) + (dCdP * PRESSURE)
    #     t_exp = tensor + (dCdT_m2 * (TEMP - 300) * 10**-3) + (dCdP * PRESSURE)
    #     print("------------------------------------------------------------------------------------------------")
    #     print("Elastic stiffness tensors using Quantum Computation Density Functional Theory code")
    #     print(t_dft)
    #     print("------------------------------------------------------------------------------------------------")
    #     print("Elastic stiffness tensors using previous experimental results")
    #     print(t_exp)
    #     print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")