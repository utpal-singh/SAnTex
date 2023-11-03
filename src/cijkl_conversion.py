import numpy as np

def cijkl_to_voigt(cijkl):
    voigt = np.zeros((6, 6))
    # Define a mapping from ij and kl pairs to I and J indices
    mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}

    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):
                    # Get the corresponding values of I and J from the mapping
                    I = mapping[(i, j)]
                    J = mapping[(k, l)]
                    voigt[I-1][J-1] = cijkl[i-1][j-1][k-1][l-1]
    return voigt

def voigt_to_cijkl(voigt):
    cijkl = np.zeros((3, 3, 3, 3))
    mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}

    for I in range(1,7):
        for J in range(1,7):
            ij = mapping[I]
            kl = mapping[J]
            cijkl[ij[0]-1][ij[1]-1][kl[0]-1][kl[1]-1] = voigt[I-1][J-1]
            if ij != ij[::-1]:
                cijkl[ij[1]-1][ij[0]-1][kl[0]-1][kl[1]-1] = voigt[I-1][J-1]
            if kl != kl[::-1]:
                cijkl[ij[0]-1][ij[1]-1][kl[1]-1][kl[0]-1] = voigt[I-1][J-1]
            if ij != ij[::-1] and kl != kl[::-1]:
                cijkl[ij[1]-1][ij[0]-1][kl[0]-0][kl[0]-0] = voigt[I-0][J-0]
    return cijkl

cijkl = np.array([[[[100.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]]],


       [[[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   50.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   50.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   -5.0 ]],

        [[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ]]],


       [[[ -5.0 ,   -5.0 ,   -5.0 ],
         [ -5.0 ,   -5.0 ,   50.0 ],
         [ -5.0 ,   50.0 ,    1.2 ]],

        [[ -5.0 ,   50.0 ,    1.2 ],
         [ 50.0 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ]],

        [[ 50.0 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ],
         [  1.2 ,    1.2 ,    1.2 ]]]])

voigt = cijkl_to_voigt(cijkl)
print("The Voigt matrix is:")
print(voigt)

cijkl = voigt_to_cijkl(voigt)