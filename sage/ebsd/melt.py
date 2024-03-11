import numpy as np
from ..tensor.tensor import Tensor


def calcMelttensor():
    voigt_melt = np.array([[16.1,   15.967,  15.967,   0.,      0,   0.   ],
                                [ 15.967, 16.1,   15.967,    0.,      0,   0.   ],
                                [ 15.967,  15.967,  16.1,    0.,     0,    0.   ],
                                [  0.,      0.,      0.,     0.01,    0.,      0],
                                [  0,   0,  0,    0.,     0.01,    0.   ],
                                [  0.,      0.,     0.,      0,   0.,     0.01 ]])



    tensor = Tensor()

    tensor_melt = tensor.voigt_to_tensor(voigt_melt)
    return tensor_melt