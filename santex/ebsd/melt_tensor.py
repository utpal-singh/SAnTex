import numpy as np
from ..tensor.tensor import Tensor

def calc_melt_tensor():
    """
    Returns the melt stiffness tensor.
    
    Parameters: 
    - None
    
    Returns:
    - np.ndarray - melt stifness tensor
    
    Reference 
    Lee, A. L., Walker, A. M., Lloyd, G. E., & Torvela, T. (2017). 
    Modeling the impact of melt on seismic properties during mountain building. 
    Geochemistry, Geophysics, Geosystems, 18(3), 1090–1110. https://doi.org/10.1002/2016GC006705
    """
    voigt_melt = np.array([[16.1,   15.967,  15.967,   0.,      0,   0.   ],
                                [ 15.967, 16.1,   15.967,    0.,      0,   0.   ],
                                [ 15.967,  15.967,  16.1,    0.,     0,    0.   ],
                                [  0.,      0.,      0.,     0.01,    0.,      0],
                                [  0,   0,  0,    0.,     0.01,    0.   ],
                                [  0.,      0.,     0.,      0,   0.,     0.01 ]])



    tensor = Tensor()

    tensor_melt = tensor.voigt_to_tensor(voigt_melt)
    return tensor_melt, voigt_melt
