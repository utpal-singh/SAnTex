# `santex.tensor` — Tensor Operations

Convert between Voigt notation and full 4th-order tensor, rotate tensors by Euler angles.

```python
from santex.tensor import Tensor
t = Tensor()
```

## Class `Tensor`

### `voigt_to_tensor(cij: np.ndarray) → np.ndarray`

Convert a 6 × 6 Voigt matrix to a 3 × 3 × 3 × 3 stiffness tensor.

```python
import numpy as np
cij = np.array([
    [320.5,  68.15,  71.6, 0, 0, 0],
    [ 68.15, 196.5,  76.8, 0, 0, 0],
    [ 71.6,   76.8, 233.5, 0, 0, 0],
    [  0,      0,     0,  64, 0, 0],
    [  0,      0,     0,   0,77, 0],
    [  0,      0,     0,   0, 0,78.7],
])
cijkl = t.voigt_to_tensor(cij)   # shape (3,3,3,3)
```

### `tensor_to_voigt(cijkl: np.ndarray) → np.ndarray`

Convert a 3 × 3 × 3 × 3 tensor back to 6 × 6 Voigt form.

```python
cij_back = t.tensor_to_voigt(cijkl)   # shape (6,6)
```

### `rotate_tensor(cijkl: np.ndarray, alpha: float, beta: float, gamma: float) → np.ndarray`

Rotate a 4th-order tensor by Bunge ZXZ Euler angles (degrees).

The rotation follows:

$$C'_{ijkl} = R_{ia}\,R_{jb}\,R_{kc}\,R_{ld}\,C_{abcd}$$

where **R** = R_Z(α) · R_X(β) · R_Z(γ).

```python
alpha, beta, gamma = 10, 45, 30   # degrees
cijkl_rot = t.rotate_tensor(cijkl, alpha, beta, gamma)
cij_rot = t.tensor_to_voigt(cijkl_rot)
```

## Voigt notation index map

| Tensor indices | Voigt index |
|----------------|-------------|
| 11 | 1 |
| 22 | 2 |
| 33 | 3 |
| 23, 32 | 4 |
| 13, 31 | 5 |
| 12, 21 | 6 |

The full conversion is:

$$C_{IJ} = C_{ijkl} \qquad \text{using the index map above}$$

with the factor-of-2 convention handled automatically.

## Rotation matrices

**Z-axis rotation** by θ:

$$R_Z(\theta) = \begin{pmatrix}\cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1\end{pmatrix}$$

**X-axis rotation** by φ:

$$R_X(\phi) = \begin{pmatrix}1 & 0 & 0 \\ 0 & \cos\phi & -\sin\phi \\ 0 & \sin\phi & \cos\phi\end{pmatrix}$$

Full Bunge ZXZ: **R** = R_Z(α) · R_X(β) · R_Z(γ)
