import unittest
import numpy as np
from santex import Material, Anisotropy, Tensor, Isotropy

class TestMaterialAnisotropy(unittest.TestCase):
    def setUp(self):
        # Anisotropy and Materiual test
        self.material = Material()
        self.rock = ["Forsterite", "Diopside", "Enstatite"]
        self.fraction = [0.6, 0.25, 0.15]
        self.expected_tensor = np.array([
            [281.828, 75.0585, 74.522, 0., 2.095, 0.],
            [75.0585, 185.4755, 74.368, 0., 2.445, 0.],
            [74.522, 74.368, 226.12, 0., 13.1425, 0.],
            [0., 0., 0., 65.984, 0., 2.75],
            [2.095, 2.445, 13.1425, 0., 71.4465, 0.],
            [0., 0., 0., 2.75, 0., 74.602]
        ])
        self.expected_density = 3177.1809637264096
        self.expected_maxvp = 9420.564966931908
        
        # Tensor test
        self.tensor = Tensor()
        self.cij_forsterite = np.array([
            [320.5, 68.15, 71.6, 0, 0, 0],
            [68.15, 196.5, 76.8, 0, 0, 0],
            [71.6, 76.8, 233.5, 0, 0, 0],
            [0, 0, 0, 64, 0, 0],
            [0, 0, 0, 0, 77, 0],
            [0, 0, 0, 0, 0, 78.7]
        ])
        self.expected_cijkl_forsterite = np.array([
            [[[320.5, 0., 0.], [0., 68.15, 0.], [0., 0., 71.6]],
             [[0., 78.7, 0.], [78.7, 0., 0.], [0., 0., 0.]],
             [[0., 0., 77.], [0., 0., 0.], [77., 0., 0.]]],
            [[[0., 78.7, 0.], [78.7, 0., 0.], [0., 0., 0.]],
             [[68.15, 0., 0.], [0., 196.5, 0.], [0., 0., 76.8]],
             [[0., 0., 0.], [0., 0., 64.], [0., 64., 0.]]],
            [[[0., 0., 77.], [0., 0., 0.], [77., 0., 0.]],
             [[0., 0., 0.], [0., 0., 64.], [0., 64., 0.]],
             [[71.6, 0., 0.], [0., 76.8, 0.], [0., 0., 233.5]]]
        ])

        # Isotropy test
        self.isotropy = Isotropy()

    def test_isotropy_calculate_seismic_properties(self):
        density, aks, amu, vp, vs, vbulk, akt = self.isotropy.calculate_seismic_properties(
            'Forsterite', temperature=2000, pressure=2, return_vp_vs_vbulk=True, return_aktout=True)
        expected_density = 3034.1119267366676
        expected_vp = 7.5781436531609305
        expected_vs = 4.294565839495191
        expected_vbulk = 5.730428867439992

        self.assertAlmostEqual(density, expected_density, places=5)
        self.assertAlmostEqual(vp, expected_vp, places=5)
        self.assertAlmostEqual(vs, expected_vs, places=5)
        self.assertAlmostEqual(vbulk, expected_vbulk, places=5)

if __name__ == '__main__':
    unittest.main()
