import unittest
import numpy as np
from santex.anisotropy import Anisotropy
from santex.material import Material
from santex.tensor import Tensor
from santex.isotropy import Isotropy


class TestMaterialAnisotropy(unittest.TestCase):
    def setUp(self):
        # Anisotropy and Materiual test

        """
        Test seismic property calculations for Forsterite at specified P-T conditions. 
        based on 
        1. Elastic moduli, pressure derivatives, and temperature derivatives of single-crystal olivine and single-crystal forsterite
        10.1029/JB074i025p05961
        2. Faccenda, Manuele, et al. "ECOMAN: an open-source package for geodynamic and seismological modelling of mechanical anisotropy." Solid Earth 15.10 (2024): 1241-1264.
        """
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
        """
        Test seismic property calculations for Forsterite at specified P-T conditions.

        Expected values are based on:
        Hacker, B.R., & Abers, G.A. (2004). 
        Subduction Factory 3: An Excel worksheet and macro for calculating the 
        densities, seismic wave speeds, and H2O contents of minerals and rocks 
        at pressure and temperature. 
        Geochemistry, Geophysics, Geosystems, 5(1). https://doi.org/10.1029/2003GC000614
        """

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

    def test_seismic_anisotropy(self):
        self.cij_forsterite = np.array([
            [320.5, 68.15, 71.6, 0, 0, 0],
            [68.15, 196.5, 76.8, 0, 0, 0],
            [71.6, 76.8, 233.5, 0, 0, 0],
            [0, 0, 0, 64, 0, 0],
            [0, 0, 0, 0, 77, 0],
            [0, 0, 0, 0, 0, 78.7]
        ])
        self.average_density = 3355
        self.anisotropy = Anisotropy(self.cij_forsterite*10**9, self.average_density)
        self.anisotropy_values = self.anisotropy.anisotropy_values()

        self.assertAlmostEqual(self.anisotropy_values['maxvp'], 9773.8, places=0)
        self.assertAlmostEqual(self.anisotropy_values['minvp'], 7653.1, places=0)
        self.assertAlmostEqual(self.anisotropy_values['maxvs1'], 5459.3, places=0)
        self.assertAlmostEqual(self.anisotropy_values['minvs1'], 4790.8, places=0)
        self.assertAlmostEqual(self.anisotropy_values['maxvs2'], 4832.7, places=0)
        self.assertAlmostEqual(self.anisotropy_values['minvs2'], 4367.6, places=0)
        
        



if __name__ == '__main__':
    unittest.main()
