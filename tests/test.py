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

        Expected values are calculated from:
        Hacker, B.R., & Abers, G.A. (2004). 
        Subduction Factory 3: An Excel worksheet and macro for calculating the 
        densities, seismic wave speeds, and H2O contents of minerals and rocks 
        at pressure and temperature. 
        Geochemistry, Geophysics, Geosystems, 5(1). https://doi.org/10.1029/2003GC000614
        """

        density, aks, amu, vp, vs, vbulk, akt = self.isotropy.calculate_seismic_properties(
            'Forsterite', temperature=2000, pressure=2, return_vp_vs_vbulk=True, return_aktout=True)
        expected_density = 3034.11
        expected_vp = 7.58
        expected_vs = 4.29
        expected_vbulk = 5.73

        self.assertAlmostEqual(density, expected_density, places=0)
        self.assertAlmostEqual(vp, expected_vp, places=0)
        self.assertAlmostEqual(vs, expected_vs, places=0)
        self.assertAlmostEqual(vbulk, expected_vbulk, places=0)

    def test_seismic_anisotropy(self):
        """
        Test seismic anisotropy for a single crystal of Forsterite with MTEX

        The expected results in MTEX can be generated from:

        cs_tensor = crystalSymmetry('mmm',[4.7646,10.2296,5.9942],...
        'x||a','z||c','mineral','Olivine');

        M = [[320.5  68.15  71.6     0     0     0];...
            [ 68.15  196.5  76.8     0     0     0];...
            [  71.6   76.8 233.5     0     0     0];...
            [   0      0      0     64     0     0];...
            [   0      0      0      0    77     0];...
            [   0      0      0      0     0  78.7]];

            % Define density (g/cm3)
            rho=3.355;

            % Define tensor object in MTEX
            % Cij -> Cijkl - elastic stiffness tensor
            C = stiffnessTensor(M,cs_tensor,'density',rho);

        [vp,vs1,vs2,pp,ps1,ps2] = C.velocity('harmonic');

        max(vp)
        max(vs1)
        min(vp)
        min(vs1)

        """

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

        self.assertAlmostEqual(self.anisotropy_values['maxvp']/1000, 9.7738, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvp']/1000, 7.6531, places=2)
        self.assertAlmostEqual(self.anisotropy_values['maxvs1']/1000, 5.4593, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvs1']/1000, 4.7908, places=2)
        self.assertAlmostEqual(self.anisotropy_values['maxvs2']/1000, 4.8327, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvs2']/1000, 4.3676, places=2)


    def test_ebsd(self):
        """
        Test seismic anisotropy and ebsd from a free data Forsterite.ctf available from MTEX,
        The expected results in MTEX can be generated from https://mtex-toolbox.github.io/CPOSeismicProperties.html
        """

        import numpy as np
        from santex.ebsd import EBSD
        from santex.tensor import Tensor
        from santex.anisotropy import Anisotropy
        from santex.material import Material

        ebsd = EBSD("../notebooks/Forsterite.ctf")
        df = ebsd.get_ebsd_data()
        df = ebsd.filter_by_phase_number_list(df = df, phase_list = [4, 5, 6, 7])
        material_instance = Material()
        rho_Fo = material_instance.load_density("Forsterite")
        rho_diop = material_instance.load_density("Diopside")
        rho_ens = material_instance.load_density("Enstatite")
        cij_Fo = material_instance.voigt_high_PT('Forsterite')
        cij_ens = material_instance.voigt_high_PT('Enstatite')
        cij_diop = material_instance.voigt_high_PT('Diopside')
        cij = [cij_Fo, cij_ens, cij_diop]
        density = [rho_Fo, rho_ens, rho_diop]
        forsterite = ebsd.get_euler_angles(phase = 1, data=df)
        enstatite = ebsd.get_euler_angles(phase = 2, data=df)
        diopside = ebsd.get_euler_angles(phase = 3, data=df)
        euler_angles = [forsterite, enstatite, diopside]
        average_tensor, average_density = ebsd.get_anisotropy_for_ebsd(cij, euler_angles, density)
        anis = Anisotropy(average_tensor*10**9, average_density)
        self.anisotropy_values = anis.anisotropy_values()

        self.assertAlmostEqual(self.anisotropy_values['maxvp']/1000, 8.8533, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvp']/1000, 7.8634, places=2)
        self.assertAlmostEqual(self.anisotropy_values['maxvs1']/1000,5.0728, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvs1']/1000,4.7689, places=2)
        self.assertAlmostEqual(self.anisotropy_values['maxvs2']/1000,4.8876, places=2)
        self.assertAlmostEqual(self.anisotropy_values['minvs2']/1000,4.6546, places=2)


        
        



if __name__ == '__main__':
    unittest.main()
