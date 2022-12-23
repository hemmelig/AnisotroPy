# -*- coding: utf-8 -*-
"""
Testing the 2-layer anisotropic modelling module, which is validated against
results computed with MSAT.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import unittest

import numpy as np

import anisotropy.materials as materials
import anisotropy.utils.errors as errors


class TestMaterials(unittest.TestCase):
    def test_load_materials(self):
        print("")
        print("="*60)
        print("Testing the function 'load'...")
        # Test importing valid/invalid materials
        material = materials.load("olivine")
        print("\tTest 1 - test material has correct type...")
        self.assertIsInstance(material, materials.Material)
        print("\t\t ...tests pass!")

        expected_C = np.array(
            [[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
             [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
             [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
             [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
             [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
             [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]])
        print("\tTest 2 - test elastic moduli are as expected...")
        self.assertTrue(np.allclose(material.C, expected_C))
        print("\t\t ...tests pass!")

        print("\tTest 3 - test density is as expected...")
        self.assertEqual(material.rho, 3.355)
        print("\t\t ...tests pass!")

        print("\tTest 4 - test invalid material ID raises an error...")
        self.assertRaises(errors.InvalidMaterialID,
                          materials.load, "not_a_material")
        print("\t\t ...tests pass!")

        print("All tests for the function 'load' have passed!")
        print("="*60)

    def test_voigt_contraction(self):
        # Test the Voigt contraction matrix is correct
        # - Load material
        material = materials.load("olivine")

    def test_voigt_reuss_hill_averaging(self):
        print()
        print("="*60)
        print("Testing the suite of Voigt-Reuss-Hill averaging functions...")
        material1 = materials.isotropic_C(vp=5.0, vs=3.0, rho=4.0)
        material2 = materials.isotropic_C(la=5.0, mu=4.0, rho=4.0)
        vrh_C, _, voigt_C, reuss_C = materials.voigt_reuss_hill_average(
            [material1, material2],
            [0.2, 0.8])
        print("\tTest 1 - testing the Voigt averaging function using\n"
              "\t\t two simple isotropic materials...")
        expected_C = np.array(
            [[30.4,  9.6,  9.6,  0.0,  0.0,  0.0],
             [ 9.6, 30.4,  9.6,  0.0,  0.0,  0.0],
             [ 9.6,  9.6, 30.4,  0.0,  0.0,  0.0],
             [ 0.0,  0.0,  0.0, 10.4,  0.0,  0.0],
             [ 0.0,  0.0,  0.0,  0.0, 10.4,  0.0],
             [ 0.0,  0.0,  0.0,  0.0,  0.0, 10.4]]
        )
        self.assertTrue(np.allclose(voigt_C, expected_C))
        print("\t\t ...tests pass!")

        print("\tTest 2 - testing the Reuss averaging function using\n"
              "\t\t two simple isotropic materials...")
        a, b, c = 15.729147, 5.999418, 4.864865
        expected_C = np.array(
            [[a, b, b, 0, 0, 0],
             [b, a, b, 0, 0, 0],
             [b, b, a, 0, 0, 0],
             [0, 0, 0, c, 0, 0],
             [0, 0, 0, 0, c, 0],
             [0, 0, 0, 0, 0, c]]
        )
        self.assertTrue(np.allclose(reuss_C, expected_C))
        print("\t\t ...tests pass!")

        print("\tTest 3 - testing the VRH averaging function using\n"
              "\t\t two simple isotropic materials...")
        a, b, c = 46.129147/2, 15.599418/2, 15.264865/2
        expected_C = np.array(
            [[a, b, b, 0, 0, 0],
             [b, a, b, 0, 0, 0],
             [b, b, a, 0, 0, 0],
             [0, 0, 0, c, 0, 0],
             [0, 0, 0, 0, c, 0],
             [0, 0, 0, 0, 0, c]]
        )
        self.assertTrue(np.allclose(vrh_C, expected_C))
        print("\t\t ...tests pass!")

        print("\tTest 4 - testing a suite of cases where one\n\t\t expects an"
              " exception to be raised...")
        self.assertRaises(ValueError, materials.voigt_reuss_hill_average,
                          [material1, material2], [0.3, 0.8])
        print("\t\t ...tests pass!")

        print("All tests for the Voigt-Reuss-Hill averaging have passed!")
        print("="*60)

    def test_isotropic_elastic_tensor(self):
        print()
        print("="*60)
        print("Testing the function 'isotropic_C'...")
        template = ("\tTest {number} - building an isotropic elastic tensor "
                    "from\n\t\t {arguments}...")

        # Test the method for creating an isotropic elastic tensor with various
        # input types
        # 1. With Vp and Vs
        print(template.format(number=1, arguments="Vp and Vs"))
        isotropic_C = materials.isotropic_C(vp=5.0, vs=3.0)
        expected_C = np.array(
            [[59.25, 16.59, 16.59,  0.00,  0.00,  0.00],
             [16.59, 59.25, 16.59,  0.00,  0.00,  0.00],
             [16.59, 16.59, 59.25,  0.00,  0.00,  0.00],
             [ 0.00,  0.00,  0.00, 21.33,  0.00,  0.00],
             [ 0.00,  0.00,  0.00,  0.00, 21.33,  0.00],
             [ 0.00,  0.00,  0.00,  0.00,  0.00, 21.33]]
        )
        expected_rho = 2.37
        self.assertTrue(isotropic_C.is_isotropic)
        self.assertTrue(np.allclose(isotropic_C.C, expected_C))
        self.assertEqual(isotropic_C.rho, expected_rho)
        print("\t\t ...tests pass!")

        # 2. With Vp, Vs, and rho
        print(template.format(number=2, arguments="Vp, Vs, and rho"))
        isotropic_C = materials.isotropic_C(vp=5.0, vs=3.0, rho=4.0)
        expected_C = np.array(
            [[100,  28,  28,  0,  0,  0],
             [ 28, 100,  28,  0,  0,  0],
             [ 28,  28, 100,  0,  0,  0],
             [  0,   0,   0, 36,  0,  0],
             [  0,   0,   0,  0, 36,  0],
             [  0,   0,   0,  0,  0, 36]]
        )
        self.assertTrue(isotropic_C.is_isotropic)
        self.assertTrue(np.allclose(isotropic_C.C, expected_C))
        print("\t\t ...tests pass!")

        # 3. With the Lamé parameters, lambda and mu
        print(template.format(
            number=3,
            arguments="the Lame parameters, la and mu, and rho"))
        isotropic_C = materials.isotropic_C(la=5.0, mu=4.0, rho=4.0)
        expected_C = np.array(
            [[13.0,  5.0,  5.0, 0.0, 0.0, 0.0],
             [ 5.0, 13.0,  5.0, 0.0, 0.0, 0.0],
             [ 5.0,  5.0, 13.0, 0.0, 0.0, 0.0],
             [ 0.0,  0.0,  0.0, 4.0, 0.0, 0.0],
             [ 0.0,  0.0,  0.0, 0.0, 4.0, 0.0],
             [ 0.0,  0.0,  0.0, 0.0, 0.0, 4.0]]
        )
        self.assertTrue(isotropic_C.is_isotropic)
        self.assertEqual(isotropic_C.lame_coefficients[0], 5.0)
        self.assertEqual(isotropic_C.lame_coefficients[1], 4.0)
        self.assertTrue(np.allclose(isotropic_C.C, expected_C))
        print("\t\t ...tests pass!")

        # 4. With the bulk (K) and shear (G) moduli
        print(template.format(
            number=4,
            arguments="the elastic moduli, K and G, and rho"))
        isotropic_C = materials.isotropic_C(K=8.0, G=6.0, rho=4.0)
        expected_C = np.array(
            [[16.0,  4.0,  4.0, 0.0, 0.0, 0.0],
             [ 4.0, 16.0,  4.0, 0.0, 0.0, 0.0],
             [ 4.0,  4.0, 16.0, 0.0, 0.0, 0.0],
             [ 0.0,  0.0,  0.0, 6.0, 0.0, 0.0],
             [ 0.0,  0.0,  0.0, 0.0, 6.0, 0.0],
             [ 0.0,  0.0,  0.0, 0.0, 0.0, 6.0]]
        )
        self.assertTrue(isotropic_C.is_isotropic)
        self.assertEqual(isotropic_C.bulk_modulus, 8.0)
        self.assertEqual(isotropic_C.shear_modulus, 6.0)
        self.assertTrue(np.allclose(isotropic_C.C, expected_C))
        print("\t\t ...tests pass!")

        # 5. Test for correctly raised Exceptions
        print("\tTest 5 - testing a suite of cases where one\n\t\t expects an"
              " exception to be raised...")
        self.assertRaises(errors.MissingDensityValue,
                          materials.isotropic_C, K=8.0, G=6.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, K=3.0, rho=3.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, vp=3.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, la=3.0, rho=3.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, K=3.0, la=3.0, rho=3.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, K=3.0, vs=2.0, rho=2.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, mu=3.0, vp=2.0)
        print("\t\t ...tests pass!")

        print("All tests for the function 'isotropic_C' have passed!")
        print("="*60)

    def test_hexagonal1_C(self):

        vp0 = 6
        vs0 = 3
        ani = 50
        rho = 5
        incs = np.array([0, 30, 60, 90])
        azis = [0] * len(incs)

        hex1_material = materials.hexagonal1_C(vp0, vs0, ani, rho)
        _, _, vs2, _, _ = hex1_material.phase_velocities(incs, azis)

        # Theta is meassured dwon from the z axis.
        inc = (90 - incs) * np.pi / 180

        # The Hexagonal1 definition of anisotropy
        db = vs0 * ani / 100.0
        LL = rho * (vs0 + db / 2.0) ** 2.0
        NN = rho * (vs0 - db / 2.0) ** 2.0

        # Auld (1973) Acoustic waves in solid. c66 and c33 are here the unique axes
        vsh = 1/(np.sqrt(rho / (LL * np.sin(inc) ** 2 + NN * np.cos(inc) ** 2)))

        self.assertTrue(np.allclose(vs2, vsh))  # this is working

    def test_birch_law(self):
        print()
        print("="*60)
        print("Testing the function '_vp2rho'...")
        # Test the conversion function of Birch law
        print("\tTest 1 - test conversion of Vp to density...")
        from anisotropy.materials.core import _vp2rho
        vp = 3.0
        expected_rho = 1.73

        self.assertEqual(_vp2rho(vp), expected_rho)
        print("\t\t ...tests pass!")
        print("All tests for the function '_vp2rho' have passed!")
        print("="*60)

    def test_phase_velocities(self):
        print()
        print("="*60)
        print("Testing the function 'phase_velocities'...")
        print("\tTest 1 - testing against anisotropic values")
        print("\treported for single crystal antigorite...")
        from anisotropy.materials.database import load
        from anisotropy.materials.core import isotropic_C
        atg = load("antigorite")

        # Table 1 of Bezacier et al.
        # Run no. 10 (note inconsicenty in angles in their table)
        vp, vs1, vs2, _, _ = atg.phase_velocities(90, 0)
        self.assertAlmostEqual(vp[0], 6.08, delta=0.01)
        self.assertAlmostEqual(vs1[0], 2.642, delta=0.01)
        self.assertAlmostEqual(vs2[0], 2.563, delta=0.015)

        # Run no. 906
        vp, vs1, vs2, _, _ = atg.phase_velocities(0, 68.71)
        self.assertAlmostEqual(vp[0], 8.758, delta=0.001)
        self.assertAlmostEqual(vs1[0], 5.107, delta=0.001)
        self.assertAlmostEqual(vs2[0], 2.451, delta=0.002)

        # Now isotropic velocities
        print("\tTest 2 - testing against isotropic values...")
        iatg = atg.isotropic
        
        # Make equal length samples of focal sphere the silly way
        incs = []
        azis = []
        for inc in range(0, 90, 30):
            for azi in range(0, 360, 30):
                incs.append(inc)
                azis.append(azi)

        # Voigt averages of Tab. 3 of Bezacier et al.
        vp, vs1, vs2, _, _ = iatg.phase_velocities(incs, azis)
        for p, s1, s2 in zip(vp, vs1, vs2):
            self.assertAlmostEqual(p, 7.3, delta=0.01)
            self.assertAlmostEqual(s1, 4.29, delta=0.01)
            self.assertAlmostEqual(s2, 4.29, delta=0.01)

        # More isotropic velocities
        print("\tTest 3 - testing against isotropic hill-averages...")
        iatg = isotropic_C(K=atg.bulk_modulus, G=atg.shear_modulus, rho=atg.rho) 

        # Hill averages of Tab. 3 of Bezacier et al.
        vp, vs1, vs2, _, _ = iatg.phase_velocities(incs, azis)
        for p, s1, s2 in zip(vp, vs1, vs2):
            self.assertAlmostEqual(p, 6.76, delta=0.01)
            self.assertAlmostEqual(s1, 3.83, delta=0.01)
            self.assertAlmostEqual(s2, 3.83, delta=0.01)
            
            # Make sure no complex values were returned
            self.assertIsInstance(p, float)
            self.assertIsInstance(s1, float)
            self.assertIsInstance(s2, float)

        print("\t\t ...tests pass!")
        print("All tests for the function 'phase_velocities' have passed!")
        print("="*60)

if __name__ == "__main__":
    unittest.main()
