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
        # Test importing valid/invalid materials
        material = materials.load("olivine")
        self.assertIsInstance(material, materials.Material)

        C = np.array(
            [[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
             [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
             [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
             [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
             [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
             [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]])

        self.assertTrue(np.all(material.C == C))
        self.assertEqual(material.rho, 3.355)
        self.assertRaises(errors.InvalidMaterialID,
                          materials.load, "not_a_material")

    def test_voigt_contraction(self):
        # Test the Voigt contraction matrix is correct
        # - Load material
        material = materials.load("olivine")

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
                          materials.isotropic_C, K=3.0, la=03.0, rho=3.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, K=3.0, vs=2.0, rho=2.0)
        self.assertRaises(errors.InsufficientElasticInformation,
                          materials.isotropic_C, mu=3.0, vp=2.0)
        print("\t\t ...tests pass!")

        print("All tests for the function 'isotropic_C' have passed!")
        print("="*60)

    def test_birch_law(self):
        # Test the conversion function of Birch law
        from anisotropy.materials.core import _vp2rho
        vp = 3.0
        expected_rho = 1.73

        self.assertEqual(_vp2rho(vp), expected_rho)

if __name__ == "__main__":
    unittest.main()
