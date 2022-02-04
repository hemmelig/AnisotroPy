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
        # Try importing an example material
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


if __name__ == "__main__":
    unittest.main()
