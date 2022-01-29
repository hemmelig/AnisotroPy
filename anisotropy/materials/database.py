# -*- coding: utf-8 -*-
"""
Module containing a collection of published elastic material descriptions, in
terms of elastic stiffness tensor, C, and bulk density, rho.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from anisotropy.materials import Material


def load(material_id):
    """
    Utility function that serves up a material based on a string identifier.

    Parameters
    ----------
    material_id : str
        Unique material identifier.

    Returns
    -------
    material : `anisotropy.Material` object
        Material loaded from database.

    """

    if material_id == "olivine":
        return _olivine()


def _olivine():
    """
    Returns the Cijkl and density for olivine.

    """

    C = np.array([[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
                  [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
                  [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
                  [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
                  [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
                  [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]])

    rho = 3.355

    return Material(C, rho)
