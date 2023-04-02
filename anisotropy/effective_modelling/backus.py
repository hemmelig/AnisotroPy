# -*- coding: utf-8 -*-
"""
Modelling the effective elastic properties of a material consisting of repeated
alternating layers of isotropic materials.

For more information, please read:

    Backus, G.E., 1962. Long-wave elastic anisotropy produced by horizontal
    layering. Journal of Geophysical Research, 67(11), pp.4427-4440.

If you use this code, please cite this article.

:copyright:
    2023, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from anisotropy.materials import Material


def model(material1, material2, relative_fraction):
    """
    Calculate the averaged elastic parameters as described in Backus (1962).
    The variable naming system used is consistent with the elastic coefficients in the
    original paper.

    L = C44 = 1/⟨1/μ⟩
    M = C66 = ⟨μ⟩
    R = 1 / C33 = ⟨θ/μ⟩
    S = ⟨θ * μ⟩
    T = ⟨θ⟩

    These map into the entries of the effective elastic matrix, C'ij, as:

    C'11 = 4M - 4S + (1-2T)^2/R
    C'12 = 2M - 4S + (1-2T)^2/R
    C'13 = (1-2T)/R
    C'33 = 1/R
    C'44 = L
    C'66 = M = (C'11 - C'12)/2

    Parameters
    ----------
    material1 : `anisotropy.Material` object
        An isotropic material that comprises layer 1, whose elastic properties
        are fully described by its stiffness tensor and density.
    material2 : `anisotropy.Material` object
        An isotropic material that comprises layer 2, whose elastic properties
        are fully described by its stiffness tensor and density.
    relative_fraction :
        Volumetric fraction of layer 1.

    Returns
    -------
    effective_material : `anisotropy.Material` object
        The composite material with effective elastic properties.

    """

    f1 = relative_fraction
    f2 = 1 - relative_fraction

    # Calculate the dimensionless parameter, theta, for each material.
    # Theta is the square of the Vs/Vp ratio.
    t1 = material1.C[3, 3] / material1.C[0, 0]
    t2 = material2.C[3, 3] / material2.C[0, 0]
    m1 = material1.C[3, 3]
    m2 = material2.C[3, 3]

    # Effective elastic parameters as defined in Backus (1962)
    L = 1 / (f1 / m1 + f2 / m2)
    M = f1 * m1 + f2 * m2
    R = f1 * t1 / m1 + f2 * t2 / m2
    S = f1 * t1 * m1 + f2 * t2 * m2
    T = f1 * t1 + f2 * t2

    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = 4 * M - 4 * S + (1 - 2 * T) ** 2 / R
    C[0, 1] = C[1, 0] = 2 * M - 4 * S + (1 - 2 * T) ** 2 / R
    C[0, 2] = C[1, 2] = C[2, 0] = C[2, 1] = (1 - 2 * T) / R
    C[2, 2] = 1 / R
    C[3, 3] = C[4, 4] = L
    C[5, 5] = M

    rho_bulk = f1 * material1.rho + f2 * material2.rho

    return Material(C, rho_bulk)
