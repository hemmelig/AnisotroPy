# -*- coding: utf-8 -*-
"""
Modelling the effective elastic properties of a material consisting of an
isotropic host matrix with unidirectionally aligned isotropic spheroid
inclusions.

For more information, please read:

    Tandon, G.P. and Weng, G.J., 1984. The effect of aspect ratio of inclusions
    on the elastic properties of unidirectionally aligned composites. Polymer
    composites, 5(4), pp.327-333.

If you use this code, please cite this article.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from anisotropy.materials import Material


def model(matrix, inclusions, alpha, fraction_inc):
    """
    Calculate the averaged elastic parameters as described in Tandon and Weng
    (1984), hereafter TW84.

    Their theory describes the effective elastic properties of a medium
    consisting of an isotropic host matrix with unidirectionally aligned
    isotropic spheroid inclusions.

    Parameters
    ----------
    matrix : `anisotropy.Material` object
        An isotropic material that forms the host matrix for the inclusions,
        whose elastic properties are fully described by its stiffness tensor
        and density.
    inclusions : `anisotropy.Material` object
        An isotropic material that occupies the inclusions, whose elastic
        properties are fully described by its stiffness tensor and density.
    alpha : float
        Aspect ratio of inclusions.
    fraction_inc : float
        Volumetric fraction of the inclusions within the composite material.

    Returns
    -------
    effective_material : `anisotropy.Material` object
        The composite material with effective elastic properties.

    """

    # Weighted average density
    rho_bulk = (1. - fraction_inc)*matrix.rho + fraction_inc*inclusions.rho

    la, mu = matrix.lame_coefficients
    lai, mui = inclusions.lame_coefficients

    # Young's modulus for the matrix
    E0 = mu*(3.*la + 2.*mu) / (la + mu)

    # Poisson's ratio of the matrix
    nu = la / (2.*(la + mu))

    # D1, D2, and D3 from TW84 (just before equation 18)
    D1 = 1. + 2.*(mui - mu) / (lai - la)
    D2 = (la + 2.*mu) / (lai - la)
    D3 = la / (lai - la)

    # Simplifying terms
    a1 = 1. - nu
    a2 = 1. - 2.*nu
    a3 = alpha**2.
    a4 = a3 - 1.

    # g and g' terms from TW84
    # g  = spheroidal inclusions (alpha > 1)
    # g' = disc-like inclusions (alpha < 1)
    if alpha > 1.:
        g = (alpha / a4**1.5)*(alpha*a4**0.5 - np.arccosh(alpha))
    elif alpha < 1.:
        g = (alpha / (-a4)**1.5)*(np.arccos(alpha) - alpha*(-a4)**0.5)

    # Eshelby's Sijkl tensor
    s11 = (a2 + (3.*a3 - 1.)/a4 - (a2 + (3.*a3)/a4)*g) / (2.*a1)
    s22 = s33 = (3.*a3) / (8.*a1*a4) + ((a2 - 9./(4.*a4))*g) / (4.*a1)
    s23 = s32 = (a3 / (2.*a4) - (a2 + 3./(4.*a4))*g) / (4.*a1)
    s21 = s31 = (-a3) / (2.*a1*a4) + (3.*a3 / a4 - a2)*g / (4.*a1)
    s12 = s13 = (a2 + 1./a4)/(-2.*a1) + (a2 + 3./(2.*a4))*g / (2.*a1)
    s44 = (a3 / (2.*a4) + (a2 - 3./(4.*a4))*g) / (4.*a1)
    s55 = s66 = (a2 - (a3 + 1.)/a4 - 0.5*(a2 - (3.*(a3 + 1.))/a4)*g) / (4.*a1)

    # B term from TW84 (after equation 17)
    B1 = fraction_inc*D1 + D2 + (1. - fraction_inc)*(D1*s11 + 2.*s21)
    B2 = fraction_inc + D3 + (1. - fraction_inc)*(D1*s12 + s22 + s23)
    B3 = fraction_inc + D3 + (1. - fraction_inc)*(s11 + (1. + D1)*s21)
    B4 = fraction_inc*D1 + D2 + (1. - fraction_inc)*(s12 + D1*s22 + s23)
    B5 = fraction_inc + D3 + (1. - fraction_inc)*(s12 + s22 + D1*s23)

    # A terms from TW84 (after equation 20)
    A1 = D1*(B4 + B5) - 2.*B2
    A2 = (1. + D1)*B2 - (B4 + B5)
    A3 = B1 - D1*B3
    A4 = (1. + D1)*B1 - 2.*B3
    A5 = (1. - D1) / (B4 - B5)
    A = 2.*B2*B3 - B1*(B4 + B5)

    # Longitudinal Young's modulus E11
    E11 = E0 / (1. + fraction_inc*(A1 + 2.*nu*A2)/A)

    # Transverse Young's modulus E22
    E22 = E0 / (1. + fraction_inc*(-2.*nu*A3 + (1. - nu)*A4
                                   + (1. + nu)*A5*A) / (2.*A))

    # In-plane shear modulus mu12
    mu12 = 1. + fraction_inc / ((mu / (mui - mu)) + 2.*(1. - fraction_inc)*s66)
    mu12 *= mu

    # Out-plane shear modulus mu23
    mu23 = 1. + fraction_inc / ((mu / (mui - mu)) + 2.*(1. - fraction_inc)*s44)
    mu23 *= mu

    # Equation for nu31 (Sayers, 1992 (36))
    nu31 = nu - (fraction_inc*(nu*(A1 + 2.*nu*A2) + (A3 - nu*A4))
                 / (A + fraction_inc*(A1 + 2.*nu*A2)))

    # Plane-strain bulk modulus K23
    numer = (1. + nu)*a2
    denom = (1. - nu*(1. + 2.*nu31)
             + fraction_inc*(2.*(nu31 - nu)*A3
                  + (1. - nu*(1. + 2.*nu31))*A4) / A)
    K0 = la + mu  # Plane-strain bulk modulus of matrix
    K23 = K0*numer / denom
    nu12_sq = E11/E22 - (1./mu23 + 1./K23)*E11/4.

    # Calculate C per Sayers (1992) equations (24)--(29), rotated into a frame
    # in which the inclusions are orientated vertically
    C = np.array(
        [[E11 + 4.*nu12_sq*K23, 2.*nu31*K23, 2.*nu31*K23,    0,    0,    0],
         [         2.*nu31*K23,  mu23 + K23, -mu23 + K23,    0,    0,    0],
         [         2.*nu31*K23, -mu23 + K23,  mu23 + K23,    0,    0,    0],
         [                   0,           0,           0, mu23,    0,    0],
         [                   0,           0,           0,    0, mu12,    0],
         [                   0,           0,           0,    0,    0, mu12]])

    return Material(C, rho_bulk)
