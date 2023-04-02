# -*- coding: utf-8 -*-
"""
Modelling the effective elastic properties of a material consisting of an
isotropic host matrix with unidirectionally aligned cracks.

For more information, please read:

    Hudson, J.A., 1980, September. Overall properties of a cracked solid. In
    Mathematical Proceedings of the Cambridge Philosophical Society (Vol. 88, No. 2,
    pp. 371-384). Cambridge University Press.

    Hudson, J.A., 1981. Wave speeds and attenuation of elastic waves in material
    containing cracks. Geophysical Journal International, 64(1), pp.133-150.

    Hudson, J.A., 1986. A higher order approximation to the wave propagation constants
    for a cracked solid. Geophysical Journal International, 87(1), pp.265-274.

If you use this code, please cite these articles.

:copyright:
    2023, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from anisotropy.materials import Material


def model(matrix, inclusions, aspect_ratio, crack_density):
    """
    Calculate the effective elastic properties of a material comprised of an isotropic
    matrix hosting penny-shaped cracks (ellipsoidal inclusions) that may contain a
    different isotropic material. The effective anisotropy is controlled by the crack
    density and the aspect ratio of the inclusions.

    The theory is valid when the crack density (= Na^2 / nu) << 1, where N is the number
    of cracks of radius a in a volume nu.

    The expressions used describe a system that is transversely isotropic about x1, as
    written in Crampin (1984). The original expressions in Hudson (1982) correspond to a
    system that is transversely isotropic about x3. They give the most general
    expressions, but can be tuned to model various different crack types:

    Dry cracks: Set vpi, vsi, and rhoi to 0.
    Water-filled: Set vsi to 0 and vpi*rhoi = 2.25x10E9

    Parameters
    ----------
    matrix : `anisotropy.Material` object
        An isotropic material that forms the host matrix for the inclusions, whose
        elastic properties are fully described by its stiffness tensor and density.
    inclusions : `anisotropy.Material` object
        An isotropic material that occupies the inclusions, whose elastic properties are
        fully described by its stiffness tensor and density.
    aspect_ratio : float
        Aspect ratio of inclusions.
    crack_density : float
        Volumetric fraction of the inclusions within the composite material.

    Returns
    -------
    effective_material : `anisotropy.Material` object
        The composite material with effective elastic properties.

    """

    if crack_density > 0.1:
        print(
            (
                "The expressions in Hudson (1982) are only valid for epsilon"
                " << 0.1. Continuing anyway..."
            )
        )

    la, mu = matrix.lame_coefficients
    _, mui = inclusions.lame_coefficients
    inclusions_K = inclusions.bulk_modulus

    # Calculate theta = la + 2*mu, as per Backus (1962)
    th = la + 2 * mu

    K = ((inclusions_K + 4 * mui / 3) / (np.pi * aspect_ratio * mu)) * th / (la + mu)
    M = ((4 * mui) / (np.pi * aspect_ratio * mu)) * th / (3 * la + 4 * mu)
    U1 = (4 / 3) * (th / (la + mu)) / (1 + K)
    U3 = (16 / 3) * (th / (3 * la + 4 * mu)) / (1 + M)

    # Calculate the first-order correction from Hudson (1981)
    C1 = np.array(
        [
            [U1 * th**2, la * th * U1, la * th * U1, 0, 0, 0],
            [la * th * U1, U1 * la**2, U1 * la**2, 0, 0, 0],
            [la * th * U1, U1 * la**2, U1 * la**2, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, U3 * mu**2, 0],
            [0, 0, 0, 0, 0, U3 * mu**2],
        ]
    )
    C1 *= -crack_density / mu

    # Calculate the second-order correction from Hudson (1982)
    q = 15 * (la / mu) ** 2 + 28 * (la / mu) + 28
    X = 2 * mu * (3 * la + 8 * mu) / th
    C2 = np.array(
        [
            [th * q * U1**2, la * q * U1**2, la * q * U1**2, 0, 0, 0],
            [
                la * q * U1**2,
                q * (la * U1) ** 2 / th,
                q * (la * U1) ** 2 / th,
                0,
                0,
                0,
            ],
            [
                la * q * U1**2,
                q * (la * U1) ** 2 / th,
                q * (la * U1) ** 2 / th,
                0,
                0,
                0,
            ],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, X * U3**2, 0],
            [0, 0, 0, 0, 0, X * U3**2],
        ]
    )
    C2 *= crack_density**2 / 15

    C = matrix.C + C1 + C2
    rho_bulk = (1 - crack_density) * matrix.rho + crack_density * inclusions.rho

    return Material(C, rho_bulk)
