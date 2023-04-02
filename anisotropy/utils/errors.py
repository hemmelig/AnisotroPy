# -*- coding: utf-8 -*-
"""
Module that supplies custom exceptions for the AnisotroPy package.

:copyright:
    2023, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""


class MissingDensityValue(Exception):
    """
    Custom exception to handle case when the user requests an isotropic material
    without providing a density value.

    """

    def __init__(self):
        super().__init__(
            "Must provide a density when constructing an isotropic material."
        )


class InsufficientElasticInformation(Exception):
    """
    Custom exception to handle case when the user attempts to build an isotropic
    material without provide sufficient elastic moduli. These can be in the form
    of the Lame coefficients, lambda and mu, the isotropic body wave velocities,
    vp and vs, or the elastic bulk and shear moduli, K and G.

    """

    def __init__(self):
        super().__init__(
            "Insufficient elastic information - consult the documentation!"
        )


class InvalidMaterialID(Exception):
    """
    Custom exception to handle case when the user requests the elastic
    properties for a material that is not included in the material database.

    """

    def __init__(self, material_id):
        super().__init__(f"Material with ID {material_id} not in database.")
