# -*- coding: utf-8 -*-
"""
Module that supplies custom exceptions for the AnisotroPy package.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

class InvalidMaterialID(Exception):
    """
    Custom exception to handle case when the user requests the elastic
    properties for a material that is not included in the material database.

    """

    def __init__(self, mID):
        msg = (f"Material with ID {mID!r} not in database.")
        super().__init__(msg)
