# -*- coding: utf-8 -*-
"""
The :mod:`anisotropy.materials` module provides a suite of standard elastic
materials and means of interacting with them.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from .core import (Material, voigt_reuss_hill_average, voigt_average,
                   reuss_average, isotropic_C)
from .database import load
