# -*- coding: utf-8 -*-
"""
The :mod:`anisotropy.effective_modelling` module provides a suite of functions
to model the effective anisotropic properties of composite materials. This
includes:
    * Hudson's model for a material containing penny-shaped cracks
    * Tandon and Weng's model for a material containing ellipsoidal inclusions
    * Silver and Savage's model for N-layer modelling
    * Backus' model for a material comprised of alternating layers of elastic
      material

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from .n_layers import (ElasticLayer, _calculate_misfit, _silver_savage_94,
                       _plot_effective_splitting)
