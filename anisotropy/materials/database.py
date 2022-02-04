# -*- coding: utf-8 -*-
"""
Module containing a collection of published elastic material descriptions, in
terms of elastic stiffness tensor, C, and bulk density, rho.

All values are reported in units of GPa and g/cm^3.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from anisotropy.materials import Material
import anisotropy.utils.errors as errors


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

    try:
        C, rho = globals()[f"_{material_id}"]()
        return Material(np.asarray(C), rho)
    except KeyError:
        raise errors.InvalidMaterialID(material_id)


def _antigorite():
    """
    Elastic constants of the mineral antigorite (GPa) in Voigt notation.
    From:
      Bezacier, L., Reynard, B., Bass, J.D., Sanchez-Valle, C. and Van de
      Moortèle, B., 2010. Elasticity of antigorite, seismic detection of
      serpentinites, and anisotropy in subduction zones. Earth and Planetary
      Science Letters, 289(1-2), pp.198-208.

    """

    C = [[208.4,  66.2, 15.9,   0.0,  2.4,   0.0],
         [ 66.2, 201.6,  5.0,   0.0, -4.4,   0.0],
         [ 15.9,   5.0, 96.7,   0.0,  2.5,   0.0],
         [  0.0,   0.0,  0.0,  17.4,  0.0, -13.1],
         [  2.4,  -4.4,  2.5,   0.0, 18.3,   0.0],
         [  0.0,   0.0,  0.0, -13.1,  0.0,  65.0]]

    rho = 2.620

    return C, rho


def _biotite():
    """
    Elastic constants of the mineral biotite (GPa) in Voigt notation.
    From:
      Aleksandrov and Ryzhova (1986).

    """

    C = [[186.0,  32.4, 11.6, 0.0, 0.0,  0.0],
         [ 32.4, 186.0, 11.6, 0.0, 0.0,  0.0],
         [ 11.6,  11.6, 54.0, 0.0, 0.0,  0.0],
         [  0.0,   0.0,  0.0, 5.8, 0.0,  0.0],
         [  0.0,   0.0,  0.0, 0.0, 5.8,  0.0],
         [  0.0,   0.0,  0.0, 0.0, 0.0, 76.8]]

    rho = 2.800

    return C, rho


def _blueschist_felsic():
    """
    Elastic constants of felsic blueschist (GPa) in Voigt notation.
    From:
      Cao, Y., Jung, H. and Song, S., 2013. Petro‐fabrics and seismic properties
      of blueschist and eclogite in the North Qilian suture zone, NW China:
      Implications for the low‐velocity upper layer in subducting slab,
      trench‐parallel seismic anisotropy, and eclogite detectability in the
      subduction zone. Journal of Geophysical Research: Solid Earth, 118(6),
      pp.3037-3058.

    """

    C = [[149.85,  38.70,  32.59, -0.15, -1.00, -0.19],
         [ 38.70, 163.55,  30.03,  1.05, -1.81, -1.78],
         [ 32.59,  30.03, 121.62,  0.22, -0.95, -0.13],
         [ -0.15,   1.05,   0.22, 48.03, -0.63, -1.14],
         [ -1.00,  -1.81,  -0.95, -0.63, 48.62, -0.01],
         [ -0.19,  -1.78,  -0.13, -1.14, -0.01, 58.42]]

    rho = 2.970

    return C, rho


def _blueschist_mafic():
    """
    Elastic constants of mafic blueschist (GPa) in Voigt notation.
    From:
      Cao, Y., Jung, H. and Song, S., 2013. Petro‐fabrics and seismic properties
      of blueschist and eclogite in the North Qilian suture zone, NW China:
      Implications for the low‐velocity upper layer in subducting slab,
      trench‐parallel seismic anisotropy, and eclogite detectability in the
      subduction zone. Journal of Geophysical Research: Solid Earth, 118(6),
      pp.3037-3058.

    """

    C = [[190.79,  62.28,  52.94, -0.44,  4.68,  0.60],
         [ 62.28, 218.38,  53.10, -0.87,  1.57,  0.28],
         [ 52.94,  53.10, 158.04, -0.44,  2.66, -0.35],
         [ -0.44,  -0.87,  -0.44, 60.86, -0.29,  1.86],
         [  4.68,   1.57,   2.66, -0.29, 58.94, -0.20],
         [  0.60,   0.28,  -0.35,  1.86, -0.20, 69.63]]

    rho = 3.190

    return C, rho


def _clinopyroxene_92():
    """
    Elastic constants of the mineral clinopyroxene (GPa) in Voigt notation.
    From:
      Bhagat, S.S., Bass, J.D. and Smyth, J.R., 1992. Single‐crystal elastic
      properties of omphacite‐C2/c by Brillouin spectroscopy. Journal of
      Geophysical Research: Solid Earth, 97(B5), pp.6843-6848.

    """

    C = [[257.3,  85.9,  76.2,  0.0,  7.1,  0.0],
         [ 85.9, 216.2,  71.8,  0.0, 13.3,  0.0],
         [ 76.2,  71.8, 260.2,  0.0, 33.7,  0.0],
         [  0.0,   0.0,   0.0, 80.2,  0.0, 10.2],
         [  7.1,  13.3,  33.7,  0.0, 70.6,  0.0],
         [  0.0,   0.0,   0.0, 10.2,  0.0, 85.8]]

    rho = 3.327

    return C, rho


def _clinopyroxene_98():
    """
    Elastic constants of the mineral clinopyroxene (GPa) in Voigt notation.
    From:
      Collins, M.D. and Brown, J.M., 1998. Elasticity of an upper mantle
      clinopyroxene. Physics and chemistry of minerals, 26(1), pp.7-13.

    """

    C = [[237.8,  83.5,  80.0,  0.0,  9.0,  0.0],
         [ 83.5, 183.6,  59.9,  0.0,  9.5,  0.0],
         [ 80.0,  59.9, 229.5,  0.0, 48.1,  0.0],
         [  0.0,   0.0,   0.0, 76.5,  0.0,  8.4],
         [  9.0,   9.5,  48.1,  0.0, 73.0,  0.0],
         [  0.0,   0.0,   0.0,  8.4,  0.0, 81.6]]

    rho = 3.190

    return C, rho


def _dolomite():
    """
    Elastic constants of dolomite mineral (GPa) in Voigt notation.
    From:
      Humbert, P. and Plicque, F., 1972. Elastic properties of monocrystalline
      rhombohedral carbonates-Calcite, Magnesite, Dolomite. Comptes rendus
      Hebdomadaires des Seances de L Academie des Sciences Serie B, 275(11),
      p.391.

    """

    C = [[205.0,  71.0,  57.4, -19.5,  13.7,   0.0],
         [ 71.0, 205.0,  57.4,  19.5, -13.7,   0.0],
         [ 57.4,  57.4, 113.0,   0.0,   0.0,   0.0],
         [-19.5,  19.5,   0.0,  39.8,   0.0, -13.7],
         [ 13.7, -13.7,   0.0,   0.0,  39.8, -19.5],
         [  0.0,   0.0,   0.0, -13.7, -19.5,  67.0]]

    rho = 2.840

    return C, rho


def _eclogite_foliated():
    """
    Elastic constants of foliated eclogite rock (GPa) in Voigt notation.
    From:
      Cao, Y., Jung, H. and Song, S., 2013. Petro‐fabrics and seismic properties
      of blueschist and eclogite in the North Qilian suture zone, NW China:
      Implications for the low‐velocity upper layer in subducting slab,
      trench‐parallel seismic anisotropy, and eclogite detectability in the
      subduction zone. Journal of Geophysical Research: Solid Earth, 118(6),
      pp.3037-3058.

    """

    C = [[203.45,  67.76,  64.47,  0.08,  1.90, -0.40],
         [ 67.76, 220.58,  63.65,  0.46,  0.59,  0.06],
         [ 64.47,  63.65, 189.75,  0.13,  0.95, -0.20],
         [  0.08,   0.46,   0.13, 66.32, -0.27,  0.73],
         [  1.90,   0.59,   0.95, -0.27, 65.77, -0.02],
         [ -0.40,   0.06,  -0.20,  0.73, -0.02, 70.75]]

    rho = 3.300

    return C, rho


def _eclogite_massive():
    """
    Elastic constants of massive eclogite rock (GPa) in Voigt notation.
    From:
      Cao, Y., Jung, H. and Song, S., 2013. Petro‐fabrics and seismic properties
      of blueschist and eclogite in the North Qilian suture zone, NW China:
      Implications for the low‐velocity upper layer in subducting slab,
      trench‐parallel seismic anisotropy, and eclogite detectability in the
      subduction zone. Journal of Geophysical Research: Solid Earth, 118(6),
      pp.3037-3058.

    """

    C = [[238.85,  82.01,  81.44,  0.30, -0.02,  0.50],
         [ 82.01, 242.12,  81.11, -0.66,  0.33,  0.12],
         [ 81.44,  81.11, 235.57, -0.28,  0.22,  0.31],
         [  0.30,  -0.66,  -0.28, 78.72,  0.27,  0.00],
         [ -0.02,   0.33,   0.22,  0.27, 78.37,  0.25],
         [  0.50,   0.12,   0.31,  0.00,  0.25, 77.91]]

    rho = 3.490

    return C, rho


def _epidote():
    """
    Elastic constants of the mineral epidote (GPa) in Voigt notation.
    From:
      Aleksandrov et al. (1974).

    """

    C = [[211.5,  65.6,  43.2,  0.0,  -6.5,  0.0],
         [ 65.6, 239.0,  43.6,  0.0, -10.4,  0.0],
         [ 43.2,  43.6, 202.1,  0.0, -20.0,  0.0],
         [  0.0,   0.0,   0.0, 39.1,   0.0, -2.3],
         [ -6.5, -10.4, -20.0,  0.0,  43.4,  0.0],
         [  0.0,   0.0,   0.0, -2.3,   0.0, 79.5]]

    rho = 3.465

    return C, rho


def _garnet():
    """
    Elastic constants of the mineral garnet (GPa) in Voigt notation.
    From:
      Babuška, V., Fiala, J., Kumazawa, M., Ohno, I. and Sumino, Y., 1978.
      Elastic properties of garnet solid-solution series. Physics of the Earth
      and Planetary Interiors, 16(2), pp.157-176.

    """

    C = [[306.2, 112.5, 112.5,  0.0,  0.0,  0.0],
         [112.5, 306.2, 112.5,  0.0,  0.0,  0.0],
         [112.5, 112.5, 306.2,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 92.7,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 92.7,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 92.7]]

    rho = 3.660

    return C, rho


def _glaucophane():
    """
    Elastic constants of the mineral glaucophane  (GPa) in Voigt notation.
    From:
      Bezacier, L., Reynard, B., Bass, J.D., Wang, J. and Mainprice, D., 2010.
      Elasticity of glaucophane, seismic velocities and anisotropy of the
      subducted oceanic crust. Tectonophysics, 494(3-4), pp.201-210.

    """

    C = [[122.3,  45.7,  37.20,  0.0,  2.30,  0.0],
         [ 45.7, 231.5,  74.90,  0.0, -4.80,  0.0],
         [ 37.2,  74.9, 254.60,  0.0, -2.37,  0.0],
         [  0.0,   0.0,   0.00, 79.6,  0.00,  8.9],
         [  2.3,  -4.8,  -2.37,  0.0, 52.80,  0.0],
         [  0.0,   0.0,   0.00,  8.9,  0.00, 51.2]]

    rho = 3.070

    return C, rho


def _harzburgite():
    """
    Elastic constants of harzburgite rock (GPa) in Voigt notation.
    From:
      Covey‐Crump, S.J., Schofield, P.F., Stretton, I.C., Knight, K.S. and
      Ismaïl, W.B., 2003. Using neutron diffraction to investigate the elastic
      properties of anisotropic rocks: results from an olivine+ orthopyroxene
      mylonite. Journal of Geophysical Research: Solid Earth, 108(B2).

    Note:
      This is the Voigt-averaged elastic stiffness tensor for the Oman
      Harzburgite. The paper also calculates the Reuss-averaged elastic
      stiffness tensor.

    """

    C = [[226.50,  75.34,  74.73, -0.27, -2.99,  1.85],
         [ 75.34, 242.80,  73.68, -3.60, -1.91,  4.14],
         [ 74.73,  73.68, 230.00, -4.36, -4.27, -0.27],
         [ -0.27,  -3.60,  -4.36, 80.75,  1.81, -2.19],
         [ -2.99,  -1.91,  -4.27,  1.81, 76.94, -1.88],
         [  1.85,   4.14,  -0.27, -2.19, -1.88, 79.15]]

    rho = 3.200

    return C, rho


def _hornblende():
    """
    Elastic constants of the mineral hornblende (GPa) in Voigt notation.
    From:
      Aleksandrov and Ryzhova (1986).

    """

    C = [[116.0,  49.9,  61.4,  0.0,  4.3,  0.0],
         [ 49.9, 159.7,  65.5,  0.0, -2.5,  0.0],
         [ 61.4,  65.5, 191.6,  0.0, 10.0,  0.0],
         [  0.0,   0.0,   0.0, 57.4,  0.0, -6.2],
         [  4.3,  -2.5,  10.0,  0.0, 31.8,  0.0],
         [  0.0,   0.0,   0.0, -6.2,  0.0, 36.8]]

    rho = 3.200

    return C, rho


def _jadeite():
    """
    Elastic constants of the mineral jadeite (GPa) in Voigt notation.
    From:
      Kandelin, J. and Weidner, D.J., 1988. The single-crystal elastic
      properties of jadeite. Physics of the Earth and Planetary Interiors,
      50(3), pp.251-260.

    """

    C = [[274.0,  94.0,  71.0,  0.0,  4.0,  0.0],
         [ 94.0, 253.0,  82.0,  0.0, 14.0,  0.0],
         [ 71.0,  82.0, 282.0,  0.0, 28.0,  0.0],
         [  0.0,   0.0,   0.0, 88.0,  0.0, 13.0],
         [  4.0,  14.0,  28.0,  0.0, 65.0,  0.0],
         [  0.0,   0.0,   0.0, 13.0,  0.0, 94.0]]

    rho = 3.330

    return C, rho


def _lawsonite():
    """
    Elastic constants of the mineral lawsonite (GPa) in Voigt notation.
    From:
      Sinogeikin, S.V., Schilling, F.R. and Bass, J.D., 2000. Single crystal
      elasticity of lawsonite. American Mineralogist, 85(11-12), pp.1834-1837.

    """

    C = [[226.0,  69.0,  65.0,  0.0,  0.0,  0.0],
         [ 69.0, 214.0,  82.0,  0.0,  0.0,  0.0],
         [ 65.0,  82.0, 259.0,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 65.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 60.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 17.0]]

    rho = 3.090

    return C, rho


def _lherzolite():
    """
    Elastic constants of lherzolite rock (GPa) in Voigt notation.
    From:
      Peselnick, L., Nicolas, A. and Stevenson, P.R., 1974. Velocity anisotropy
      in a mantle peridotite from the Ivrea Zone: Application to upper mantle
      anisotropy. Journal of Geophysical Research, 79(8), pp.1175-1182.

    """

    C = [[187.40,  63.71,  63.87,  0.78,  2.02, -3.20],
         [ 63.71, 211.25,  64.50, -3.07,  0.87, -5.78],
         [ 63.87,  64.50, 190.03,  0.38,  2.38, -0.12],
         [  0.78,  -3.07,   0.38, 67.92, -2.12,  1.60],
         [  2.02,   0.87,   2.38, -2.12, 63.12, -0.55],
         [ -3.20,  -5.78,  -0.12,  1.60, -0.55, 66.83]]

    rho = 3.270

    return C, rho


def _lizardite_atom():
    """
    Elastic constants of the mineral lizardite (GPa) in Voigt notation.
    Derived from atomistic calculations in:
      Auzende, A.L., Pellenq, R.M., Devouard, B., Baronnet, A. and Grauby, O.,
      2006. Atomistic calculations of structural and elastic properties of
      serpentine minerals: the case of lizardite. Physics and chemistry of
      minerals, 33(4), pp.266-275.

    """

    C = [[229.0800,  89.0440, 13.5580, -0.0001,  4.6025,  0.0001],
         [ 89.0440, 229.0800, 13.5570, -0.0001, -4.6016,  0.0001],
         [ 13.5580,  13.5570, 45.8380, -0.0001,  0.0015,  0.0001],
         [ -0.0001,  -0.0001, -0.0001, 12.7650, -0.0001, -4.4598],
         [  4.6025,  -4.6016,  0.0015, -0.0001, 12.7740,  0.0001],
         [  0.0001,   0.0001,  0.0001, -4.4598,  0.0001, 70.0166]]

    rho = 2.5155

    return C, rho


def _lizardite():
    """
    Elastic constants of mineral lizardite (GPa) in Voigt notation.
    Derived from Density Functional Theory in:
      Reynard, B., Hilairet, N., Balan, E. and Lazzeri, M., 2007. Elasticity of 
      serpentines and extensive serpentinization in subduction zones.
      Geophysical Research Letters, 34(13).

    """

    C = [[245.0,  50.0, 31.0,  0.0,  0.0,  0.0],
         [ 50.0, 245.0, 31.0,  0.0,  0.0,  0.0],
         [ 31.0,  31.0, 23.0,  0.0,  0.0,  0.0],
         [  0.0,   0.0,  0.0, 11.6,  0.0,  0.0],
         [  0.0,   0.0,  0.0,  0.0, 11.6,  0.0],
         [  0.0,   0.0,  0.0,  0.0,  0.0, 97.5]]

    rho = 2.610

    return C, rho


def _muscovite():
    """
    Elastic constants of the mineral muscovite (GPa) in Voigt notation.
    From:
      Vaughan, M.T. and Guggenheim, S., 1986. Elasticity of muscovite and its
      relationship to crystal structure. Journal of Geophysical Research: Solid
      Earth, 91(B5), pp.4657-4664.

    """

    C = [[181.0,  48.8,  25.6,   0.0, -14.2,  0.0],
         [ 48.8, 178.4,  21.2,   0.0,   1.1,  0.0],
         [ 25.6,  21.2,  58.6,   0.0,   1.0,  0.0],
         [  0.0,   0.0, -14.2,  16.5,   0.0, -5.2],
         [-14.2,   1.1,   1.0,   0.0,  19.5,  0.0],
         [  0.0,   0.0,   0.0,  -5.2,   0.0, 72.0]]

    rho = 2.834

    return C, rho


def _olivine():
    """
    Elastic constants of the mineral olivine (GPa) in Voigt notation.
    From:
      Abramson, E.H., Brown, J.M., Slutsky, L.J. and Zaug, J., 1997. The elastic
      constants of San Carlos olivine to 17 GPa. Journal of Geophysical
      Research: Solid Earth, 102(B6), pp.12253-12263.

    """

    C = [[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
         [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
         [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]]

    rho = 3.355

    return C, rho


def _orthopyroxene():
    """
    Elastic constants of the mineral orthopyroxene (GPa) in Voigt notation.
    From:
      Chai, M., Brown, J.M. and Slutsky, L.J., 1997. The elastic constants of an
      aluminous orthopyroxene to 12.5 GPa. Journal of Geophysical Research:
      Solid Earth, 102(B7), pp.14779-14785.

    """

    C = [[236.9,  79.6,  63.2,  0.0,  0.0,  0.0],
         [ 79.6, 180.5,  56.8,  0.0,  0.0,  0.0],
         [ 63.2,  56.8, 230.4,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 84.3,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 79.4,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 80.1]]

    rho = 3.304

    return C, rho


def _plagioclase_64():
    """
    Elastic constants of the mineral plagioclase (GPa) in Voigt notation.
    From:
      Ryzhova (1964).

    """

    C = [[ 81.8,  39.3,  40.7,  0.0, -9.0,  0.0],
         [ 39.3, 145.0,  34.1,  0.0, -7.9,  0.0],
         [ 40.7,  34.1, 133.0,  0.0,-18.5,  0.0],
         [  0.0,   0.0,   0.0, 17.7,  0.0, -0.8],
         [ -9.0,  -7.9, -18.5,  0.0, 31.2,  0.0],
         [  0.0,   0.0,   0.0, -0.8,  0.0, 33.3]]

    rho = 2.700

    return C, rho


def _plagioclase_06():
    """
    Elastic constants of the mineral plagioclase (GPa) in Voigt notation.
    From:
      Brown et al. (2006).

    """

    C = [[69.90,  33.240,  31.56,   5.28, -2.46,  -0.720],
         [33.24, 183.280,   7.53,   5.31, -7.60,  -0.423],
         [31.56,   7.530, 175.65, -17.48,  5.86, -11.290],
         [ 5.28,   5.310, -17.48,  26.93, -3.94,  -6.560],
         [-2.46,  -7.600,   5.86,  -3.94, 26.91,   0.980],
         [-0.72,  -0.423, -11.29,  -6.56,  0.98,  33.390]]

    rho = 2.700

    return C, rho


def _quartz():
    """
    Elastic constants of the mineral quartz (GPa) in Voigt notation.
    From:
      Lakshanov et al. (2007)

    """

    C = [[86.9,   7.6,  12.0,  17.8,   0.0,   0.0],
         [ 7.6,  86.9,  12.0, -17.8,   0.0,   0.0],
         [12.0,  12.0, 106.4,   0.0,   0.0,   0.0],
         [17.8, -17.8,   0.0,  59.5,   0.0,   0.0],
         [ 0.0,   0.0,   0.0,   0.0,  59.5, -17.8],
         [ 0.0,   0.0,   0.0,   0.0, -17.8,  39.6]]

    rho = 2.649

    return C, rho


def _serpentinite_37():
    """
    Elastic constants of serpentinite rock sample HPS-M (GPa) in Voigt notation.
    From:
      Watanabe, T., Shirasugi, Y., Yano, H. and Michibayashi, K., 2011. Seismic
      velocity in antigorite-bearing serpentinite mylonites. Geological Society,
      London, Special Publications, 360(1), pp.97-112.
    Mineralogy:
        ``Ol`` (57.7%), ``Atg`` (36.9%), ``Trm`` (4.5%), ``Mgt`` (1.1%)
    Note:
      This is the Voigt-averaged elastic stiffness tensor for the HPS-M
      sample. The paper also calculates the Reuss-averaged and Voigt-Reuss-Hill-
      averaged elastic stiffness tensors.

    """

    C = [[205.52,  66.36,  62.29, -0.10, -1.48,  3.86],
         [ 66.36, 195.79,  65.23, -0.37,  0.20,  1.54],
         [ 62.29,  65.23, 193.30, -1.78, -0.24,  0.83],
         [ -0.10,  -0.37,  -1.78, 66.17,  1.47, -0.57],
         [ -1.48,   0.20,  -0.24,  1.47, 64.70, -0.84],
         [  3.86,   1.54,   0.83, -0.57, -0.84, 67.83]]

    rho = 3.000

    return C, rho


def _serpentinite_80():
    """
    Elastic constants of serpentinite rock sample HKB-B (GPa) in Voigt notation.
    From:
      Watanabe, T., Shirasugi, Y., Yano, H. and Michibayashi, K., 2011. Seismic
      velocity in antigorite-bearing serpentinite mylonites. Geological Society,
      London, Special Publications, 360(1), pp.97-112.
    Mineralogy:
        ``Ol`` (12.0%), ``Atg`` (80.2%), ``Mgt`` (7.8%)
    Note:
      This is the Voigt-averaged elastic stiffness tensor for the HPS-M
      sample. The paper also calculates the Reuss-averaged and Voigt-Reuss-Hill-
      averaged elastic stiffness tensors.

    """

    C = [[192.25,  49.35,  41.70, -4.55,  8.04,  9.78],
         [ 49.35, 156.90,  42.36, -6.91,  0.71,  1.84],
         [ 41.70,  42.36, 141.62, -4.28,  1.11,  0.19],
         [ -4.55,  -6.91,  -4.28, 53.48,  0.01, -0.06],
         [  8.04,   0.71,   1.11,  0.01, 51.91, -3.72],
         [  9.78,   1.84,   0.19, -0.06, -3.72, 59.13]]

    rho = 2.800

    return C, rho


def _zoisite():
    """
    Elastic constants of the mineral zoisite (GPa) in Voigt notation.
    From:
      Mao, Z., Jiang, F. and Duffy, T.S., 2007. Single-crystal elasticity of
      zoisite Ca2Al3Si3O12 (OH) by Brillouin scattering. American Mineralogist,
      92(4), pp.570-576.

    """

    C = [[279.8,  94.7,  88.7,  0.0,  0.0,  0.0],
         [ 94.7, 249.2,  27.5,  0.0,  0.0,  0.0],
         [ 88.7,  27.5, 209.4,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 51.8,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 81.4,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 66.3]]

    rho = 3.343

    return C, rho
