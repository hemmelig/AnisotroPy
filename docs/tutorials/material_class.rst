Introduction to the `anisotropy.materials.Material` class
=========================================================
In this tutorial, we will cover the basics of the `~anisotropy.materials.Material` class. This class performs the task of representing an arbitrarily (an)isotropic elastic material and provides methods for computing a number of elastic moduli e.g. the bulk and shear moduli, group and phase velocities, and performing arbitrary rotations.

As detailed in the general background on materials and elasticity, a linearly elastic material can be defined completely by an elastic stiffness tensor, :math:`C_{ijkl}` (a 4th-rank tensor), and a density, :math:`\rho`. In practice, these are often measured in the lab and reported in publications, though it is possible to derive theoretical single-crystal values. **FIND A GOOD LEAPING OFF POINT FOR THOSE INTERESTED IN THIS**

A material can also be provided with an assortment of metadata, including: a material symmetry (to be implemented), a material ID, and a reference + DOI (Digital Object Identifier; to be implemented) if the values are derived from literature.

Creating a new `Material`
-------------------------
To create a new `Material` object, you need to provide a stiffness tensor in the reduced Voigt form and a density. Here is an example for creating a `Material` object representing single-crystal olivine:

::

    from anisotropy.materials import Material

    # Moduli and density are from Abramson et al., 1997
    C = [[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
         [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
         [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]]
    rho = 3.355

    material = Material(C, rho)
    print(material)

We now have an object whose various elastic properties we can query. The `print` call returns a neat summary of the material.

Loading an existing `Material`
------------------------------
There is a large body of work dedicated to the experimental determination of the elastic moduli of minerals (single-crystals) and rocks (aggregates). `AnisotroPy` provides a database of published material descriptions that can be simply loaded:

::

    from anisotropy.materials import load

    material = load("olivine")
    print(material)

As above, the `print` call returns a neat summary of the material. We can compare two materials by simply using the `==` operator

::

    from anisotropy.materials import load, Material

    C = [[320.5,  68.1,  71.6,  0.0,  0.0,  0.0],
         [ 68.1, 196.5,  76.8,  0.0,  0.0,  0.0],
         [ 71.6,  76.8, 233.5,  0.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0, 64.0,  0.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0, 77.0,  0.0],
         [  0.0,   0.0,   0.0,  0.0,  0.0, 78.7]]
    rho = 3.355

    material1 = Material(C, rho)
    material2 = load("olivine")

    print(material1 == material2)

Exploring the elastic properties of a `Material`
------------------------------------------------
Phase velocities
################
For example, we can calculate the seismic phase speeds for a plane wave propagating along the the crystallographic "a" ([100]), "b" ([010]), and "c" ([001]) axes:

::

    from anisotropy.materials import load

    material = load("olivine")

    # Calculate phase speeds in olivine for a plane wave propagating along:
    # the "a" ([100]) axis
    a_vp, a_vs1, a_vs2, *_ = material.phase_speeds(0, 0)

    # the "b" ([010]) axis
    b_vp, b_vs1, b_vs2, *_ = material.phase_speeds(0, 90)

    # the "c" ([001]) axis
    c_vp, c_vs1, c_vs2, *_ = material.phase_speeds(90, 0)

    print("The phase speeds (in km/s) for a plane wave propagating along the a/b/c axes:")
    print(f" Vp: a - {a_vp:5.3f}; b - {b_vp:5.3f}; c - {c_vp:5.3f}")
    print(f"Vs1: a - {a_vs1:5.3f}; b - {b_vs1:5.3f}; c - {c_vs1:5.3f}")
    print(f"Vs2: a - {a_vs2:5.3f}; b - {b_vs2:5.3f}; c - {c_vs2:5.3f}")

Bulk and shear moduli
#####################
The bulk and shear moduli of an isotropic material are measures of how resistant that material is to compression and shearing, respectively.

The bulk modulus, often denoted as K or B, is defined as the ratio of the infinitesimal pressure increase to the resulting relative decrease of the volume. This property appears in the isotropic equation for the P-wave velocity.

The shear modulus, often denoted as G, S, or μ, is defined as the ratio of shear stress to the shear strain. This property appears in the isotropic equations for both the P- and S-wave velocities.

Summary
-------
In the next tutorial, we will cover how to build materials that represent aggregates of multiple constitutive materials, for example for the mantle.
