Mixing `anisotropy.materials.Material`s
=======================================
In reality, most of the Earth is composed of a mixture of different elastic materials. We can approximate the elastic properties of these 'aggregates' using some simple averaging relationships for different proportions of the corresponding stiffness matrices.

For an aggregate of different minerals, Voigt (1928) proposed the averaging—over all possible lattice orientations—of the relations expressing the stress in a single crystal in terms of the given strain, assuming that the strain is uniform throughout the aggregate. The Voigt average elastic moduli for a composite material comprised of *i* different minerals are described mathematically by

.. math:: M_{V} =  \sum_{i}^{N}f_{i}M_{i},

where :math:`M_{i}` and :math:`f_{i}` are the elastic moduli and fractional proportion for each constitutive mineral, respectively. Conversely, Reuß (1929) proposed the averaging of the relations expressing the strain in terms of the given stress, assuming that the stress is uniform throughout an aggregate. Thus, the Reuß average elastic moduli for a composite material comprised of *i* different minerals are described mathematically by

.. math:: \frac{1}{M_{R}} =  \sum_{i}^{N}f_{i}\frac{1}{M_{i}},

where again :math:`M_{i}` and :math:`f_{i}` are the elastic moduli and fractional proportion for each constitutive mineral, respectively. Hill (1952) subsequently demonstrated, by appealing to conservation of energy arguments, that these to aggregate averages are in fact the maximum and minimum bounds on the elastic moduli of the composite material, respectively. In the same work, they proposed the arithmetic mean of the two, :math:`M_{VRH} = \frac{M_{V} + M_{R}}{2}`, as an approximation that more closely matched the experimental observations, known as the Voigt-Reuß-Hill average.

Voigt, Reuß, and Voigt-Reuß-Hill averaging in AnisotroPy
--------------------------------------------------------
AnisotroPy provides functions for calculating the effective elastic tensors for aggregates. These functions accept lists of `anisotropy.materials.Material` objects and the corresponding volume fractions. This can be performed as simply

::

    from anisotropy.materials import load, Material, voigt_reuss_hill_average

    material1 = load("olivine")
    material2 = load("clinopyroxene_98")

    C_vrh, rho_agg, C_voigt, C_reuss = voigt_reuss_hill_average(
        [material1, material2],
        [0.7, 0.3]    
    )

    material_vrh = Material(C_vrh, rho_agg)

Aligning and watering down single-crystal tensors
-------------------------------------------------
Initially aligned against blah. Macroscopic anisotropy often arises due to *preferential* alignment, not full alignment. Hence, anisotropy observed are effectively lower than theoretical values. To approximate these two features, we can rotate the tensor to lie along a certain axis, and water down a single-crystal elastic tensor with a greater proportion of its isotropic component, using the ideas set out by Browaeys and Chevrot (2004).

The elastic tensor of any material can be decomposed into its isotropic, hexagonal, tetragonal, orthorhombic, or monoclinic component via a property of the Material class:

::

    from anisotropy.materials import load

    olivine = load("olivine")

    isotropic = olivine.isotropic
    hexagonal = olivine.hexagonal
    tetragonal = olivine.tetragonal
    orthorhombic = olivine.orthorhombic
    monoclinic = olivine.monoclinic

These are `Material` objects. The isotropic portion, for example, can then be combined with the original olivine using the averaging methods outlined above to produce an aggregate material that represents a preferentially aligned olivine mass. For example:

::

    from anisotropy.materials import load, Material

    olivine = load("olivine")
    isotropic = olivine.isotropic

    C_vrh, rho_agg, *_ = voigt_reuss_hill_average(
        [isotropic, olivine],
        [0.7, 0.3]    
    )
    bulk_medium = Material(C_vrh, rho_agg)

References
----------
Browaeys, J.T. and Chevrot, S., 2004. Decomposition of the elastic tensor and geophysical applications. *Geophysical Journal International, 159(2*)*, pp.667-678.

Hill, R., 1952. The elastic behaviour of a crystalline aggregate. *Proceedings of the Physical Society. Section A, 65(5)*, p.349. (and references therein)
