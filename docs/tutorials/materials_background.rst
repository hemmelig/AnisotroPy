General background on materials and elasticity
==============================================
.. note:: This *précis* covers the fundamentals of materials, deformation, and the elastic equations. For a more complete discussion, read Karato (), Aki and Richards (1998), and Marsden and Hughes (1983)

In this tutorial, we will cover the basics of the Material class. This class performs the task of representing an arbitrarily (an)isotropic elastic material and provides methods for computing a number of elastic moduli e.g. the bulk and shear moduli.

A basic material is defined by an elastic stiffness tensor, :math:`C_{ijkl}` (a 4th-rank tensor), and a density, :math:`\rho`.

A material can also be provided with an assortment of metadata, including: a material ID, and a reference + DOI (Digital Object Identifier).

Stress, strain, and elasticity
------------------------------
The concepts of stress and strain are key to understanding deformation, and in turn seismic anisotropy.

Wave propagation in arbitrarily (an)isotropic elastic media
-----------------------------------------------------------
.. note:: For detailed discussions of the theory of seismic wave propagation in anisotropic media, on which this material is based, see e.g. Babuška and Cara (1991) or Silver (1996).

Within the *linear* elastic theory of continuum mechanics, seismic anisotropy enters the elastodynamic equations of motion through the fourth-order stiffness tensor, composed of 81 independent elastic moduli. This tensor formally relates the applied stress, :math:`\sigma_{ij}`, to the resulting deformation of an elastic body, :math:`\epsilon_{kl}`, via Hooke's law

.. math:: \sigma_{ij} = C_{ijkl}\epsilon_{kl}.

The 81 independent elastic moduli can be reduced to 36 through a number of symmetry arguments. Firstly,

.. math:: C_{ijkl} = C_{jikl},

due to the symmetry of :math:`\sigma_{ij}`. Secondly,

.. math:: C_{ijkl} = C_{ijlk},

due to the symmetry of :math:`\epsilon_{kl}`. And finally, these 36 independent moduli can be reduced further by the constraint

.. math:: C_{ijkl} = C_{klij},

which can be established from thermodynamic arguments (e.g. Nye, 1985). Assuming a perfectly elastic material, all of the energy put in while deforming it can be recovered by allowing the material to return to its equilibrium position. This means that the integral of the strain energy density associated with a given strain is path independent, which is only true if the above equality holds true. Overall, the number of independent elastic moduli is reduced from 81 to 21. Hence, the general 3x3x3x3 stiffness tensor can be expressed in a contracted "Voigt" form as the symmetric 6x6 matrix

.. math:: c_{mn} = 
  \begin{bmatrix}
  C_{1111} & C_{1122} & C_{1133} & C_{1123} & C_{1131} & C_{1112} \\ 
  C_{2211} & C_{2222} & C_{2233} & C_{2223} & C_{2231} & C_{2212} \\ 
  C_{3311} & C_{3322} & C_{3333} & C_{3323} & C_{3331} & C_{3312} \\ 
  C_{2311} & C_{2322} & C_{2333} & C_{2323} & C_{2331} & C_{2312} \\ 
  C_{3111} & C_{3122} & C_{3133} & C_{3123} & C_{3131} & C_{3112} \\ 
  C_{1211} & C_{1222} & C_{1233} & C_{1223} & C_{1231} & C_{1212} \\ 
  \end{bmatrix}.

Symmetry systems
----------------
Isotropic
#########
An isotropic solid is a simple hyperelastic material whose symmetry group at each point is equal to the full special orthogonal group :math:`SO(3)` (i.e. the 3-D rotation group). This means that stretching of the material about any given point in any orientation (equivalent to stretching the material along a fixed axis *after* performing an arbitrary 3-D rotation) leads to the same strain energy being produced. From this property it can be shown that the elastic tensor relative to a stress free equilibrium configuration takes the form

.. math:: C_{ijkl} = \lambda\delta_{ij}\delta_{kl} + \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{kj}),

where :math:`\lambda` and :math:`\mu` are the Lamé coefficients. The corresponding Voigt-contracted form stiffness matrix is

.. math:: c_{mn} = 
  \begin{bmatrix}
  \lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\ 
  \lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\ 
  \lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\ 
  0 & 0 & 0 & \mu & 0 & 0 \\ 
  0 & 0 & 0 & 0 & \mu & 0 \\ 
  0 & 0 & 0 & 0 & 0 & \mu \\ 
  \end{bmatrix}.

Creating a new `Material`
-------------------------
To create a new `Material` object, you need to provide a stiffness tensor in the reduced Voigt form. Due to the symmetry constraints, outlined briefly above

References
----------
Babuška, V. and Cara, M. (1991). *Seismic anisotropy in the Earth*, volume 10. Springer Science & Business Media.

Karato

Silver, P. G. (1996). Seismic anisotropy beneath the continents: Probing the depths of geology. *Annual Review of Earth and Planetary Sciences*, 24(1):385–432.