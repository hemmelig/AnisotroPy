# -*- coding: utf-8 -*-
"""
Module for performing all things related to materials and their elastic
properties. This includes the elastic stiffness and compliance tensors, various
moduli, etc.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

import anisotropy.utils as utils


C_iso = np.array([[237.5533,  78.4733,  78.4733,  0.0000,  0.0000,  0.0000],
                  [ 78.4733, 237.5533,  78.4733,  0.0000,  0.0000,  0.0000],
                  [ 78.4733,  78.4733, 237.5533,  0.0000,  0.0000,  0.0000],
                  [  0.0000,   0.0000,   0.0000, 79.5400,  0.0000,  0.0000],
                  [  0.0000,   0.0000,   0.0000,  0.0000, 79.5400,  0.0000],
                  [  0.0000,   0.0000,   0.0000,  0.0000,  0.0000, 79.5400]])

class Material:
    """
    Small class to encapsulate information that defines a material with an
    elastic stiffness tensor, C, and density, rho, and a suite of useful
    functions to query its properties.

    Within the linear elastic theory of continuum mechanics, seismic anisotropy
    enters the elastodynamic equations of motion through the fourth-order
    stiffness tensor, composed of 81 independent elastic moduli. This tensor
    formally relates the applied stress, σ_ij , to the resulting deformation of
    an elastic body, ε_kl, via Hooke’s law.

    Attributes
    ----------
    C : array-like of float, shape(6, 6)
        Stiffness tensor for the material, in GPa.
    rho : float
        Density of material, in g/cm^3.
    is_isotropic : bool
        Tests whether the material is isotropic.
    lame_coefficients : [float, float]
        Calculates Lame coefficients for the material.
    bulk_modulus : float
        Calculates the bulk modulus, K, for the material.

    Methods
    -------
    phase_velocities
        Calculates the phase velocities for the material at requested azimuth/
        inclination pairs.
    group_velocities
        Calculates the group velocities, which are the phase velocities
        projected onto the group velocity directions, for the material at
        requested azimuth/inclination pairs.
    rotate
        Arbitrarily rotates the materials stiffness tensor in 3 dimensions.

    TO-DO
    -----
    - Provide methods to transform between the Voigt-form and full tensor
      representations of the elasticity tensors
    - Finish documentation
    - Implement method to calculate the isotropic (symmetric) portion of a
      arbitrarily anisotropic stiffness tensor

    """

    _VOIGT_CONTRACTION_MATRIX = np.array([[0, 5, 4],
                                          [5, 1, 3],
                                          [4, 3, 2]])
    
    isotropic = C_iso

    def __init__(self, C, rho):
        """Instantiate the Material object."""

        self.C = C
        self.rho = rho

    def __str__(self):
        """Pretty string representation of the materials stiffness tensor, C."""

        C = self.C
        str_ = ("C_ijkl =\n\n"
                f"  {C[0,0]: 8.3f} {C[0, 1]: 8.3f} {C[0, 2]: 8.3f} "
                f"{C[0, 3]: 8.3f} {C[0, 4]: 8.3f} {C[0, 5]: 8.3f}\n"
                f"  {C[1,0]: 8.3f} {C[1, 1]: 8.3f} {C[1, 2]: 8.3f} "
                f"{C[1, 3]: 8.3f} {C[1, 4]: 8.3f} {C[1, 5]: 8.3f}\n"
                f"  {C[2,0]: 8.3f} {C[2, 1]: 8.3f} {C[2, 2]: 8.3f} "
                f"{C[2, 3]: 8.3f} {C[2, 4]: 8.3f} {C[2, 5]: 8.3f}\n"
                f"  {C[3,0]: 8.3f} {C[3, 1]: 8.3f} {C[3, 2]: 8.3f} "
                f"{C[3, 3]: 8.3f} {C[3, 4]: 8.3f} {C[3, 5]: 8.3f}\n"
                f"  {C[4,0]: 8.3f} {C[4, 1]: 8.3f} {C[4, 2]: 8.3f} "
                f"{C[4, 3]: 8.3f} {C[4, 4]: 8.3f} {C[4, 5]: 8.3f}\n"
                f"  {C[5,0]: 8.3f} {C[5, 1]: 8.3f} {C[5, 2]: 8.3f} "
                f"{C[5, 3]: 8.3f} {C[5, 4]: 8.3f} {C[5, 5]: 8.3f}\n"
                f"\nDensity = {self.rho}")

        return str_

    def phase_velocities(self, inclination, azimuth=np.arange(0, 361, 1)):
        """
        Calculate the phase velocities for the 6x6 elasticity matrix C along
        the direction (ψ, θ), in degrees, and return P-wave velocity Vp, the
        fast and slow shear wave velocities, Vs1 and Vs2, the polarisation of
        the fast shear wave, pol, and the shear wave velocity anisotropy, swa.

        Parameters
        ----------
        inclination : list of float or float
            Inclination, θ, from the x1-x2 plane towards x3, in degrees.
        azimuth : list of float or float, optional
            Azimuth, ψ, from x1 towards x2, in degrees.

        Returns
        -------
        vp : list of float
            P-wave phase velocity in km/s.
        vs1 : list of float
            Faster S-wave phase velocity in km/s.
        vs2 : list of float
            Slower S-wave phase velocity in km/s.
        φφ : list of float
            Polarisation direction of the faster S-wave.
        shear_wave_anisotropy : list of float
            Calculated shear wave anisotropy.

        Raises
        ------
        ValueError
            When vectors provided are not of equal length.

        """

        # Run tests to see if the angles are scalars or lists
        if np.isscalar(azimuth):
            azimuth = [azimuth]
        if np.isscalar(inclination):
            inclination = np.ones(len(azimuth))*inclination
        if len(azimuth) != len(inclination):
            print("Must provide vectors of azimuth, ψ, and inclination, θ, of "
                  "equal length.")
            raise TypeError

        vp, vs1, vs2, φφ = [], [], [], []
        for ψ, θ in zip(azimuth, inclination):
            x = utils.azinc2vec(ψ, θ)

            # Make the Christoffel matrix
            M = self._christoffel(x)

            # Determine eigenvalues and eigenvectors corresponding to the
            # velocities and polarisation vectors.
            eigenvalues, eigenvectors = np.linalg.eig(M)
            ip = np.argmax(eigenvalues)
            is2 = np.argmin(eigenvalues)
            is1 = 3 - ip - is2
            phase_velocities = np.sqrt(eigenvalues/self.rho)

            vp.append(phase_velocities[ip])
            vs1.append(phase_velocities[is1])
            vs2.append(phase_velocities[is2])
            φφ.append(
                self._calculate_polarisation(ψ, θ, x, eigenvectors[:, is1]))

        shear_wave_anisotropy = 200*(np.subtract(vs1, vs2)/np.add(vs1, vs2))

        return vp, vs1, vs2, φφ, shear_wave_anisotropy

    def group_velocities(self, inclination, azimuth=np.arange(0, 361, 1)):
        """
        Calculate the group velocities for the 6x6 elasticity matrix C along
        the direction (ψ, θ), in degrees, and return P-wave velocity Vp, and
        the fast and slow shear wave velocities, Vs1 and Vs2.

        !! CURRENTLY APPEARS TO NOT BE WORKING !!

        Parameters
        ----------
        inclination : list of float or float
            Inclination, θ, from the x1-x2 plane towards x3, in degrees.
        azimuth : list of float or float, optional
            Azimuth, ψ, from x1 towards x2, in degrees.

        Returns
        -------
        gvp : list of float
            P-wave group velocity in km/s.
        gvs1 : list of float
            Faster S-wave group velocity in km/s.
        gvs2 : list of float
            Slower S-wave group velocity in km/s.

        """

        # Run tests to see if the angles are scalars or lists
        if np.isscalar(azimuth):
            azimuth = [azimuth]
        if np.isscalar(inclination):
            inclination = np.ones(len(azimuth))*inclination
        if len(azimuth) != len(inclination):
            print("Must provide vectors of azimuth, ψ, and inclination, θ, of "
                  "equal length.")
            raise TypeError

        gvp, gvs1, gvs2 = [], [], []
        for ψ, θ in zip(azimuth, inclination):
            x = utils.azinc2vec(ψ, θ)

            # Make the Christoffel matrix
            M = self._christoffel(x)

            # Determine eigenvalues and eigenvectors corresponding to the
            # velocities and polarisation vectors.
            eigenvalues, eigenvectors = np.linalg.eig(M)
            ip = np.argmax(eigenvalues)
            is2 = np.argmin(eigenvalues)
            is1 = 3 - ip - is2
            phase_velocities = np.sqrt(eigenvalues/self.rho)

            # Calculate group velocities, which are phase velocities projected
            # onto group velocity directions.
            vp, xp = phase_velocities[ip], eigenvectors[:, ip]
            pp = xp / vp
            vs1, xs1 = phase_velocities[is1], eigenvectors[:, is1]
            ps1 = xs1 / vs1
            vs2, xs2 = phase_velocities[is2], eigenvectors[:, is2]
            ps2 = xs2 / vs2
            p = np.array([pp, ps1, ps2]).T

            group_v = np.zeros(3)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            m = self._VOIGT_CONTRACTION_MATRIX[i, j]
                            n = self._VOIGT_CONTRACTION_MATRIX[k, l]
                            group_v[i] += self.C[m, n]*p[i, k]*x[j]*x[l]

            gvp.append(group_v[0])
            gvs1.append(group_v[1])
            gvs2.append(group_v[2])

        return gvp, gvs1, gvs2

    def rotate(self, alpha, beta, gamma, order=[1, 2, 3], mode="extrinsic"):
        """Helper function to perform an arbitrary 3-D rotation of C."""
        R = self._rotate_3d(alpha, beta, gamma, order, mode)
        Q = self._build_bond_matrix(R)

        return np.dot(Q, np.dot(self.C, Q.T))

    def _rotate_3d(self, alpha, beta, gamma, order=[1, 2, 3], mode="extrinsic"):
        """
        Contructs a 3x3 matrix that specifies a 3-D rotation from the Euler
        angles alpha, beta, and gamma.

        Parameters
        ----------
        alpha : float
            The yaw of the rotation (i.e. azimuth).
        beta : float
            The pitch of the rotation (i.e. the inclination)
        gamma : float
            The roll of the rotation (i.e. the wibbly wobbly side to sidey)
        order : list of int
            Defines the order of rotation. Important, since rotations are not
            generally commutative.
        mode : {"instrinsic", "extrinsic"}
            Form of rotation.

        Returns
        -------
        R : array-like of float, shape(3, 3)
            Composite rotation matrix.

        """

        if mode == "intrinsic":
            r_z, r_y, r_x = np.deg2rad([alpha, beta, gamma])
        elif mode == "extrinsic":
            r_z, r_y, r_x = np.deg2rad([gamma, beta, alpha])

        R_z = np.array([[ np.cos(r_z),  np.sin(r_z), 0],
                        [-np.sin(r_z),  np.cos(r_z), 0],
                        [           0,            0, 1]])

        R_y = np.array([[np.cos(r_y), 0, -np.sin(r_y)],
                        [          0, 1,            0],
                        [np.sin(r_y), 0,  np.cos(r_y)]])

        R_x = np.array([[1,            0,           0],
                        [0,  np.cos(r_x), np.sin(r_x)],
                        [0, -np.sin(r_x), np.cos(r_x)]])
        Rs = [R_x, R_y, R_z]

        R = np.dot(Rs[order[2]-1], np.dot(Rs[order[1]-1], Rs[order[0]-1]))

        return R
        
    def _build_bond_matrix(self, R):
        """Construct Bond-like 6x6 representation of a 3x3 rotation matrix."""
        return np.array(
            [[    R[0,0]**2,     R[0,1]**2,     R[0,2]**2,             2*R[0,1]*R[0,2],             2*R[0,0]*R[0,2],             2*R[0,0]*R[0,1]],
             [    R[1,0]**2,     R[1,1]**2,     R[1,2]**2,             2*R[1,1]*R[1,2],             2*R[1,0]*R[1,2],             2*R[1,0]*R[1,1]],
             [    R[2,0]**2,     R[2,1]**2,     R[2,2]**2,             2*R[2,1]*R[2,2],             2*R[2,0]*R[2,2],             2*R[2,0]*R[2,1]],
             [R[1,0]*R[2,0], R[1,1]*R[2,1], R[1,2]*R[2,2], R[1,1]*R[2,2]+R[1,2]*R[2,1], R[1,0]*R[2,2]+R[1,2]*R[2,0], R[1,0]*R[2,1]+R[1,1]*R[2,0]],
             [R[0,0]*R[2,0], R[0,1]*R[2,1], R[0,2]*R[2,2], R[0,1]*R[2,2]+R[0,2]*R[2,1], R[0,0]*R[2,2]+R[0,2]*R[2,0], R[0,0]*R[2,1]+R[0,1]*R[2,0]],
             [R[0,0]*R[1,0], R[0,1]*R[1,1], R[0,2]*R[1,2], R[0,1]*R[1,2]+R[0,2]*R[1,1], R[0,0]*R[1,2]+R[0,2]*R[1,0], R[0,0]*R[1,1]+R[0,1]*R[1,0]]])

    def _christoffel(self, x):
        """
        Calculate the 3x3 matrix, M, that appears in the Christoffel equation
        for a plane wave in a generally (an)isotropic material.

        Parameters
        ----------
        x : array of floats
            The 3-vector normal to the slowness surface, i.e. parallel to the
            direction of propagation.

        Returns
        -------
        M : 3x3 array of floats
            The Christoffel matrix, M.

        """

        M = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        m = self._VOIGT_CONTRACTION_MATRIX[i, j]
                        n = self._VOIGT_CONTRACTION_MATRIX[k, l]
                        M[i, k] += self.C[m, n]*x[j]*x[l]

        return M

    def _calculate_polarisation(self, inclination, azimuth, x1, x2):
        """
        Calculates the projection of the fast shear wave polarisation onto the
        wavefront plane.

        Parameters
        ----------
        inclination : list of float or float
            Inclination, θ, from the x1-x2 plane towards x3, in degrees.
        azimuth : list of float or float, optional
            Azimuth, ψ, from x1 towards x2, in degrees.
        x1 : list of float
            Vector pointing along (ψ, θ).
        x2 : list of float
            Polarisation vector of the fast shear wave.

        Returns
        -------
        φ : float
            Projection of the polarisation of the fast shear wave onto the
            wavefront plane.

        """

        plane = np.cross(x1, np.cross(x1, x2))
        plane /= np.linalg.norm(plane)

        rotated_vec = utils.rotate2xy(plane.copy(), azimuth, inclination)

        φ = np.rad2deg(np.arctan2(rotated_vec[1], rotated_vec[2]))
        φ = φ + 180 if φ < -90 else φ
        φ = φ - 180 if φ > 90 else φ

        return φ

    @property
    def is_isotropic(self):
        """Tests whether the material is isotropic. Returns a boolean."""
        # TO BE IMPLEMENTED
        return True

    @property
    def lame_coefficients(self):
        """Returns the Lame coefficients for the material."""
        # Material must be isotropic for these to be meaningful!
        if not self.is_isotropic:
            raise ValueError
        mu = self.C[3, 3]
        la = self.C[0, 0] - 2*mu

        return la, mu

    @property
    def bulk_modulus(self):
        """Returns the bulk modulus for the material."""
        # Material must be isotropic for this to be meaningful!
        if not self.is_isotropic:
            raise ValueError
        la, mu = self.lame_coefficients

        return la + 2*mu/3

def voigt_average(stiffness_tensors, volume_fractions):
    """
    Calculates the Voigt average of a set material phases, described by their
    elastic stiffness tensors and densities, and constitutive phase volume
    fractions.

    The Voigt average is defined as the weighted arithmetic mean of elastic
    stiffness tensors, that is:

                    C_v = f1*C_1 + f2*C_2 + ... + fn*C_n

    For more information, read:

        Hill, R., 1952. The elastic behaviour of a crystalline aggregate.
        Proceedings of the Physical Society. Section A, 65(5), p.349.

    Parameters
    ----------
    stiffness_tensors : list of array-like of float, shape(n, 6, 6)
        Stiffness tensors of constitutive phases for which to calculate the
        Reuss average.
    volume_fractions : list of float
        Fraction of each constitutive phase in the composite material.

    Returns
    -------
    C_v : array-like of float, shape(6, 6)
        Voigt average elastic stiffness tensor, C, for the composite material.

    """

    # Validate volume fraction sums to unity
    if np.sum(volume_fractions) != 1:
        raise ValueError("Volume fractions must sum to unity!")

    return np.add(*[C*volume_fraction
                    for C, volume_fraction
                    in zip(stiffness_tensors, volume_fractions)])


def reuss_average(stiffness_tensors, volume_fractions):
    """
    Calculates the Reuss average of a set material phases, described by their
    elastic stiffness tensors and densities, and constitutive phase volume
    fractions.

    The Reuss average is defined as the inverse of the weighted arithmetic mean
    of elastic compliance tensors (i.e. the inverse of the elastic stiffness
    tensors), that is:

                    C_r = __________________1_________________
                          f1*1/C_1 + f2*1/C_2 + ... + fn*1/C_n

    For more information, read:

        Hill, R., 1952. The elastic behaviour of a crystalline aggregate.
        Proceedings of the Physical Society. Section A, 65(5), p.349.

    Parameters
    ----------
    stiffness_tensors : list of array-like of float, shape(n, 6, 6)
        Stiffness tensors of constitutive phases for which to calculate the
        Reuss average.
    volume_fractions : list of float
        Fraction of each constitutive phase in the composite material.

    Returns
    -------
    C_r : array-like of float, shape(6, 6)
        Reuss average elastic stiffness tensor, C, for the composite material.

    """

    # Validate volume fraction sums to unity
    if np.sum(volume_fractions) != 1:
        raise ValueError("Volume fractions must sum to unity!")

    C_r = np.add(*[np.linalg.inv(C)*volume_fraction
                   for C, volume_fraction
                   in zip(stiffness_tensors, volume_fractions)])

    return np.linalg.inv(C_r)


def voigt_reuss_hill_average(materials, volume_fractions):
    """
    Calculates the Voigt-Reuss-Hill average of a set material phases, described
    by their elastic stiffness tensors and densities, and constitutive phase
    volume fractions.

    This was defined as the arithmetic mean of the Voigt average (arithmetic
    mean of stiffness tensors weighted by volume fraction) and the Reuss
    average (arithmetic mean of compliance tensors weighted by volume
    fraction).

    For more information, read:

        Hill, R., 1952. The elastic behaviour of a crystalline aggregate.
        Proceedings of the Physical Society. Section A, 65(5), p.349.

    Parameters
    ----------
    materials : list of `anisotropy.Material`
        Constitutive phases for which to calculate the Voigt-Reuss-Hill average.
    volume_fractions : list of float
        Fraction of each constitutive phase in the composite material.

    Returns
    -------
    C_vrh : array-like of float, shape(6, 6)
        Voigt-Reuss-Hill average elastic stiffness tensor, C, for the composite
        material.
    rho_average : float
        Density of composite material.
    C_v : array-like of float, shape(6, 6)
        Voigt average elastic stiffness tensor, C, for the composite material.
    C_r : array-like of float, shape(6, 6)
        Reuss average elastic stiffness tensor, C, for the composite material.

    """

    # Validate volume fraction sums to unity
    if np.sum(volume_fractions) != 1:
        raise ValueError("Volume fractions must sum to unity!")

    C_v = voigt_average([material.C for material in materials],
                         volume_fractions)
    C_r = reuss_average([material.C for material in materials],
                         volume_fractions)
    C_vrh = (C_v + C_r) / 2
    rho_vrh = np.average([material.rho for material in materials],
                          weights=volume_fractions)

    return C_vrh, rho_vrh, C_v, C_r


def isotropic_C(vp=None, vs=None, rho=None, la=None, mu=None, K=None, G=None):
    """
    Calculate the coefficients of the stiffness tensor expressed in the Voigt
    matrix form.

    Voigt's form does not preserve the sum of the squares, which in the case of
    Hooke's law has geometric significance.

    Parameters
    ----------
    vp : float, optional
        Isotropic P wave velocity. Units of km/s.
    vs : float, optional
        Isotropic S wave velocity. Units of km/s.
    rho : float, optional
        Bulk density of material. Units of g/cm^3
    la : float, optional
        Lamé's first parameter. Units of m^2/s^2.
    mu : float, optional
        Lamé's second parameter - the shear modulus.
    K : float, optional
        Bulk modulus.
    G : float, optional
        Shear modulus.

    Returns
    -------
    C : 6x6 array of floats
        Isotropic elastic stiffness tensor in Voigt matrix form.

    """

    C = np.zeros((6, 6))

    if vp is not None and vs is not None:
        C[0, 0] = vp**2
        C[3, 3] = vs**2
        C[0, 1] = C[0, 0] - 2*C[3, 3]
        if rho is not None:
            C *= rho
        else:
            C *= _vp2rho(vp)
    elif la is not None and mu is not None:
        C[0, 0] = la + 2*mu
        C[3, 3] = mu
        C[0, 1] = la
    elif K is not None and G is not None:
        C[0, 0] = K + 4*G/3
        C[3, 3] = G
        C[0, 1] = C[0, 0] - 2*C[3, 3]
    else:
        print("No arguments provided.")

    C[1, 1] = C[2, 2] = C[0, 0]
    C[4, 4] = C[5, 5] = C[3, 3]
    C[1, 0] = C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C[0, 1]

    return C


def _vp2rho(vp):
    """
    Convert from Vp to density using Birch's law:

    rho = 0.32*Vp + 0.77

    Parameters
    ----------
    vp : float
        P wave velocity.

    Returns
    -------
    rho : float
        Density using empirical law.

    """

    return 0.32*vp + 0.77
