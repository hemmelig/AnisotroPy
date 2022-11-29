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

import itertools

import numpy as np

import anisotropy.utils as utils
import anisotropy.utils.errors as errors


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
        Calculates Lamé coefficients for the material.
    bulk_modulus : float
        Calculates the bulk modulus, K, for the material.
    shear_modulus : float
        Calculates the shear modulus, G, for the material.

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
    - Finish documentation
    - Implement method to calculate the isotropic (symmetric) portion of a
      arbitrarily anisotropic stiffness tensor
    - Remove olivine specific isotropic component

    """

    _VOIGT_CONTRACTION_MATRIX = np.array([[0, 5, 4],
                                          [5, 1, 3],
                                          [4, 3, 2]])
    
    isotropic = C_iso

    def __init__(self, C, rho, material_id="", reference=""):
        """Instantiate the Material object."""

        self.C = C
        self.rho = rho
        self.id = material_id
        self.reference = reference

    def __str__(self):
        """Pretty string representation of the materials stiffness tensor, C."""

        C = self.C
        str_ = (
            f"Material ID - {self.id}\n\n"
            "C_ijkl =\n\n"
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
            f"\nDensity (g/cm^3) = {self.rho}\n\n"
            "Reference\n"
            f"{self.reference}"
        )

        return str_

    def __eq__(self, other):
        """Magic method for rich comparison '=='"""

        if not isinstance(other, Material):
            raise TypeError(
                "Comparison only valid between anisotropy.materials.Material "
                "objects."
            )

        pred1 = np.allclose(self.C, other.C)
        pred2 = self.rho == other.rho

        return bool(pred1 and pred2)

    def phase_velocities(self, inclination, azimuth=None):
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
        elif azimuth is None:
            azimuth = np.arange(0, 361, 1)
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
            # abs to ensure real velocity
            phase_velocities = abs(np.sqrt(eigenvalues/self.rho))

            vp.append(phase_velocities[ip])
            vs1.append(phase_velocities[is1])
            vs2.append(phase_velocities[is2])
            φφ.append(
                self._calculate_polarisation(θ, ψ, x, eigenvectors[:, is1]))

        shear_wave_anisotropy = 200*(np.subtract(vs1, vs2)/np.add(vs1, vs2))

        return vp, vs1, vs2, φφ, shear_wave_anisotropy

    def group_velocities(self, inclination, azimuth=None):
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
        elif azimuth is None:
            azimuth = np.arange(0, 361, 1)
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
        """
        Arbitrary 3-D rotation of C.

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
        """
        R = self._rotate_3d(alpha, beta, gamma, order, mode)
        Q = self._build_bond_matrix(R)
        self.C = np.dot(Q, np.dot(self.C, Q.T))

        return self

    def plot_velocity(
        self,
        which="vp",
        incs=[],
        azis=[],
        incr=(None, None),
        minmax=True,
        cmap="inferno_r",
        ax=None,
        fig_kw={"figsize": (3, 3), "constrained_layout": True},
    ):
        """
        Plot seismic velocity on stereonet

        Parameters:
        -----------
            which: (str)
                which velocity to plot: vp, vs1, vs2, vpvs1, or vpvs2
            incs: (list)
                Limited inclinations (degree from x1-x2-plane to x3) to plot
            azis: (list)
                Limited azimuths (degree up from x1 to x2) to plot
            incr: (2*tuple)
                Increments in inclination and azimuth
            minmax: (bool)
                Show minimum / maximum values in plot
            cmap: (str)
                Color map to use to fraw colors
            ax: (matpltlib.Axes)
                Draw axis in external figure
            fig_kw: (dict)
                kwargs passed to created figure (ignored when ax is used)
        """

        # import matplotlib.pyplot as mp
        import mplstereonet as ms
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        def _xys(mr, ir, maz, iaz):
            """
            Return recatangle corners from midpoint (m) and increment (i) in radial
            (r) and azimuthal (az) direction
            """
            x1 = maz - iaz / 2
            x2 = maz + iaz / 2
            y1 = mr - ir / 2
            y2 = mr + ir / 2
            xs = [x1, x1, x2, x2]
            ys = [y1, y2, y2, y1]
            return xs, ys

        def _inc(arr):
            """Smallest increment within array"""
            try:
                ad = abs(np.diff(arr))
                return min(ad[ad > 0])
            except ValueError:
                print("Could not determine increment of fixed point spacing")
                return 0

        _label = {
            "vp": "$V_P$ (km/s)",
            "vs1": "$V_{{{S1}}}$ (km/s)",
            "vs2": "$V_{{{S2}}}$ (km/s)",
            "vpvs1": "$V_P/V_{{{S1}}}$",
            "vpvs2": "$V_P/V_{{{S2}}}$",
        }

        docont = False
        if not any(incs) or not any(azis):
            # Evenly sample sphere for triangulation
            # https://stackoverflow.com/questions/
            # 9600801/evenly-distributing-n-points-on-a-sphere
            n = 500
            idx = np.arange(0, n, dtype=float) + 0.5
            incs = ((np.pi - np.arccos(idx / n)) * 180 / np.pi) % 90
            azis = ((np.pi * (1 + 5**0.5) * idx) * 180 / np.pi) % 360

            # Circumfer for nicer triangulation
            azis = np.append(azis, np.linspace(0, 360, n // 10))
            incs = np.append(incs, np.zeros(n // 10))
            docont = True
        else:
            incs = np.array(incs)
            azis = np.array(azis)

        lons, lats = ms.line(incs, azis)

        vp, vs1, vs2, _, _ = self.phase_velocities(incs, azis)

        if which == "vp":
            dats = vp
        elif which == "vs1":
            dats = vs1
        elif which == "vs2":
            dats = vs2
        elif which == "vpvs1":
            dats = np.divide(vp, vs1)
        elif which == "vpvs2":
            dats = np.divide(vp, vs2)
        else:
            msg = f"Unknown 'which': {which}"
            raise ValueError(msg)

        imin = np.argmin(dats)
        imax = np.argmax(dats)

        # Get the colors right
        norm = Normalize(vmin=min(dats), vmax=max(dats))
        mapper = ScalarMappable(norm=norm, cmap=cmap)

        # Do the plot
        if ax:
            # Make sure it's a stereonet
            fig = ax.get_figure()
            spec = ax.get_subplotspec()
            ax.set_axis_off()
            del ax
            ax = fig.add_subplot(spec, projection="stereonet")
        else:
            fig, ax = ms.subplots(1, 1, **fig_kw)

        ax.set_azimuth_ticks([])

        if docont:
            # contour the half sphere
            ax.tricontourf(lons, lats, dats, norm=norm, cmap=cmap)
        else:
            # fill only sampled segments
            # find angle increments
            dinc = incr[0]
            if not dinc:
                dinc = _inc(incs)
            daz = incr[1]
            if not daz:
                daz = _inc(azis)

            for inc, az, dat in zip(incs, azis, dats):
                xs, ys = _xys(inc, dinc, az, daz)
                lon, lat = ms.line(ys, xs)
                ax.fill(lon, lat, color=mapper.cmap(norm(dat)))

        # Write out minimum and maximum
        ax.plot(lons[imin], lats[imin], "o", mec="black", mfc="none")
        ax.plot(lons[imax], lats[imax], "x", mec="black", mfc="none")

        cb = fig.colorbar(mapper, ax=ax, orientation="horizontal", label=_label[which])
        ticks = [dats[imin], dats[imin] + (dats[imax] - dats[imin])/2, dats[imax]]
        cb.set_ticks(ticks)
        cb.set_ticklabels(["{:.2f}".format(t) for t in ticks])

        return fig

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
        r00, r01, r02 = R[0, :]
        r10, r11, r12 = R[1, :]
        r20, r21, r22 = R[2, :]
        return np.array(
            [[ r00**2,  r01**2,  r02**2,       2*r01*r02,       2*r00*r02,       2*r00*r01],
             [ r10**2,  r11**2,  r12**2,       2*r11*r12,       2*r10*r12,       2*r10*r11],
             [ r20**2,  r21**2,  r22**2,       2*r21*r22,       2*r20*r22,       2*r20*r21],
             [r10*r20, r11*r21, r12*r22, r11*r22+r12*r21, r10*r22+r12*r20, r10*r21+r11*r20],
             [r00*r20, r01*r21, r02*r22, r01*r22+r02*r21, r00*r22+r02*r20, r00*r21+r01*r20],
             [r00*r10, r01*r11, r02*r12, r01*r12+r02*r11, r00*r12+r02*r10, r00*r11+r01*r10]])

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

        try:
            φ = np.rad2deg(np.arctan2(rotated_vec[1], rotated_vec[2]))
        except TypeError:
            # rotated_vec may be complex for isotropic media
            φ = np.nan
        φ = φ + 180 if φ < -90 else φ
        φ = φ - 180 if φ > 90 else φ

        return φ

    @property
    def is_isotropic(self):
        """Tests whether the material is isotropic. Returns a boolean."""
        C = self.C
        # Build predicates
        pred1 = bool(C[0, 1] == C[0, 2] == C[1, 2])
        pred2 = bool(C[3, 3] == C[4, 4] == C[5, 5])
        pred3 = bool(C[0, 0] == C[1, 1] == C[2, 2] == C[0, 1] + 2*C[3, 3])

        return bool((pred1 and pred2 and pred3))

    @property
    def isotropic(self):
        """Returns the isotropic component of elastic material."""
        return Material(decompose_C(self)["isotropic"], self.rho)

    @property
    def hexagonal(self):
        """Returns the hexagonal component of elastic material."""
        return Material(decompose_C(self)["hexagonal"], self.rho)

    @property
    def tetragonal(self):
        """Returns the tetragonal component of elastic material."""
        return Material(decompose_C(self)["tetragonal"], self.rho)

    @property
    def orthorhombic(self):
        """Returns the orthorombic component of elastic material."""
        return Material(decompose_C(self)["orthorhombic"], self.rho)

    @property
    def monoclinic(self):
        """Returns the monoclinic component of elastic material."""
        return Material(decompose_C(self)["monoclinic"], self.rho)

    @property
    def lame_coefficients(self):
        """Returns the Lamé coefficients for the material."""
        # Material must be isotropic for these to be meaningful!
        if not self.is_isotropic:
            raise ValueError
        mu = self.C[3, 3]
        la = self.C[0, 0] - 2*mu

        return la, mu

    @property
    def bulk_modulus(self):
        """Returns the bulk modulus for the material."""
        # If the material is isotropic, simply use the Lamé coefficients
        if self.is_isotropic:
            la, mu = self.lame_coefficients
            return la + 2*mu/3

        # Calculate Voigt's bulk modulus
        C = self.C.copy()
        K_voigt = (C[0, 0] + C[1, 1] + C[2, 2]
                   + 2*(C[0, 1] + C[0, 2] + C[1, 2])) / 9

        # Calculate Reuss' bulk modulus
        S = np.linalg.inv(C)
        K_reuss = 1 / (S[0, 0] + S[1, 1] + S[2, 2]
                       + 2*(S[0, 1] + S[0, 2] + S[1, 2]))

        return (K_voigt + K_reuss) / 2

    @property
    def shear_modulus(self):
        """Returns the shear modulus for the material."""
        # If the material is isotropic, simply use the Lamé coefficients
        if self.is_isotropic:
            _, mu = self.lame_coefficients
            return mu

        # Calculate Voigt's shear modulus
        C = self.C.copy()
        G_voigt = (C[0, 0] + C[1, 1] + C[2, 2]
                   - (C[0, 1] + C[1, 2] + C[0, 2])
                   + 3*(C[3, 3] + C[4, 4] + C[5, 5])) / 15

        # Calculate Reuss' shear modulus
        S = np.linalg.inv(C)
        G_reuss = 15 / (4*(S[0, 0] + S[1, 1] + S[2, 2])
                        - 4*(S[0, 1] + S[0, 2] + S[1, 2])
                        + 3*(S[3, 3] + S[4, 4] + S[5, 5]))

        return (G_voigt + G_reuss) / 2

    @property
    def C_tensor(self):
        """Returns the elastic stiffness tensor in full 4th rank tensor form."""
        tensor = np.zeros((3, 3, 3, 3), dtype=float)
        for i, j, k, l in itertools.product(range(3), repeat=4):
            m = self._VOIGT_CONTRACTION_MATRIX[i, j]
            n = self._VOIGT_CONTRACTION_MATRIX[k, l]
            tensor[i, j, k, l] = self.C[m, n]

        return tensor


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

    if rho is None:
        if vp is not None:
            rho = _vp2rho(vp)
            # Emit a warning here
        else:
            raise errors.MissingDensityValue

    if vp is not None and vs is not None:
        C[0, 0] = vp**2
        C[3, 3] = vs**2
        C[0, 1] = vp**2 - 2*vs**2
        C *= rho
    elif la is not None and mu is not None:
        C[0, 0] = la + 2*mu
        C[3, 3] = mu
        C[0, 1] = la
    elif K is not None and G is not None:
        C[0, 0] = K + 4*G/3
        C[3, 3] = G
        C[0, 1] = K - 2*G/3
    else:
        raise errors.InsufficientElasticInformation

    C[1, 1] = C[2, 2] = C[0, 0]
    C[4, 4] = C[5, 5] = C[3, 3]
    C[1, 0] = C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C[0, 1]

    return Material(C, rho, material_id="isotropic material")


def hexagonal1_C(vp, vs, ani, rho):
    """
    Voigt's representation of simplified 1-parameter hexagonal symmetry, with anisotropy:

    ani = (V∥ - V⊥) / V · 100%

    where V∥ and V⊥ are the seismic velocities parallel and perpendicular to the
    symmetry axis and V is the average seismic velocity. P- and S-wave
    anisotropy are equal and pure elliptical (Levin and Park, 1997, GJI). The
    anisotropy axis is parallel to x1.

    Parameters
    ----------
    vp : float
        Isotropic P wave velocity. Units of km/s.
    vs : float
        Isotropic S wave velocity. Units of km/s.
    ani : float
        anisotropy% = (V∥ - V⊥) / V · 100%
    rho : float
        Bulk density of material. Units of g/cm^3

    Returns
    -------
    C : 6x6 array of floats
        Isotropic elastic stiffness tensor in Voigt matrix form.

    """

    C = np.zeros((6, 6))

    da = vp * ani / 100.0
    db = vs * ani / 100.0
    AA = rho * (vp - da / 2.0) ** 2.0
    CC = rho * (vp + da / 2.0) ** 2.0
    LL = rho * (vs + db / 2.0) ** 2.0
    NN = rho * (vs - db / 2.0) ** 2.0
    AC = rho * vp**2
    FF = -LL + np.sqrt(
        (2.0 * AC) ** 2 - 2.0 * AC * (AA + CC + 2.0 * LL) + (AA + LL) * (CC + LL)
    )

    C[2, 2] = AA
    C[1, 1] = AA
    C[0, 0] = CC

    C[2, 1] = AA - 2 * NN
    C[1, 2] = AA - 2 * NN

    C[2, 0] = FF
    C[0, 2] = FF
    C[1, 0] = FF
    C[0, 1] = FF

    C[5, 5] = LL
    C[4, 4] = LL
    C[3, 3] = NN

    return Material(C, rho, material_id="hexagonal material")


def decompose_C(material, symmetry="all"):
    """
    Decomposes an elastic tensor after the formulation set out in Browaeys and
    Chevrot, 2004. They propose a decomposition of the elastic tensor by
    representing it as triclinic elastic vector, X, before transforming it via
    a cascade of projections into a sum of vectors belonging to the different
    symmetry classes.

    Parameters
    ----------
    material : `anisotropy.materials.Material` object
        Material to decompose into elastic tensors representing each symmetry
        class.
    symmetry : str, optional
        Specify which component to return - default is to return a dictionary
        containing all decompositions.

    Returns
    -------
    decomposed_elements : dict
        Dictionary containing constitutive symmetry components as key, value
        pairs.

    """

    C = material.C.copy()

    decomposed_elements = {
        "isotropic": None,
        "hexagonal": None,
        "tetragonal": None,
        "orthorhombic": None,
        "monoclinic": None
    }
    for symmetry_class in decomposed_elements.keys():
        X = _C_tensor2vector(C)
        M = _projectors(symmetry_class)
        X_sc = np.dot(M, X)
        C_sc = _C_vector2tensor(X_sc)
        decomposed_elements[symmetry_class] = C_sc
        C -= C_sc

    return decomposed_elements


def _projectors(symmetry_class):
    """
    Utility function for serving the required projection matrices for each
    symmetry class

    Parameters
    ----------
    symmetry_class : str
        Name of symmetry class required.

    Returns
    -------
    M : array-like of float, shape(21, 21)

    """

    rt2 = np.sqrt(2)
    M = np.zeros((21, 21))

    if symmetry_class == "isotropic":
        M[0:3, 0:3] = 3/15
        M[0:3, 3:6] = rt2/15
        M[0:3, 6:9] = 2/15

        M[3:6, 0:3] = rt2/15
        M[3:6, 3:6] = 4/15
        M[3:6, 6:9] = -rt2/15

        M[6:9, 0:3] = 2/15
        M[6:9, 3:6] = -rt2/15
        M[6:9, 6:9] = 1/5

    if symmetry_class == "hexagonal":
        M[0:2, 0:2] = 3/8
        M[0:2, 5] = M[5, 0:2] = 1/(4*rt2)
        M[0:2, 8] = M[8, 0:2] = 1/4
        M[2, 2] = 1.
        M[3:5, 3:5] = M[6:8, 6:8] = M[8, 8] = 1/2
        M[5, 5] = 3/4
        M[5, 8] = M[8, 5] = -1/(2*rt2)

    if symmetry_class == "tetragonal":
        M[2, 2] = M[5, 5] = M[8, 8] = 1.
        M[0:2, 0:2] = M[3:5, 3:5] = M[6:8, 6:8] = 1/2

    if symmetry_class == "orthorhombic":
        np.fill_diagonal(M, 1)
        M[9:, 9:] = 0

    if symmetry_class == "monoclinic":
        np.fill_diagonal(M, 1)
        M[:, 10:12] = M[:, 13:15] = M[:, 16:18] = M[:, 19:21] = 0

    return M

def _C_tensor2vector(C):
    """
    Convert an elastic tensor, C, to an elastic vector, X, as set out in
    Equation 2.2 of Browaeys and Chevrot, 2004.

    Parameters
    ----------
    C : array-like of float, shape(6, 6)
        Elastic tensor for the material, in GPa, to be converted.

    Returns
    -------
    X : array-like of float, shape(21)
        Elastic vector representation of the elastic tensor, C, in GPa.

    """

    rt2 = np.sqrt(2)
    X = np.zeros(21)
    X[0:3] = C[0, 0], C[1, 1], C[2, 2]
    X[3:6] = rt2*C[1, 2], rt2*C[0, 2], rt2*C[0, 1]
    X[6:9] = 2*C[3, 3], 2*C[4, 4], 2*C[5, 5]
    X[9:12] = 2*C[0, 3], 2*C[1, 4], 2*C[2, 5]
    X[12:15] = 2*C[2, 3], 2*C[0, 4], 2*C[1, 5]
    X[15:18] = 2*C[1, 3], 2*C[2, 4], 2*C[0, 5]
    X[18:21] = 2*rt2*C[4, 5], 2*rt2*C[3, 5], 2*rt2*C[3, 4]

    return X


def _C_vector2tensor(X):
    """
    Convert an elastic vector, X, to an elastic vector, C, as set out in
    Equation 2.2 of Browaeys and Chevrot, 2004.

    Parameters
    ----------
    X : array-like of float, shape(21)
        Elastic vector representation of the elastic tensor, C, in GPa.

    Returns
    -------
    C : array-like of float, shape(6, 6)
        Elastic tensor for the material, in GPa, to be converted.

    """

    rt2 = np.sqrt(2)
    rt22 = rt2*2
    C = np.array([
        [    X[0], X[5]/rt2, X[4]/rt2,     X[9]/2,    X[13]/2,    X[17]/2],
        [X[5]/rt2,     X[1], X[3]/rt2,    X[15]/2,    X[10]/2,    X[14]/2],
        [X[4]/rt2, X[3]/rt2,     X[2],    X[12]/2,    X[16]/2,    X[11]/2],
        [  X[9]/2,  X[15]/2,  X[12]/2,     X[6]/2, X[20]/rt22, X[19]/rt22],
        [ X[13]/2,  X[10]/2,  X[16]/2, X[20]/rt22,     X[7]/2, X[18]/rt22],
        [ X[17]/2,  X[14]/2,  X[11]/2, X[19]/rt22, X[18]/rt22,     X[8]/2]
    ])

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
