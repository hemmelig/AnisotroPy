# -*- coding: utf-8 -*-
"""
Modelling of the effective seismic anisotropy of a system of N arbitrarily
anisotropic layers, based either on the analytical system of equations outlined
by Silver and Savage (1994) or by directly applying a sequence of splitting
operators to the first-derivative of a Gaussian wavelet.

For more information, please read:

    Silver, P.G. and Savage, M.K., 1994. The interpretation of shear-wave
    splitting parameters in the presence of two anisotropic layers. Geophysical
    Journal International, 119(3), pp.949-963.

If you use this code, please cite this article.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import copy

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from anisotropy.materials import Material, voigt_reuss_hill_average


π = np.pi


class ElasticLayer:
    """
    Small class to encapsulate information that defines a layer composed of an
    (an)isotropic material and provide a simple interface to a suite of useful
    functions to query its properties.

    Attributes
    ----------
    material : `anisotropy.Material` object
        (An)isotropic material, whose elastic properties are fully described by
        its stiffness tensor and density.
    thickness : float
        Thickness of the elastic layer, in km.
    ψ_a : float
        Azimuthal alignment of crystalline lattice a-axis of the material.
    dip : float, optional
        Angle of dip of the layer, relative to horizontal, in degrees.
        Defaults to 0.
    fractional_alignment : float
        Degree of mixing between completely aligned material and its isotropic
        component. Must be between 0 and 1.

    Methods
    -------
    traversal_distance
        Calculate the distance travelled within a dipping layer.
    effective_dt
        Calculate the effective delay time for a ray passing through the layer.
    effective_fast
        Calculate the effective fast orientation direction for a ray passing
        through the layer.

    """

    def __init__(self, material, thickness, ψ_a, dip=0., fractional_alignment=0.3):
        """Instantiate the Layer object."""

        self.thickness = thickness
        self.dip = dip
        self.ψ_a = ψ_a
        self.fractional_alignment = fractional_alignment

        self.material = copy.deepcopy(material)

    def traversal_distance(self, angle_of_incidence, azimuth=np.arange(0, 361, 1)):
        """
        Calculates effective distance that a ray traverses when passing through
        the layer at some given angle of incidence

        Parameters
        ----------
        angle_of_incidence : float
            Angle of incidence of a ray with the surface, relative to vertical,
            in degrees.
        azimuth : float or list of floats, optional
            Azimuth(s) of ray(s) at which to calculate traversal distances.
            Defaults to the full range of azimuths at 1 degree increments.

        Returns
        -------
        distance : float
            Effective distance travelled within the layer, in km.

        """

        dip, azimuth = np.deg2rad(self.dip), np.deg2rad(azimuth)

        apparent_dip = np.rad2deg(np.arctan(-np.tan(dip)*np.cos(azimuth)))

        γ = np.deg2rad(90 - apparent_dip - angle_of_incidence)
        β = np.deg2rad(90 + apparent_dip)

        return self.thickness * np.sin(β) / np.sin(γ)

    def effective_dt(self, angle_of_incidence, azimuth=np.arange(0, 361, 1)):
        """
        Calculates the effective delay times of a ray traversing the layer with
        a given angle of incidence with the surface, as measured from the
        vertical.

        Parameters
        ----------
        angle_of_incidence : float
            Angle of incidence of a ray with the surface, relative to vertical,
            in degrees.
        azimuth : float or list of floats, optional
            Azimuth(s) of ray(s) at which to calculate traversal distances.
            Defaults to the full range of azimuths at 1 degree increments.

        Returns
        -------
        effective_dt : `numpy.ndarray` of float
            Effective delay times.

        """

        inclination = 90 - angle_of_incidence
        distance_in_layer = self.traversal_distance(angle_of_incidence, azimuth)
        phase_velocities = self.material.phase_velocities(inclination, azimuth)

        return distance_in_layer*(1/np.array(phase_velocities[2])
                                  - 1/np.array(phase_velocities[1]))

    def effective_fast(self, angle_of_incidence, azimuth=np.arange(0, 361, 1)):
        """
        Calculates the effective orientation of polarisation of the fast shear
        wave traversing the layer with a given angle of incidence with the
        surface, as measured from the vertical.

        Parameters
        ----------
        angle_of_incidence : float
            Angle of incidence of a ray with the surface, relative to vertical,
            in degrees.
        azimuth : float or list of floats, optional
            Azimuth(s) of ray(s) at which to calculate traversal distances.
            Defaults to the full range of azimuths at 1 degree increments.

        Returns
        -------
        effective_fast : `numpy.ndarray` of float
            Effective orientation of polarisation of the fast shear wave.

        """

        inclination = 90 - angle_of_incidence
        phase_velocities = self.material.phase_velocities(inclination, azimuth)

        return _unwind_angle(azimuth + phase_velocities[3])

    @property
    def material(self):
        """Get and set the constitutive material of the layer."""

        return self._material

    @material.setter
    def material(self, value):
        iso_olivine = Material(value.isotropic, value.rho)
        value.C = value.rotate(0, -self.dip, self.ψ_a, order=[3, 2, 1])

        C_vrh, rho_vrh, *_ = voigt_reuss_hill_average(
            [iso_olivine, value],
            [1-self.fractional_alignment, self.fractional_alignment])

        self._material = Material(C_vrh, rho_vrh)


def fit_2layer_model(observations, material, angle_of_incidence):
    """
    Perform a gridsearch over the 4 free parameters that define a 2 layer
    model, namely (T_upper, φ_upper) and (T_lower, φ_lower)

    TODO
      - Write tests
      - Testing the free parameters in increments of 5/10 km (thicknesses)
        and 2.5/5 degrees for azimuth? Already doing some additional model
        sweeping in the _calculate_misfit function.
      - Allow the code to smartly account for the phase? e.g. calculate the
        angle of incidence for different phases at different distances - we
        are searching over how the phase interacts with the model, so makes
        sense to include this information.

    Parameters
    ----------
    observations : pandas.DataFrame object
        This should either be a DataFrame or perhaps a custom data format
        develop as part of the overall package. This would mean that results
        generated with other splitting codes could be used by simply writing
        a parser that converts the outputs of the external codes into the
        custom data file format.
    material : `anisotropy.Material` object
        (An)isotropic material, whose elastic properties are fully described by
        its stiffness tensor and density.
    angle_of_incidence : float
        Angle of incidence of a ray with the surface, relative to vertical, in
        degrees.

    """

    T1s, T2s = np.arange(50, 150, 25), np.arange(50, 150, 25)
    a1s, a2s = np.arange(-90, 95, 10), np.arange(-90, 95, 10)
    misfits = []
    print("Performing misfit calculation...")
    for T1 in T1s:
        for a1 in a1s:
            print(f"   ...Layer 1 - thickness = {T1} km, ψ_a = {a1}...")
            layer1 = ElasticLayer(material, thickness=T1, ψ_a=a1, fractional_alignment=0.3888)
            for T2 in T2s:
                for a2 in a2s:
                    print(f"      ...Layer 2 - thickness = {T2} km, ψ_a = {a2}...")
                    layer2 = ElasticLayer(material, thickness=T2, ψ_a=a2, fractional_alignment=0.3888)

                    layers = [layer2, layer1]
                    misfits.append(_calculate_misfit(observations, angle_of_incidence, layers))

    return misfits


def _calculate_misfit(observations, layers, angle_of_incidence):
    """
    Calculate the misfit between measured values of (φ, δt) and a trial N layer
    model. This is based on the equations set out in:

        Liddell, M.V., Bastow, I., Darbyshire, F., Gilligan, A. and Pugh, S.,
        2017. The formation of Laurentia: Evidence from shear wave splitting.
        Earth and Planetary Science Letters, 479, pp.170-178.

        Merry, T.A., Bastow, I.D., Kounoudis, R., Ogden, C.S., Bell, R.E. and
        Jones, L., 2021. The influence of the North Anatolian Fault and a
        fragmenting slab architecture on upper mantle seismic anisotropy in the
        eastern Mediterranean. Geochemistry, Geophysics, Geosystems, 22(9),
        p.e2021GC009896.

    Parameters
    ----------
    observations : pandas.DataFrame object
        This should either be a DataFrame or perhaps a custom data format
        develop as part of the overall package. This would mean that results
        generated with other splitting codes could be used by simply writing
        a parser that converts the outputs of the external codes into the
        custom data file format.
    layers : list of `anisotropy.ElasticLayer` object
        Layers in the model.
    angle_of_incidence : float
        Angle of incidence of a ray with the surface, relative to vertical, in
        degrees.

    Returns
    -------
    M : float
        Calculated misfit for the set of measurements and the given model.

    """

    φφ, δδt, dφφ, dδδt, backazimuths = observations
    Γφ, Γδt = [], []
    for φ, δt, dφ, dδt, backazimuth in zip(φφ, δδt, dφφ, dδδt, backazimuths):
        azimuths = np.arange(backazimuth-5, backazimuth+6, 1)
        effective_results = model(0.125, layers, angle_of_incidence, azimuths)
        effective_φφ, effective_δδt = effective_results

        Γφx = np.array([abs(φ - effective_φ) - dφ
                        for effective_φ in effective_φφ])
        Γφx[Γφx < 0] = 0
        Γφ.append(min(Γφx))

        Γδtx = np.array([abs(δt - effective_δt) - dδt
                         for effective_δt in effective_δδt])
        Γδtx[Γδtx < 0] = 0
        Γδt.append(min(Γδtx))

    rms_φ = np.sqrt(np.mean([Γφ_i**2 for Γφ_i in Γφ]))
    rms_δt = np.sqrt(np.mean([Γδt_i**2 for Γδt_i in Γδt]))
    print(rms_φ, rms_δt)

    M = rms_φ/np.std(φφ) + rms_δt/np.std(δδt)

    return M


def model(frequency, angle_of_incidence, layers, azimuth=np.arange(0, 361, 1),
          method="silver_savage"):
    """
    Calculate the effective splitting of a ray passing through N arbitrarily
    anisotropic layers.

    Parameters
    ----------
    frequency : int or float
        Frequency of effective wavelet for purposes of modelling.
    angle_of_incidence : float
        Angle of incidence of a ray with the surface, relative to vertical, in
        degrees.
    layers : list of `anisotropy.ElasticLayer` object
        Layers in the model.
    azimuth : float or list of floats, optional
        Azimuth(s) of ray(s) at which to calculate traversal distances. Defaults
        to the full range of azimuths at 1 degree increments.
    method : str, {"silver_savage", "gaussian_wavelet"}
        Choose the method of modelling - either analytical or forward modelling.

    Returns
    -------
    effective_results : list of list of float
        Effective fast direction and delay times for the N-layer system.

    """

    if method == "silver_savage":
        return _silver_savage_94(frequency, angle_of_incidence, layers, azimuth)


def _silver_savage_94(frequency, angle_of_incidence, layers, azimuth=np.arange(0, 361, 1)):
    """
    Calculates the effective splitting for N layers using the method set out in
    Silver and Savage (1994):

        Silver, P.G. and Savage, M.K., 1994. The interpretation of shear-wave
        splitting parameters in the presence of two anisotropic layers.
        Geophysical Journal International, 119(3), pp.949-963.

    Parameters
    ----------
    frequency : int or float
        Frequency of effective wavelet for purposes of modelling.
    angle_of_incidence : float
        Angle of incidence of a ray with the surface, relative to vertical, in
        degrees.
    layers : list of `anisotropy.ElasticLayer` object
        Layers in the model.
    azimuth : float or list of floats, optional
        Azimuth(s) of ray(s) at which to calculate traversal distances. Defaults
        to the full range of azimuths at 1 degree increments.

    Returns
    -------
    effective_φ : float
        Effective φ value for the N-layer system, in degrees.
    effective_δt : float
        Effective δt value for the N-layer system, in seconds.

    """

    # Check for a single layer
    if not isinstance(layers, list):
        print("Must provide a list of `anisotropy.ElasticLayer`.")
        raise ValueError

    effective_φφ, effective_δδt = [], []
    for az in azimuth:
        φφ = [layer.effective_fast(angle_of_incidence, az) for layer in layers]
        δδt = [layer.effective_dt(angle_of_incidence, az) for layer in layers]
        φφ, δδt = _aggregate_layers(φφ, δδt)
            
        # Drop any layers with δt=0
        zero_δt = np.where(δδt == 0)
        φφ, δδt = np.delete(φφ, zero_δt), np.delete(δδt, zero_δt)

        # Unwind fast directions
        φφ = _unwind_angle(φφ)

        # Process
        ω = 2.0 * π * frequency
        θθ = (ω / 2) * δδt
        αα = (2*π/180) * (φφ - az)

        S = np.prod(np.cos(θθ))
        Cc = S * np.sum(np.tan(θθ) * np.cos(αα))
        Cs = S * np.sum(np.tan(θθ) * np.sin(αα))

        ap, ap_perp = 0, 0
        for i in range(len(φφ)-1):
            for j in range(i+1, len(φφ)):
                ap += np.tan(θθ[i])*np.tan(θθ[j])*np.cos(αα[i] - αα[j])
                ap_perp += np.tan(θθ[i])*np.tan(θθ[j])*np.sin(αα[i] - αα[j])
        ap = S * (1 - ap)
        ap_perp *= S
        αa = np.arctan((ap_perp**2 + Cs**2) / (ap_perp*ap + Cs*Cc))
        θa = np.arctan(ap_perp / (Cs*np.cos(αa) - Cc*np.sin(αa)))

        effective_φ = _unwind_angle(az + (αa * 90 / π))
        effective_δt = (2 / ω) * θa

        # Handle zero δts
        if effective_δt < 0:
            effective_φ = _unwind_angle(effective_φ + 90)
            effective_δt = abs(effective_δt)

        effective_φφ.append(effective_φ)
        effective_δδt.append(effective_δt)

    return effective_φφ, effective_δδt


def _plot_effective_splitting(effective_φφ, effective_δδt, layers, angle_of_incidence, observations=None):
    """
    Plots the results of an N-layer model as a function of back-azimuth (source
    polarisation).

    Parameters
    ----------
    effective_φφ : numpy.ndarray of float
        Effective φ values for the N-layer system, in degrees.
    effective_δδt : numpy.ndarray of float
        Effective δt values for the N-layer system, in seconds.
    layers : list of `anisotropy.ElasticLayer` object
        Layers in the model.
    angle_of_incidence : float
        Angle of incidence of a ray with the surface, relative to vertical, in
        degrees.
    observations : pandas.DataFrame object, optional
        This should either be a DataFrame or perhaps a custom data format
        develop as part of the overall package. This would mean that results
        generated with other splitting codes could be used by simply writing
        a parser that converts the outputs of the external codes into the
        custom data file format.

    """

    fig, axes = plt.subplots(2, figsize=(16, 12), constrained_layout=True)
    source_polarisations = np.arange(0, 361, 1)

    ax = axes[0]
    effective_φφ_sequences = _wrap_φφ(effective_φφ)
    for effective_φφ_sequence in effective_φφ_sequences:
        ax.plot(source_polarisations, effective_φφ_sequence, color="k", lw=1)
    for layer, c, label in zip(layers, ["#99d8c9", "#f768a1"], ["Upper", "Lower"]):
        fasts = layer.effective_fast(angle_of_incidence, source_polarisations)
        ax.plot(source_polarisations, fasts, color=c, lw=1, linestyle="--",
                label=f"{label} layer")

    # Beautify
    ax.set_xlim([0, 360])
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_xticklabels([])
    ax.set_ylim([-90, 90])
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ytickrange = np.arange(-90, 91, 30)
    ax.set_yticks(ytickrange)
    ax.set_yticklabels(ytickrange)
    ax.set_ylabel("Effective φ, °")

    ax = axes[1]
    ax.plot(source_polarisations, effective_δδt, color="k", lw=1,
            label="Effective model")
    for layer, c, label in zip(layers, ["#99d8c9", "#f768a1"], ["Upper", "Lower"]):
        δt = layer.effective_dt(angle_of_incidence, source_polarisations)
        ax.plot(source_polarisations, δt, color=c, lw=1, linestyle="--",
                label=f"{label} layer")
    ax.legend(fontsize=14)

    # Beautify
    ax.set_xlim([0, 360])
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    xtickrange = np.arange(0, 361, 60)
    ax.set_xticks(xtickrange)
    ax.set_xticklabels(xtickrange)
    ax.set_ylim([0, 4])
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax.set_xlabel("Source polarisation, °")
    ax.set_ylabel("Effective delay time, s")

    # Add any observations
    if observations is not None:
        axes[0].scatter(observations["Baz."], observations["PHI"], marker="s",
                        s=25, color="k")
        axes[1].scatter(observations["Baz."], observations["TLAG"], marker="d",
                        s=25, color="k")


def _wrap_φφ(effective_φφ):
    """
    Seeks φ value "wraparounds" and produces two sequences that, when plotted,
    appear correct.

    Parameters
    ----------
    effective_φφ : numpy.ndarray of float
        Effective φ values for the N-layer system, in degrees.

    Returns
    -------
    sequences : list of lists of float
        Wrapped φ sequences ready for plotting.
    
    """

    effective_φφ = np.concatenate([effective_φφ[:-1], effective_φφ])

    # Seeking points where data "wraps around" and needs to be split
    correction, corrections = 0, []
    for i, φ in enumerate(effective_φφ[1:]):
        if (φ - effective_φφ[i]) < -90:
            correction = 180
        elif (φ - effective_φφ[i]) > 90:
            correction = 0
        corrections.append(correction)

    effective_φφ, corrections = effective_φφ[360:], corrections[359:]
    sequence1 = [sum(x) for x in zip(effective_φφ, corrections)]
    sequence2 = [x - 180 for x in sequence1]

    return [sequence1, sequence2]


def _aggregate_layers(φφ, δδt, threshold=1):
    """
    Agglomerate splitting operators which are (near) parallel or perpendicular.
    
    Parameters
    ----------
    φφ : list of floats
        Fast orientation values for each layer.
    δt : list of floats
        Delay time values for each layer.
    threshold : float, optional
        Number of degrees within which to seek layers to aggregate. Default 1.

    Returns
    -------
    φφ_a : numpy.ndarray of floats
        Fast orientation values for aggregated layers.
    δδt_a : numpy.ndarray of floats
        Delay time values for aggregated layers.

    """

    aggregate_φ, aggregate_δt = np.empty(len(φφ)), np.empty(len(φφ))
    aggregate_φ[:], aggregate_δt[:] = np.nan, np.nan
    aggregate_φ[0], aggregate_δt[0] = φφ[0], δδt[0]

    # Process φs. If the reference direction and the comparison direction are
    # (nearly) identical, sum the δts. If they differ by (nearly) 90 degrees,
    # take the difference.
    i_ref, i_agg, i = 0, 0, 1
    while i < len(φφ):
        if abs(φφ[i_ref] - φφ[i]) < threshold:
            # Add the δt
            aggregate_δt[i_agg] += δδt[i]
        elif abs(abs(φφ[i_ref] - φφ[i]) - 90) < threshold:
            # Subtract the δt
            aggregate_δt[i_agg] -= δδt[i]
        else:
            # Proceed to next element
            i_ref = i
            i_agg += 1
            aggregate_φ[i_agg] = φφ[i_ref]
            aggregate_δt[i_agg] = δδt[i_ref]
        i += 1

    nan_positions = ~np.isnan(aggregate_δt)
    φφ_a, δδt_a = aggregate_φ[nan_positions], aggregate_δt[nan_positions]

    # Check for any negative numbers
    for i in np.where(δδt_a < 0)[0]:
        φφ_a[i] += 90
        δδt_a[i] = abs(δδt_a[i])

    return φφ_a, δδt_a


def _unwind_angle(angle):
    """
    Utility function to find the equivalent angle between -90 and 90 degrees
    assuming 180 degree periodicity.

    Parameters
    ----------
    angle : float or list of floats
        Angles to be tested and brought within the range of -90 to 90 degrees.

    Returns
    -------
    angle : float or list of floats
        Angles within the range of -90 to 90 degrees.

    """

    angle = np.copy(np.array([angle]))

    # Shift p/m 180 to 0
    angle -= 180*np.fix(angle/180)

    # Bring within range -90 < angle < 90
    angle[np.where(angle <= -90)] += 180
    angle[np.where(angle > 90)] -= 180

    return angle[0]
