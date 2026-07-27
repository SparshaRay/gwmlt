"""
Physics Utilities
"""

import numpy as np
import lalsimulation
import lalinference.imrtgr.nrutils as nr

from gwmlt.config import config
from gwmlt.constants import G_SI, C_SI, MSUN_SI
from gwmlt.core.sources import JFrameSpins, LFrameSpins, BinarySystem


# L and J frame interconversion functions ---------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def convert_jframe_to_lframe(
        mass_1 : float, mass_2 : float,
        spins_jframe : JFrameSpins,
        f_ref   : float | None = None,
        phi_ref : float = 0.0
    ) -> LFrameSpins :

    """
    Convert spins from the J-frame to the L-frame.

    Parameters
    ----------
    mass_1 : float
        Mass of the first object (in solar masses).
    mass_2 : float
        Mass of the second object (in solar masses).
    spins_jframe : JFrameSpins
        Spins in the J-frame.
    f_ref : float, optional
        Reference frequency (in Hz) for spin values. 
    phi_ref : float, optional
        Orbital phase at reference frequency (in radians).
        Default is 0.0

    Returns
    -------
    LFrameSpins
        Spins in the L-frame.
    """

    if f_ref is None : f_ref = config.waveform.f_ref

    inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = (
        lalsimulation.SimInspiralTransformPrecessingNewInitialConditions(
            spins_jframe.theta_jn,
            spins_jframe.phi_jl,
            spins_jframe.tilt_1,
            spins_jframe.tilt_2,
            spins_jframe.phi_12,
            spins_jframe.a_1,
            spins_jframe.a_2,
            mass_1 * MSUN_SI,
            mass_2 * MSUN_SI,
            f_ref,
            phi_ref,
        )
    )
    return LFrameSpins(
        inclination = inclination,
        spin1x = spin1x,
        spin1y = spin1y,
        spin1z = spin1z,
        spin2x = spin2x,
        spin2y = spin2y,
        spin2z = spin2z,
    )


def convert_lframe_to_jframe(
        mass_1 : float, mass_2 : float,
        spins_lframe : LFrameSpins,
        f_ref   : float | None = None,
        phi_ref : float = 0.0
    ) -> JFrameSpins :

    """
    Convert spins from the L-frame to the J-frame.

    Parameters
    ----------
    mass_1 : float
        Mass of the first object (in solar masses).
    mass_2 : float
        Mass of the second object (in solar masses).
    spins_lframe : LFrameSpins
        Spins in the L-frame.
    f_ref : float, optional
        Reference frequency (in Hz) for spin values. 
    phi_ref : float, optional
        Orbital phase at reference frequency (in radians).
        Default is 0.0
    
    Returns
    -------
    JFrameSpins
        Spins in the J-frame.
    """

    if f_ref is None : f_ref = config.waveform.f_ref

    thetajn, phijl, s1pol, s2pol, s12_deltaphi, spin1_a, spin2_a = (
        lalsimulation.SimInspiralTransformPrecessingWvf2PE(
            spins_lframe.inclination,
            spins_lframe.spin1x,
            spins_lframe.spin1y,
            spins_lframe.spin1z,
            spins_lframe.spin2x,
            spins_lframe.spin2y,
            spins_lframe.spin2z,
            mass_1,
            mass_2,
            f_ref,
            phi_ref,
        )
    )
    return JFrameSpins(
        a_1    = spin1_a,
        a_2    = spin2_a,
        tilt_1 = s1pol,
        tilt_2 = s2pol,
        phi_12 = s12_deltaphi,
        phi_jl = phijl,
        theta_jn = thetajn,
    )


# Merger, ringdown, and remnant properties --------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# All of these are taken from https://git.ligo.org/anuj.mishra/gwmat
# Please verify these function before using them.

def get_remnant_properties(m1, m2, a1=0.0, a2=0.0, tilt1=0.0, tilt2=0.0, phi12=0.0) :

    """
    Estimates the mass and spin of the final remnant based on the binary's initial conditions. Uses NR fits. 

    Parameters
    ----------
    m1, m2 : float or ndarray
        Component masses (in solar masses).
    a1, a2 : float, optional
        Dimensionless spin magnitudes. Defaults to 0.
    tilt1, tilt2 : float, optional
        Spin tilt angles relative to orbital angular momentum (in radians). Defaults to 0.
    phi12 : float, optional
        Difference in azimuthal angle between the two spins (in radians). Defaults to 0.

    Returns
    -------
    dict
        Containing `"m_f"` (final mass in solar masses) and `"a_f"` (final dimensionless spin).
    """

    # Final mass fits do not use phi12, so we set it to zero
    phi12_zero = np.array([0.0]) if np.isscalar(phi12) else np.zeros_like(phi12)
    
    m_f = nr.bbh_average_fits_precessing(
        m1=m1, m2=m2, chi1=a1, chi2=a2, tilt1=tilt1, tilt2=tilt2, 
        phi12=phi12_zero, quantity="Mf", fits=["UIB2016", "HL2016"]
    )
    a_f = nr.bbh_average_fits_precessing(
        m1=m1, m2=m2, chi1=a1, chi2=a2, tilt1=tilt1, tilt2=tilt2, 
        phi12=phi12, quantity="af", fits=["UIB2016", "HL2016", "HBR2016"]
    )

    return {"m_f": float(m_f), "a_f": float(a_f)}


def isco_radius(mass, spin) :

    """
    Returns the equatorial Innermost Stable Circular Orbit (ISCO) radius
    for a Kerr black hole given its mass and dimensionless spin parameter.
    It depends on whether the orbit is prograde or retrograde.

    References
    ----------
    - Eq. 2.21 of Bardeen et al. (1972): https://ui.adsabs.harvard.edu/abs/1972ApJ...178..347B/abstract
    - Eq. 1 of Hanna et al. (2008): https://arxiv.org/pdf/0801.4297.pdf

    Parameters
    ----------
    mass : float
        Mass of the black hole (the remnant) in solar masses.
    spin : float
        Dimensionless spin parameter of the black hole.

    Returns
    -------
    dict
        A dictionary containing `"prograde"` and `"retrograde"` ISCO radii in solar masses.
    """

    z1 = 1 + np.cbrt(1 - spin**2) * (np.cbrt(1 + spin) + np.cbrt(1 - spin))
    z2 = np.sqrt(3 * spin**2 + z1**2)
    
    comn = np.sqrt((3 - z1) * (3 + z1 + 2 * z2))
    r_prograde   = mass * (3 + z2 - comn)
    r_retrograde = mass * (3 + z2 + comn)
    
    return {"prograde": r_prograde, "retrograde": r_retrograde}


def kerr_isco_freq(final_mass, final_spin) :

    """
    Returns the GW frequency at the ISCO for given remnant mass and spin.
    Set spin to zero to get the Schwarzschild ISCO frequency.
    
    References
    ----------
    - Eq. 4 of Hanna et al. (2008): https://arxiv.org/pdf/0801.4297.pdf

    Parameters
    ----------
    final_mass : float
        Mass of the remnant black hole (in solar masses).
    final_spin : float
        Dimensionless spin magnitude of the remnant black hole.

    Returns
    -------
    dict
        GW frequencies at ISCO for `"prograde"` and `"retrograde"` cases (in Hz).
    """

    fac = 1.0 / (2 * np.pi * final_mass * (G_SI * MSUN_SI / C_SI**3))
    radii = isco_radius(final_mass, final_spin)
    
    # f_gw = 2 * f_orb
    f_prograde   = 2 * fac * (final_spin + (radii["prograde"]   / final_mass) ** 1.5) ** -1
    f_retrograde = 2 * fac * (final_spin + (radii["retrograde"] / final_mass) ** 1.5) ** -1
    
    return {"prograde": f_prograde, "retrograde": f_retrograde}


def kerr_isco_freq_from_bbh_pars(m1, m2, a1=0.0, a2=0.0, tilt1=0.0, tilt2=0.0, phi12=0.0) :
    """
    Wrapper to get Kerr ISCO frequency directly from the initial binary parameters.
    """
    remnant = get_remnant_properties(m1, m2, a1, a2, tilt1, tilt2, phi12)
    return kerr_isco_freq(remnant["m_f"], remnant["a_f"])


def bkl_isco_freq(m1, m2) :

    """
    Mass ratio dependent GW frequency at ISCO 
    derived from an estimate of the final spin of the remnant.
    
    References
    ----------
    - Buonanno, Kidder, Lehner (2007): https://arxiv.org/abs/0709.3839
    - Eq. 5 of Hanna et al. (2008): https://arxiv.org/abs/0801.4297v2

    Parameters
    ----------
    m1, m2 : float or ndarray
        Component masses of the binary (in solar masses).

    Returns
    -------
    float or ndarray
        Estimated GW frequency at ISCO (in Hz).
    """

    q = np.minimum(m1 / m2, m2 / m1)
    f_schwarz = kerr_isco_freq(m1 + m2, 0.0)["prograde"]
    return 0.5 * f_schwarz * (1 + 2.8 * q - 2.6 * q**2 + 0.8 * q**3)


def meco_freq(total_mass, q, spin1z, spin2z) :

    """
    Returns the Minimum Energy Circular Orbit (MECO) GW frequency.
    Uses lalsimulation/lib/LALSimIMRPhenomXUtilities.c#L44 phenomenological fit 
    to hybrid minimum energy circular orbit (MECO) function.
    Uses 3.5PN hybridised with test-particle limit.
    
    References
    ----------
    - Cabero et al., PRD, 95, 064016 (2017): https://arxiv.org/abs/1602.03134

    Parameters
    ----------
    total_mass : float or ndarray
        Total mass of the binary (in solar masses).
    q : float or ndarray
        Mass ratio (m_secondary / m_primary), where q <= 1.
    spin1z, spin2z : float or ndarray
        Z-components of dimensionless spin for primary and secondary.

    Returns
    -------
    float or ndarray
        MECO frequency (in Hz).
    """

    eta = np.minimum(q / (1.0 + q)**2, 0.25)
    v_meco = np.vectorize(lalsimulation.SimIMRPhenomXfMECO)
    f_meco_dimless = v_meco(eta, spin1z, spin2z)
    return BinarySystem.from_mtot(mtot=total_mass, q=q).to_physical_frequency(f_meco_dimless)


def ringdown_freq(m1, m2, a1=0.0, a2=0.0, tilt1=0.0, tilt2=0.0, phi12=0.0) :

    """
    Get the fundamental ringdown frequency. Uses the value for the omega_220
    QNM frequency from table VIII. of Berti, Cardoso and Will.
        
    References
    ----------
    - Berti, Cardoso, Will (2005): https://arxiv.org/pdf/gr-qc/0512160.pdf

    Parameters
    ----------
    m1, m2 : float or ndarray
        Component masses (in solar masses).
    a1, a2 : float, optional
        Dimensionless spin magnitudes. Defaults to 0.
    tilt1, tilt2 : float, optional
        Spin tilt angles relative to orbital angular momentum (rad). Defaults to 0.
    phi12 : float, optional
        Difference in azimuthal angle between the two spins (rad). Defaults to 0.

    Returns
    -------
    float
        Ringdown frequency (in Hz).
    """

    remnant = get_remnant_properties(m1, m2, a1, a2, tilt1, tilt2, phi12)
    m_f, a_f = remnant["m_f"], remnant["a_f"]
    
    # Dimensionless ringdown frequency
    f_dimless = 1.5251 - 1.1568 * (1 - a_f) ** 0.1292

    hz_conversion_factor = 1.0 / (2 * np.pi * m_f * (G_SI * MSUN_SI / C_SI**3))
    return hz_conversion_factor * f_dimless


def recommended_reference_freq(total_mass, q, spin1z, spin2z, f_ref=20.0, f_start=10.0, fudge_factor=0.97) :

    """
    Recommends an appropriate reference frequency for data analysis.

    References
    ----------
    Based on : asimov/pipelines/pe-configurator/peconfigurator/proc_samples.py#432
    
    Parameters
    ----------
    total_mass : float or ndarray
        Total mass of the binary (in solar masses).
    q : float or ndarray
        Mass ratio (m_secondary / m_primary), where q <= 1.
    spin1z, spin2z : float or ndarray
        Z-components of dimensionless spin for primary and secondary.
    f_ref : float, optional
        The default desired reference frequency in Hz (default 20.0).
    f_start : float, optional
        The starting frequency of the analysis (Hz) (default 10.0).
    fudge_factor : float, optional
        Safety margin applied to the MECO frequency (default 0.97).

    Returns
    -------
    float
        A recommended analysis reference frequency (in Hz).

    Raises
    ------
    ValueError
        If the required starting frequency is higher than the recommended reference frequency.
    """

    f_meco = meco_freq(total_mass, q, spin1z, spin2z)
    scaled_meco = fudge_factor * f_meco
    recommended_ref = np.minimum(np.floor(scaled_meco), f_ref)

    if np.any(recommended_ref <= f_start):
        raise ValueError(
            f"Provided f_start ({f_start} Hz) is higher than or equal to the recommended f_ref ({recommended_ref} Hz).")
        
    return recommended_ref