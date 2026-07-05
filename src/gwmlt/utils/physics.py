"""
Physics Utilities
"""

import lalsimulation

from ..config import config
from ..constants import MSUN_SI
from ..structures.sources import BBHSystem


# L and J frame interconversion functions ---------------------------------------------------------

def _convert_jframe_to_lframe(
        mass_1 : float, mass_2 : float,
        spins_jframe : BBHSystem.JFrameSpins,
        f_ref   : float = config.waveform.f_ref,
        phi_ref : float = 0.0
    ) -> BBHSystem.LFrameSpins :

    """
    Convert spins from the J-frame to the L-frame.

    Parameters
    ----------
    mass_1 : float
        Mass of the first object (in solar masses).
    mass_2 : float
        Mass of the second object (in solar masses).
    spins_jframe : BBHSystem.JFrameSpins
        Spins in the J-frame.
    f_ref : float, optional
        Reference frequency (in Hz) for spin values. 
    phi_ref : float, optional
        Orbital phase at reference frequency (in radians).
        Default is 0.0

    Returns
    -------
    BBHSystem.LFrameSpins
        Spins in the L-frame.
    """

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
    
    return BBHSystem.LFrameSpins(
        inclination=inclination,
        spin1x=spin1x,
        spin1y=spin1y,
        spin1z=spin1z,
        spin2x=spin2x,
        spin2y=spin2y,
        spin2z=spin2z,
    )


def _convert_lframe_to_jframe(
        mass_1 : float, mass_2 : float,
        spins_lframe : BBHSystem.LFrameSpins,
        f_ref   : float = config.waveform.f_ref,
        phi_ref : float = 0.0
    ) -> BBHSystem.JFrameSpins :

    """
    Convert spins from the L-frame to the J-frame.

    Parameters
    ----------
    mass_1 : float
        Mass of the first object (in solar masses).
    mass_2 : float
        Mass of the second object (in solar masses).
    spins_lframe : BBHSystem.LFrameSpins
        Spins in the L-frame.
    f_ref : float, optional
        Reference frequency (in Hz) for spin values. 
    phi_ref : float, optional
        Orbital phase at reference frequency (in radians).
        Default is 0.0
    
    Returns
    -------
    BBHSystem.JFrameSpins
        Spins in the J-frame.
    """

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
    return BBHSystem.JFrameSpins(
        a_1    = spin1_a,
        a_2    = spin2_a,
        tilt_1 = s1pol,
        tilt_2 = s2pol,
        phi_12 = s12_deltaphi,
        phi_jl = phijl,
        theta_jn = thetajn,
    )