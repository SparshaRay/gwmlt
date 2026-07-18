"""
Physics Utilities
"""

import numpy as np
import lalsimulation

from ..config import config
from ..constants import G_SI, C_SI, MSUN_SI
from ..core.sources import JFrameSpins, LFrameSpins


# L and J frame interconversion functions ---------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def convert_jframe_to_lframe(
        mass_1 : float, mass_2 : float,
        spins_jframe : JFrameSpins,
        f_ref   : float = config.waveform.f_ref,
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
        inclination=inclination,
        spin1x=spin1x,
        spin1y=spin1y,
        spin1z=spin1z,
        spin2x=spin2x,
        spin2y=spin2y,
        spin2z=spin2z,
    )


def convert_lframe_to_jframe(
        mass_1 : float, mass_2 : float,
        spins_lframe : LFrameSpins,
        f_ref   : float = config.waveform.f_ref,
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


# Dimensionless and physical interconversion classes ----------------------------------------------
# -------------------------------------------------------------------------------------------------

class DimfulToDimless :

    """
    Class to convert physical quantities to dimensionless quantities.
    """

    def __init__(self, mass_1: float, mass_2: float) -> None :

        """
        Initialize the DimfulToDimless class.

        Parameters
        ----------
        mass_1 : float
            Mass of the first object (in solar masses).
        mass_2 : float
            Mass of the second object (in solar masses).
        """

        self.mass_1 = mass_1
        self.mass_2 = mass_2
        self.total_mass = mass_1 + mass_2

        self.T_scale = (G_SI * self.total_mass * MSUN_SI) / (C_SI**3)
        self.L_scale = (G_SI * self.total_mass * MSUN_SI) / (C_SI**2)

    def get_mass_ratio(self) -> float :
        """
        Get the mass ratio.
        """
        return self.mass_1 / self.mass_2

    def get_dimless_masses(self) -> tuple[float, float] :
        """
        Get the dimensionless masses.
        """
        return self.mass_1/self.total_mass, self.mass_2/self.total_mass
    
    def get_dimless_time(self, time: float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionful time to dimensionless time.
        """
        return time / self.T_scale
    
    def get_dimless_length(self, length: float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionful length to dimensionless length.
        """
        return length / self.L_scale
    
    def get_dimless_frequency(self, frequency: float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionful frequency to dimensionless frequency.
        """
        return frequency * self.T_scale
    
    def get_dimless_Sz_from_chiz(self, chi1z: float, chi2z: float) -> tuple[float, float] :
        """
        Convert dimensionless z spin to dimensionless z angular momenta.
        """
        S1z_bar = chi1z * (self.mass_1 / self.total_mass)**2
        S2z_bar = chi2z * (self.mass_2 / self.total_mass)**2
        return S1z_bar, S2z_bar
    

class DimlessToDimful :

    """
    Class to convert dimensionless quantities to physical quantities.
    """

    def __init__(self, total_mass : float) -> None :

        """
        Initialize the DimlessToDimful class.

        Parameters
        ----------
        total_mass : float
            Total mass of the system (in solar masses).
        """

        self.total_mass = total_mass
        self.T_scale = (G_SI * self.total_mass * MSUN_SI) / (C_SI**3)
        self.L_scale = (G_SI * self.total_mass * MSUN_SI) / (C_SI**2)

    def get_dimful_masses(self, mass_ratio : float) -> tuple[float, float] :
        """
        Get the dimensionful masses from the mass ratio.
        """
        m1 = self.total_mass * mass_ratio / (1 + mass_ratio)
        m2 = self.total_mass / (1 + mass_ratio)
        return m1, m2
    
    def get_dimful_time(self, tau : float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionless time to dimensionful time.
        """
        return tau * self.T_scale
    
    def get_dimful_length(self, length_bar : float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionless length to dimensionful length.
        """
        return length_bar * self.L_scale
    
    def get_dimful_frequency(self, frequency_bar : float | np.ndarray) -> float | np.ndarray :
        """
        Convert dimensionless frequency to dimensionful frequency.
        """
        return frequency_bar / self.T_scale
    
    def get_chiz_from_dimless_Sz(self, S1z_bar : float, S2z_bar : float, mass_ratio : float) -> tuple[float, float] :
        """
        Convert dimensionless z angular momenta to dimensionless z spin.
        """
        m1, m2 = self.get_dimful_masses(mass_ratio)
        chi1z = S1z_bar / (m1 / self.total_mass)**2
        chi2z = S2z_bar / (m2 / self.total_mass)**2
        return chi1z, chi2z