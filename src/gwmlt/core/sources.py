"""
Structs for GW Source Parameters
"""

from dataclasses import dataclass, field, InitVar
from lal import LIGOTimeGPS
import numpy as np

from gwmlt.constants import G_SI, C_SI, MSUN_SI


# ![NOTE] - Bilby and LAL/PyCBC parameter naming conventions and their equivalence :
#
# | Bilby Key           | PyCBC Key        |
# |---------------------|------------------|
# | mass_1              | mass1            |
# | mass_2              | mass2            |
# | spin_1x             | spin1x           |
# | spin_1y             | spin1y           |
# | spin_1z             | spin1z           |
# | spin_2x             | spin2x           |
# | spin_2y             | spin2y           |
# | spin_2z             | spin2z           |
# | iota                | inclination      |
# | luminosity_distance | distance         |
# | phase               | coa_phase        |
# | psi                 | polarization     |
# | geocent_time        | trigger_time     |
# | approximant         | wf_approximant   |
# | f_start             | f_lower          |


# Base classes for CBC sources --------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class JFrameSpins :

    """
    Spins and orientation in the J-frame (used by Bilby).

    Parameters
    ----------
    a_1 : float, optional
        Dimensionless spin magnitude of the primary binary. 
        Must lie in the interval [0, 1). Default is 0.0.
    a_2 : float, optional
        Dimensionless spin magnitude of the secondary binary. 
        Must lie in the interval [0, 1). Default is 0.0.
    tilt_1 : float, optional
        Zenith angle between the primary binary's spin vector and the orbital angular 
        momentum vector in radians, defined at the reference frequency. 
        Default is 0.0.
    tilt_2 : float, optional
        Zenith angle between the secondary binary's spin vector and the orbital angular 
        momentum vector in radians, defined at the reference frequency. 
        Default is 0.0.
    phi_12 : float, optional
        Difference between the azimuthal angles of the two spin vectors in the 
        plane perpendicular to the orbital angular momentum in radians. 
        Default is 0.0.
    phi_jl : float, optional
        Azimuthal angle of the orbital angular momentum L on its cone of precession 
        around the total angular momentum J in radians.
        Default is 0.0.
    theta_jn : float, optional
        Angle between the total angular momentum J and the line of sight vector in radians. 
        Default is 0.0.
    """

    a_1      : float = 0.0
    a_2      : float = 0.0

    tilt_1   : float = 0.0
    tilt_2   : float = 0.0

    phi_12   : float = 0.0
    phi_jl   : float = 0.0

    theta_jn : float = 0.0


@dataclass(frozen=True)
class LFrameSpins :

    """
    Cartesian spin components and orientation in the L0-frame.
    This matches the reference frame coordinates used by PyCBC and LALSimulation.

    Parameters
    ----------
    inclination : float, optional
        Inclination angle defined as the angle between the orbital angular momentum L
        and the line of sight vector in radians. Default is 0.0.
    spin1x : float, optional
        The x-component of the primary binary's dimensionless spin vector. Default is 0.0.
    spin1y : float, optional
        The y-component of the primary binary's dimensionless spin vector. Default is 0.0.
    spin1z : float, optional
        The z-component of the primary binary's dimensionless spin vector (aligned spin).
        Default is 0.0.
    spin2x : float, optional
        The x-component of the secondary binary's dimensionless spin vector. Default is 0.0.
    spin2y : float, optional
        The y-component of the secondary binary's dimensionless spin vector. Default is 0.0.
    spin2z : float, optional
        The z-component of the secondary binary's dimensionless spin vector (aligned spin).
        Default is 0.0.
    """

    inclination : float = 0.0

    spin1x      : float = 0.0
    spin1y      : float = 0.0
    spin1z      : float = 0.0

    spin2x      : float = 0.0
    spin2y      : float = 0.0
    spin2z      : float = 0.0


@dataclass(frozen=True)
class BinarySystem :

    """
    Base class containing parameters of any Compact Binary Coalescence (CBC) source.

    This struct contains both intrinsic parameters (masses, spins, etc.) 
    and extrinsic parameters (sky location, orientation, time, etc.) of 
    a CBC source both in L-frame and J-frame coordinates.

    Parameters
    ----------
    mass_1 : float
        Detector-frame mass of the first binary in solar masses.
        This is referred to as the primary (heavier), but the ordering is not enforced.
    mass_2 : float
        Detector-frame mass of the second binary in solar masses.
        This is referred to as the secondary (lighter), but the ordering is not enforced.
    lambda_1 : float, optional
        Dimensionless tidal deformability of the primary binary.
        Default is 0.0 (i.e. no tidal effects).
    lambda_2 : float, optional
        Dimensionless tidal deformability of the secondary binary.
        Default is 0.0 (i.e. no tidal effects).
    ra : float, optional
        Right ascension of the source in radians.
        Default is 0.0.
    dec : float, optional
        Declination of the source in radians.
        Default is 0.0.
    phase : float, optional
        Coalescence orbital phase of the binary system in radians. 
        Default is 0.0.
    psi : float, optional
        Polarization angle of the gravitational wave in radians.
        Default is 0.0.
    geocent_time : float or LIGOTimeGPS, optional
        The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds. 
        Default is 0.0.
    luminosity_distance : float, optional
        Luminosity distance to the source in megaparsec (Mpc).
        Default is 100.0.
    spins : JFrameSpins, or LFrameSpins, or None, optional
        Dataclass containing spin parameters. Default is None, which implies non-spinning binaries.
    
    Attributes
    ----------
    j_frame : JFrameSpins
        Dataclass containing J-frame spin parameters.
    l_frame : LFrameSpins
        Dataclass containing L-frame spin parameters.
    mtot : float
        Total mass of the binary system.
    mchirp : float
        Chirp mass of the binary system.
    q : float
        Mass ratio (m_secondary / m_primary).
    eta : float
        Symmetric mass ratio (m1 * m2 / M_tot^2).
    chi_eff : float
        Effective aligned spin parameter.
    chi_p : float
        Effective precession spin parameter.
    t_scale : float
        Characteristic time scale (GM_tot/c^3) in seconds.
    l_scale : float
        Characteristic length scale (GM_tot/c^2) in meters.
    dimless_masses : tuple[float, float]
        Dimensionless masses (m1/M_tot, m2/M_tot).
    dimless_spins_z : tuple[float, float]
        Dimensionless z-angular momenta.

    Methods
    -------
    from_mtot(mtot, q=None, eta=None, chi_eff=None, **kwargs)
        Alternate constructor using total mass and either mass ratio or symmetric mass ratio.
    from_mchirp(mchirp, q=None, eta=None, chi_eff=None, **kwargs)
        Alternate constructor using chirp mass and either mass ratio or symmetric mass ratio.
    to_dimless_time(time_s)
        Converts physical time in seconds to dimensionless time.
    to_physical_time(tau)
        Converts dimensionless time to physical time in seconds.
    to_dimless_length(length_m)
        Converts physical length in meters to dimensionless length.
    to_physical_length(length_bar)
        Converts dimensionless length to physical length in meters.
    to_dimless_frequency(freq_hz)
        Converts physical frequency in Hz to dimensionless frequency.
    to_physical_frequency(freq_bar)
        Converts dimensionless frequency to physical frequency in Hz.
    """

    # Parameters of the binary system ----------------------------------------------

    mass_1 : float
    mass_2 : float

    lambda_1 : float = 0.0
    lambda_2 : float = 0.0

    ra  : float = 0.0
    dec : float = 0.0

    phase : float = 0.0
    psi   : float = 0.0

    geocent_time        : float | LIGOTimeGPS = 0.0
    luminosity_distance : float = 100.0

    spins : InitVar[JFrameSpins | LFrameSpins | None] = None

    j_frame : JFrameSpins = field(init=False, repr=True)
    l_frame : LFrameSpins = field(init=False, repr=True)


    # Post init chores -------------------------------------------------------------
    
    def __post_init__(self, spins: JFrameSpins | LFrameSpins | None) -> None :

        if not isinstance(self.geocent_time, LIGOTimeGPS) :
            object.__setattr__(self, "geocent_time", LIGOTimeGPS(self.geocent_time))

        if spins is None:
            object.__setattr__(self, "j_frame", JFrameSpins())
            object.__setattr__(self, "l_frame", LFrameSpins())

        elif isinstance(spins, JFrameSpins) :
            from gwmlt.utils.physics import convert_jframe_to_lframe
            l_struct = convert_jframe_to_lframe(self.mass_1, self.mass_2, spins)
            object.__setattr__(self, "j_frame", spins)
            object.__setattr__(self, "l_frame", l_struct)

        elif isinstance(spins, LFrameSpins) :
            from gwmlt.utils.physics import convert_lframe_to_jframe
            j_struct = convert_lframe_to_jframe(self.mass_1, self.mass_2, spins)
            object.__setattr__(self, "j_frame", j_struct)
            object.__setattr__(self, "l_frame", spins)

        else :
            raise TypeError(
                f"Invalid spins type '{type(spins).__name__}'. "
                "Must be JFrameSpins, LFrameSpins, or None."
            )
    

    # Alternate factory constructors -----------------------------------------------

    @classmethod
    def from_mtot(
        cls, 
        mtot: float, 
        q: float | None = None, 
        eta: float | None = None, 
        chi_eff: float | None = None, 
        **kwargs
    ) -> "BinarySystem" :
    
        """
        Initialize a BinarySystem (or child thereof) from total mass and 
        either mass ratio (q) or symmetric mass ratio (eta).

        Parameters
        ----------
        mtot : float
            Total mass of the binary in solar masses.
        q : float, optional
            Mass ratio (m_secondary / m_primary).
            Cannot be provided simultaneously with `eta`.
        eta : float, optional
            Symmetric mass ratio. Must be <= 0.25. 
            Cannot be provided simultaneously with `q`.
        chi_eff : float, optional
            Effective aligned spin parameter. If provided without explicit `spins`, 
            it assigns equal aligned spins (spin1z = spin2z = chi_eff) to both binaries.
        **kwargs
            Explicit spins and, if required, extrinsic parameters for the system.
            See the docstring of the base class for all accepted keyword arguments.
            
        Raises
        ------
        ValueError
            If both `q` and `eta` are provided.
            If `eta` > 0.25.
            If both `chi_eff` and explicit `spins` are provided in `kwargs`.

        Note
        ----
        If both `q` and `eta` are None, the function defaults to equal masses (q=1, eta=0.25).
        """

        if q is not None and eta is not None :
            raise ValueError("Cannot provide both `q` and `eta` simultaneously.")

        if chi_eff is not None :
            if kwargs.get("spins") is not None :
                raise ValueError("Cannot provide both `chi_eff` and explicit `spins`.")
            kwargs["spins"] = LFrameSpins(spin1z=chi_eff, spin2z=chi_eff)

        if eta is not None :
            if eta > 0.25 :
                raise ValueError(f"Symmetric mass ratio (eta) cannot exceed 0.25, got {eta}.")
            m1 = 0.5 * mtot * (1.0 + (1.0 - 4.0 * eta)**0.5)
            m2 = 0.5 * mtot * (1.0 - (1.0 - 4.0 * eta)**0.5)
        else :
            # Defaults to equal masses if neither q nor eta is provided
            q_val = 1.0 if q is None else q
            m1 = mtot / (1.0 + q_val)
            m2 = q_val * m1

        return cls(mass_1=m1, mass_2=m2, **kwargs)

    @classmethod
    def from_mchirp(
        cls, 
        mchirp: float, 
        q: float | None = None, 
        eta: float | None = None, 
        chi_eff: float | None = None, 
        **kwargs
    ) -> "BinarySystem" :
    
        """
        Initialize a BinarySystem (or child thereof) from chirp mass and 
        either mass ratio (q) or symmetric mass ratio (eta).

        Parameters
        ----------
        mchirp : float
            Chirp mass of the binary in solar masses.
        q : float, optional
            Mass ratio (m_secondary / m_primary).
            Cannot be provided simultaneously with `eta`.
        eta : float, optional
            Symmetric mass ratio. Must be <= 0.25. 
            Cannot be provided simultaneously with `q`.
        chi_eff : float, optional
            Effective aligned spin parameter. If provided without explicit `spins`, 
            it assigns equal aligned spins (spin1z = spin2z = chi_eff) to both binaries.
        **kwargs
            Explicit spins and, if required, extrinsic parameters for the system.
            See the docstring of the base class for all accepted keyword arguments.
            
        Raises
        ------
        ValueError
            If both `q` and `eta` are provided.
            If `eta` > 0.25.
            If both `chi_eff` and explicit `spins` are provided in `kwargs`.

        Note
        ----
        If both `q` and `eta` are None, the function defaults to equal masses (q=1, eta=0.25).
        """
        
        if q is not None and eta is not None :
            raise ValueError("Cannot provide both `q` and `eta` simultaneously.")

        if chi_eff is not None :
            if kwargs.get("spins") is not None :
                raise ValueError("Cannot provide both `chi_eff` and explicit `spins`.")
            kwargs["spins"] = LFrameSpins(spin1z=chi_eff, spin2z=chi_eff)

        if eta is not None :
            if eta > 0.25 :
                raise ValueError(f"Symmetric mass ratio (eta) cannot exceed 0.25, got {eta}.")
            mtot = mchirp * eta**(-0.6)
            m1 = 0.5 * mtot * (1.0 + (1.0 - 4.0 * eta)**0.5)
            m2 = 0.5 * mtot * (1.0 - (1.0 - 4.0 * eta)**0.5)
        else :
            # Defaults to equal masses if neither q nor eta is provided
            q_val = 1.0 if q is None else q
            m1 = (mchirp * (1.0 + q_val)**0.2) / (q_val**0.6)
            m2 = m1 * q_val

        return cls(mass_1=m1, mass_2=m2, **kwargs)
    

    # Binary properties ------------------------------------------------------------

    @property
    def mtot(self) -> float:
        """Total mass of the binary system in solar masses."""
        return self.mass_1 + self.mass_2
    
    @property
    def mchirp(self) -> float:
        """Chirp mass of the binary system in solar masses."""
        return (self.mass_1 * self.mass_2)**0.6 / self.mtot**0.2
    
    @property
    def dimless_masses(self) -> tuple[float, float]:
        """Returns the dimensionless masses (m1/M_tot, m2/M_tot)."""
        return self.mass_1 / self.mtot, self.mass_2 / self.mtot

    @property
    def q(self) -> float:
        """Mass ratio (m_secondary / m_primary). Can be > 1 if mass_2 is heavier."""
        return self.mass_2 / self.mass_1

    @property
    def eta(self) -> float:
        """Symmetric mass ratio of the binary system."""
        return (self.mass_1 * self.mass_2) / self.mtot**2

    @property
    def t_scale(self) -> float:
        """Characteristic time scale (GM_tot/c^3) of the binary system in seconds."""
        return (G_SI * self.mtot * MSUN_SI) / (C_SI**3)

    @property
    def l_scale(self) -> float:
        """Characteristic length scale (GM_tot/c^2) of the binary system in meters."""
        return (G_SI * self.mtot * MSUN_SI) / (C_SI**2)

    @property
    def dimless_spins_z(self) -> tuple[float, float]:
        """
        Returns the dimensionless z-angular momenta (S1z_bar, S2z_bar) 
        from the L-frame aligned spins (chi_z or spin_z) and dimensionless masses.
        """
        m1_bar, m2_bar = self.dimless_masses
        s1z_bar = self.l_frame.spin1z * (m1_bar**2)
        s2z_bar = self.l_frame.spin2z * (m2_bar**2)
        return s1z_bar, s2z_bar

    @property
    def chi_eff(self) -> float:
        """Effective aligned spin (weighted average of z spins)."""
        return (self.mass_1 * self.l_frame.spin1z + self.mass_2 * self.l_frame.spin2z) / self.mtot

    @property
    def chi_p(self) -> float:

        """
        Effective precession spin parameter chi_p.
        This characterizes the contribution of spin to the orbital precession
        in terms of the components of the spin vectors perpendicular to the 
        orbital angular momentum for the two binaries.

        References: 
            - https://arxiv.org/pdf/1408.1810,
            - https://cplberry.com/2020/04/18/gw190412/

        Notes
        -----
        The chi_p is obtained in the following way :
        1. If the input mass ratio (q) is greater than 1, it is inverted to 1/q to ensure q <= 1.
        2. chi_p is then calculated as the maximum between:
            - chi_1_perp: the spin perpendicular to the orbital plane for the primary black hole.
            - chi_p_v2: a weighted combination of chi_2_perp and the mass ratio q, accounting 
            for the influence of the lighter black hole's spin on precession.
        """

        # Since we need q <= 1 here
        m_heavy = max(self.mass_1, self.mass_2)
        m_light = min(self.mass_1, self.mass_2)
        q_leq = m_light / m_heavy

        # Extract in-plane spin magnitudes based on which mass is heavier
        if self.mass_1 >= self.mass_2:
            chi_heavy_perp = (self.l_frame.spin1x**2 + self.l_frame.spin1y**2)**0.5
            chi_light_perp = (self.l_frame.spin2x**2 + self.l_frame.spin2y**2)**0.5
        else:
            chi_heavy_perp = (self.l_frame.spin2x**2 + self.l_frame.spin2y**2)**0.5
            chi_light_perp = (self.l_frame.spin1x**2 + self.l_frame.spin1y**2)**0.5

        # Calculate chi_p using the q <= 1 formulation
        chi_p_light_weighted = q_leq * ((4.0 * q_leq + 3.0) / (4.0 + 3.0 * q_leq)) * chi_light_perp
        return max(chi_heavy_perp, chi_p_light_weighted)
    

    # Dimensionless and physical interconversion methods ---------------------------

    def to_dimless_time(self, time_s: float | np.ndarray) -> float | np.ndarray:
        """Convert physical time in seconds to dimensionless time."""
        return time_s / self.t_scale
    
    def to_physical_time(self, tau: float | np.ndarray) -> float | np.ndarray:
        """Convert dimensionless time to physical time in seconds."""
        return tau * self.t_scale

    def to_dimless_length(self, length_m: float | np.ndarray) -> float | np.ndarray:
        """Convert physical length in meters to dimensionless length."""
        return length_m / self.l_scale
    
    def to_physical_length(self, length_bar: float | np.ndarray) -> float | np.ndarray:
        """Convert dimensionless length to physical length in meters."""
        return length_bar * self.l_scale

    def to_dimless_frequency(self, freq_hz: float | np.ndarray) -> float | np.ndarray:
        """Convert physical frequency in Hz to dimensionless frequency (Mf)."""
        return freq_hz * self.t_scale
    
    def to_physical_frequency(self, freq_bar: float | np.ndarray) -> float | np.ndarray:
        """Convert dimensionless frequency (Mf) to physical frequency in Hz."""
        return freq_bar / self.t_scale


# Source kinds ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BBHSystem(BinarySystem) :

    """
    Parameters of a binary black hole system (BBH System) source. Inherits from `BinarySystem`.

    This struct contains both intrinsic parameters (masses, spins, etc.) 
    and extrinsic parameters (sky location, orientation, time, etc.) of 
    a BBH source both in L-frame and J-frame coordinates.

    Parameters
    ----------
    mass_1 : float
        Detector-frame mass of the first binary in solar masses.
        This is referred to as the primary (conventionally heavier, but ordering is not enforced).
    mass_2 : float
        Detector-frame mass of the second binary in solar masses.
        This is referred to as the secondary (conventionally lighter, but ordering is not enforced).
    ra : float, optional
        Right ascension of the source in radians.
        Default is 0.0.
    dec : float, optional
        Declination of the source in radians.
        Default is 0.0.
    phase : float, optional
        Coalescence orbital phase of the binary system in radians. 
        Default is 0.0.
    psi : float, optional
        Polarization angle of the gravitational wave in radians.
        Default is 0.0.
    geocent_time : float or LIGOTimeGPS, optional
        The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds. 
        Default is 0.0.
    luminosity_distance : float, optional
        Luminosity distance to the source in megaparsec (Mpc).
        Default is 100.0.
    spins : JFrameSpins, or LFrameSpins, or None, optional
        Dataclass containing spin parameters. Default is None, which implies non-spinning binaries.
    
    Attributes
    ----------
    j_frame : JFrameSpins
        Dataclass containing J-frame spin parameters.
    l_frame : LFrameSpins
        Dataclass containing L-frame spin parameters.
    *inherited :
        Inherits all properties (mtot, mchirp, q, chi_eff, etc.) from the parent `BinarySystem` class.

    Methods
    -------
    from_mtot, from_mchirp
        Inherits all alternative factory constructors from `BinarySystem`. Returns a `BBHSystem`.
    to_dimless_*, to_physical_*
        Inherits all interconversion methods from `BinarySystem`.
    """
    
    # Enforce lambda values of 0.0 for BBH systems
    
    lambda_1: float = field(default=0.0, init=False)
    lambda_2: float = field(default=0.0, init=False)