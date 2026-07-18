"""
Structs for GW Source Parameters
"""

from dataclasses import dataclass, field, InitVar
from lal import LIGOTimeGPS


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
    """

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


# Source kinds ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BBHSystem(BinarySystem) :

    """
    Parameters of a binary black hole system (BBH System) source.

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
    """
    
    # Enforce lambda values of 0.0 for BBH systems
    
    lambda_1: float = field(default=0.0, init=False)
    lambda_2: float = field(default=0.0, init=False)