"""
Structs for GW Source Parameters
"""

from __future__ import annotations
from dataclasses import dataclass, field

from lal import LIGOTimeGPS


@dataclass(frozen=True)
class BBHSystem :

    """
    Parameters of a Binary Black Hole (BBH) source.

    This struct contains both intrinsic parameters (masses, spins) 
    and extrinsic parameters (sky location, orientation, time) of 
    a compact binary coalescence both in L-frame and J-frame coordinates.

    Parameters
    ----------
    mass_1 : float
        Detector-frame mass of the first binary in solar masses.
        This is referred to as the primary (heavier), but the ordering is not enforced.
    mass_2 : float
        Detector-frame mass of the second binary in solar masses.
        This is referred to as the secondary (lighter), but the ordering is not enforced.
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
        Luminosity distance to the source in megaparsecs (Mpc).
        Default is 100.0.
    spins : BBHSystem.JFrameSpins, or BBHSystem.LFrameSpins, or None, optional
        Struct containing spin parameters. Default is None, which implies non-spinning binaries.
    
    Attributes
    ----------
    j_frame : BBHSystem.JFrameSpins
        Struct containing J-frame spin parameters.
    l_frame : BBHSystem.LFrameSpins
        Struct containing L-frame spin parameters.
    """

    mass_1 : float
    mass_2 : float

    ra     : float = 0.0
    dec    : float = 0.0

    phase  : float = 0.0
    psi    : float = 0.0

    geocent_time        : float | LIGOTimeGPS = 0.0
    luminosity_distance : float = 100.0

    spins : "BBHSystem.JFrameSpins" | "BBHSystem.LFrameSpins" | None = field(default=None, repr=False)

    j_frame: "BBHSystem.JFrameSpins" = field(default_factory=lambda: BBHSystem.JFrameSpins(), init=False, repr=True)
    l_frame: "BBHSystem.LFrameSpins" = field(default_factory=lambda: BBHSystem.LFrameSpins(), init=False, repr=True)

    @dataclass(frozen=True)
    class JFrameSpins :

        """
        Spins and orientation in the J-frame (used by Bilby).

        Parameters
        ----------
        a_1 : float, optional
            Dimensionless spin magnitude of the primary binary. 
            Must lie in the range [0, 1]. Default is 0.0.
        a_2 : float, optional
            Dimensionless spin magnitude of the secondary binary. 
            Must lie in the range [0, 1]. Default is 0.0.
        tilt_1 : float, optional
            Zenith angle between the primary spin vector and the orbital angular 
            momentum vector in radians, defined at the reference frequency. 
            Default is 0.0.
        tilt_2 : float, optional
            Zenith angle between the secondary spin vector and the orbital angular 
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
            The x-component of the primary dimensionless spin vector. Default is 0.0.
        spin1y : float, optional
            The y-component of the primary dimensionless spin vector. Default is 0.0.
        spin1z : float, optional
            The z-component of the primary dimensionless spin vector (aligned spin).
            Default is 0.0.
        spin2x : float, optional
            The x-component of the secondary dimensionless spin vector. Default is 0.0.
        spin2y : float, optional
            The y-component of the secondary dimensionless spin vector. Default is 0.0.
        spin2z : float, optional
            The z-component of the secondary dimensionless spin vector (aligned spin).
            Default is 0.0.
        """

        inclination : float = 0.0

        spin1x      : float = 0.0
        spin1y      : float = 0.0
        spin1z      : float = 0.0

        spin2x      : float = 0.0
        spin2y      : float = 0.0
        spin2z      : float = 0.0
    

    def __post_init__(self) -> None :

        if not isinstance(self.geocent_time, LIGOTimeGPS) :
            object.__setattr__(self, "geocent_time", LIGOTimeGPS(self.geocent_time))

        if self.spins is None :
            object.__setattr__(self, "j_frame", self.JFrameSpins())
            object.__setattr__(self, "l_frame", self.LFrameSpins())

        elif isinstance(self.spins, self.JFrameSpins) :
            from ..utils.general import _convert_jframe_to_lframe
            l_struct = _convert_jframe_to_lframe(self.mass_1, self.mass_2, self.spins)
            object.__setattr__(self, "j_frame", self.spins)
            object.__setattr__(self, "l_frame", l_struct)

        elif isinstance(self.spins, self.LFrameSpins) :
            from ..utils.general import _convert_lframe_to_jframe
            j_struct = _convert_lframe_to_jframe(self.mass_1, self.mass_2, self.spins)
            object.__setattr__(self, "j_frame", j_struct)
            object.__setattr__(self, "l_frame", self.spins)

        else :
            raise TypeError(
                f"Invalid spins type '{type(self.spins).__name__}'. "
                "Must be BBHSystem.JFrameSpins, BBHSystem.LFrameSpins, or None."
            )

        object.__setattr__(self, "spins", None)