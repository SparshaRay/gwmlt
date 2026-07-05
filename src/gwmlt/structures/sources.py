"""
Structs for GW Source Parameters
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class _JFrameSpins :

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
class _LFrameSpins :

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
    geocent_time : float, optional
        The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds. 
        Default is 0.0.
    luminosity_distance : float, optional
        Luminosity distance to the source in megaparsecs (Mpc).
        Default is 100.0.
    j_frame : _JFrameSpins
        Struct containing J-frame spin parameters.
    l_frame : _LFrameSpins
        Struct containing L-frame spin parameters.
    """

    mass_1 : float
    mass_2 : float

    ra     : float = 0.0
    dec    : float = 0.0

    phase  : float = 0.0
    psi    : float = 0.0

    geocent_time        : float = 0.0
    luminosity_distance : float = 100.0

    j_frame : _JFrameSpins = _JFrameSpins()
    l_frame : _LFrameSpins = _LFrameSpins()

    # Initialization in J-frame --------------------------------
    @classmethod
    def from_jframe(
        cls,

        mass_1 : float,
        mass_2 : float,

        a_1      : float = 0.0,
        a_2      : float = 0.0,
        tilt_1   : float = 0.0,
        tilt_2   : float = 0.0,
        phi_12   : float = 0.0,
        phi_jl   : float = 0.0,
        theta_jn : float = 0.0,

        ra       : float = 0.0,
        dec      : float = 0.0,

        phase    : float = 0.0,
        psi      : float = 0.0,

        geocent_time        : float = 0.0,
        luminosity_distance : float = 100.0,

    ) -> "BBHSystem" :
        
        """
        Construct a BBHSystem specifying J-frame spin coordinates.
        Automatically populates the L-frame spin values too.

        Parameters
        ----------
        mass_1 : float
            Detector-frame mass of the first binary in solar masses.
            This is referred to as the primary (heavier), but the ordering is not enforced.
        mass_2 : float
            Detector-frame mass of the second binary in solar masses.
            This is referred to as the secondary (lighter), but the ordering is not enforced.
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
        geocent_time : float, optional
            The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds. 
            Default is 0.0.
        luminosity_distance : float, optional
            Luminosity distance to the source in megaparsecs (Mpc).
            Default is 100.0.

        Returns
        -------
        BBHSystem
        """

        j_struct = _JFrameSpins(
            a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
            phi_12=phi_12, phi_jl=phi_jl, theta_jn=theta_jn
        )

        from ..utils.physics import _convert_jframe_to_lframe
        l_struct = _convert_jframe_to_lframe(mass_1, mass_2, j_struct)

        return cls(
            mass_1=mass_1, mass_2=mass_2, ra=ra, dec=dec, phase=phase, psi=psi,
            geocent_time=geocent_time, luminosity_distance=luminosity_distance,
            j_frame=j_struct, l_frame=l_struct
        )

    # Initialization in L-frame --------------------------------
    @classmethod
    def from_lframe(
        cls,

        mass_1 : float,
        mass_2 : float,

        spin1x      : float = 0.0,
        spin1y      : float = 0.0,
        spin1z      : float = 0.0,
        spin2x      : float = 0.0,
        spin2y      : float = 0.0,
        spin2z      : float = 0.0,
        inclination : float = 0.0,

        ra          : float = 0.0,
        dec         : float = 0.0,

        phase       : float = 0.0,
        psi         : float = 0.0,

        geocent_time        : float = 0.0,
        luminosity_distance : float = 100.0,
        
    ) -> "BBHSystem" :
        
        """
        Construct a complete system specifying Cartesian L0-frame spin components natively.

        Automatically performs coordinate frame transformations under the hood to calculate and 
        populate the alternative J-frame orbital angles and magnitudes.

        Parameters
        ----------
        mass_1 : float
            Detector-frame mass of the first binary in solar masses.
            This is referred to as the primary (heavier), but the ordering is not enforced.
        mass_2 : float
            Detector-frame mass of the second binary in solar masses.
            This is referred to as the secondary (lighter), but the ordering is not enforced.
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
        geocent_time : float, optional
            The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds. 
            Default is 0.0.
        luminosity_distance : float, optional
            Luminosity distance to the source in megaparsecs (Mpc).
            Default is 100.0.

        Returns
        -------
        BBHSystem
        """

        l_struct = _LFrameSpins(
            inclination=inclination,
            spin1x=spin1x, spin1y=spin1y, spin1z=spin1z,
            spin2x=spin2x, spin2y=spin2y, spin2z=spin2z
        )

        from ..utils.physics import _convert_lframe_to_jframe
        j_struct = _convert_lframe_to_jframe(mass_1, mass_2, l_struct)

        return cls(
            mass_1=mass_1, mass_2=mass_2, ra=ra, dec=dec, phase=phase, psi=psi,
            geocent_time=geocent_time, luminosity_distance=luminosity_distance,
            j_frame=j_struct, l_frame=l_struct
        )