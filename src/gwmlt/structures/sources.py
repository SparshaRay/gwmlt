"""
Structs for GW Source Parameters
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class BBHSystem :

    """
    Parameters of a Binary Black Hole (BBH) source.

    This struct contains both intrinsic parameters (masses, spins) 
    and extrinsic parameters (sky location, orientation, time) of 
    a compact binary coalescence in J-frame coordinates.

    Parameters
    ----------
    mass_1 : float
        Detector-frame mass of the primary (heavier) compact object in solar masses.
    mass_2 : float
        Detector-frame mass of the secondary (lighter) compact object in solar masses.
    a_1 : float, optional
        Dimensionless spin magnitude of the primary compact object.
        Must lie in the range [0, 1].
        Default is 0.0.
    a_2 : float, optional
        Dimensionless spin magnitude of the secondary compact object.
        Must lie in the range [0, 1].
        Default is 0.0.
    tilt_1 : float, optional
        Zenith angle between the primary spin vector and the orbital angular momentum 
        vector in radians, defined at the reference frequency.
        Default is 0.0.
    tilt_2 : float, optional
        Zenith angle between the secondary spin vector and the orbital angular momentum 
        vector in radians, defined at the reference frequency. 
        Default is 0.0.
    phi_12 : float, optional
        Difference between the azimuthal angles of the two spin vectors in the 
        plane perpendicular to the orbital angular momentum in radians.
        Default is 0.0.
    phi_jl : float, optional
        Azimuthal angle of the orbital angular momentum on its cone of precession 
        around the total angular momentum J in radians.
        Default is 0.0.
    ra : float, optional
        Right ascension of the source in radians.
        Default is 0.0.
    dec : float, optional
        Declination of the source in radians.
        Default is 0.0.
    theta_jn : float, optional
        Angle between the total angular momentum J and the line of sight vector in radians.
        Default is 0.0.
    phase : float, optional
        Coalescence orbital phase of the binary system in radians at the reference frequency.
        Default is 0.0.
    psi : float, optional
        Polarization angle of the gravitational wave in radians.
        Default is 0.0.
    geocent_time : float, optional
        The geocentric trigger time (coalescence time at the center of the Earth) in GPS seconds.
        Default is 0.0.

    Notes
    -----
    Spins and orientation angles are defined in the J-frame convention matching Bilby parameter estimation specifications.
    Downstream wrappers will project these into L-frame alignments as required by PyCBC/LALSimulation waveform models.
    (i.e., cartesian spin components spin1x, spin1y, spin1z, spin2x, spin2y, spin2z)
    """

    mass_1 : float
    mass_2 : float

    a_1    : float = 0.0
    a_2    : float = 0.0
    tilt_1 : float = 0.0
    tilt_2 : float = 0.0
    phi_12 : float = 0.0
    phi_jl : float = 0.0

    ra  : float = 0.0
    dec : float = 0.0

    theta_jn : float = 0.0
    phase    : float = 0.0
    psi      : float = 0.0

    geocent_time : float = 0.0