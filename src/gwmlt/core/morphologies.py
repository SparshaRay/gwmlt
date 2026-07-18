"""
Structs to Define Waveform Kinds
"""

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable


# Base classes for waveform morphologies ----------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Morphology :
    """
    Base class for waveform morphologies. This is primarily used for type hinting.
    The `wf_approximant` field is a necessary attribute for any given morphology. 
    It must be set downstream or post-hoc setattr.

    Attributes
    ----------
    wf_approximant : str
        The LAL/PyCBC waveform approximant string.
    allow_precession : bool
        A boolean indicating whether precession is allowed in the waveform generation.
        If False, the in-plane spin components are set to zero.
    """
    wf_approximant   : str  = field(init=False)
    allow_precession : bool = field(init=False)


# Structural protocols

@runtime_checkable
class EccentricProtocol(Protocol) :
    """Contract for morphologies that include eccentricity."""
    eccentricity : float
    anomaly      : float

@runtime_checkable
class LensedProtocol(Protocol) :
    """Contract for morphologies that include microlensing parameters."""
    m_lens : float
    y_lens : float
    z_lens : float


# Morphology kinds --------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Quasicircular(Morphology) :

    """
    Quasi-circular binary system.
    No precession (aligned L and S), no lensing effects, and no eccentricity effects are considered.

    Attributes
    ----------
    wf_approximant : str
        The LAL/PyCBC waveform model string.
        Default is "IMRPhenomXO4a".
    allow_precession : bool
        A boolean indicating whether precession is allowed in the waveform generation.
        If False, the in-plane spin components are set to zero. Default is False.
    """

    wf_approximant   : str  = field(default="IMRPhenomXO4a", init=False)
    allow_precession : bool = field(default=False, init=False)


@dataclass(frozen=True)
class Eccentric(Morphology, EccentricProtocol) :

    """
    Non-zero orbital eccentricity binary system.
    No precession (aligned L and S), and no lensing effects are considered.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity of the binary.
        Evaluated at the initial condition for LIGO band 
        and at the eccentricity reference frequency (ecc_f_ref) for the decihertz band.
    anomaly : float, optional
        True anomaly of the binary system in radians at the initial condition.
        Default is 0.0.

    Attributes
    ----------
    wf_approximant : str
        The eccentric waveform model string.
        Default is "teobresums".
    allow_precession : bool
        A boolean indicating whether precession is allowed in the waveform generation.
        If False, the in-plane spin components are set to zero. Default is False.
    """

    eccentricity     : float
    anomaly          : float = 0.0
    wf_approximant   : str   = field(default="teobresums", init=False)
    allow_precession : bool  = field(default=False, init=False)


@dataclass(frozen=True)
class Lensed(Morphology, LensedProtocol) :

    """
    Point lens effects on quasi-circular GW signal.
    No precession (aligned L and S), and no eccentricity effects are considered.

    Parameters
    ----------
    m_lens : float
        Mass of the gravitational point lens in solar masses.
    y_lens : float
        Dimensionless impact parameter between the lens and the source.
    z_lens : float, optional
        Redshift of the lens system.
        Default is 0.0 (equivalent to treating m_lens as the redshifted lens mass).

    Attributes
    ----------
    wf_approximant : str
        The eccentric waveform model string.
        Default is "IMRPhenomXO4a".
    allow_precession : bool
        A boolean indicating whether precession is allowed in the waveform generation.
        If False, the in-plane spin components are set to zero. Default is False.
    """
    
    m_lens : float
    y_lens : float
    z_lens : float = 0.0
    wf_approximant   : str  = field(default="IMRPhenomXO4a", init=False)
    allow_precession : bool = field(default=False, init=False)


# Add more morphologies here as needed, e.g. :
# @dataclass(frozen=True)
# class LensedEccentricPrecessing(Morphology, LensedProtocol, EccentricProtocol) :
#     wf_approximant : str = field(default="teobresums", init=False)
#     eccentricity   : float
#     anomaly        : float = 0.0
#     m_lens : float
#     y_lens : float
#     z_lens : float = 0.0
#     allow_precession : bool = field(default=True, init=False)