"""
Structs to Define Waveform Type
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Quasicircular:
    """
    Quasi-circular binary system.
    No precession (aligned L and S), no lensing effects, and no eccentricity effects are considered.

    Parameters
    ----------
    wf_approximant : str, optional
        The LAL/PyCBC waveform model string.
        Default is "IMRPhenomXO4a".
    """
    wf_approximant : str = "IMRPhenomXO4a"


@dataclass(frozen=True)
class Eccentric:
    """
    Non-zero orbital eccentricity binary system.
    No precession (aligned L and S), and no lensing effects are considered.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity of the binary evaluated at the initial condition.
    anomaly : float, optional
        True anomaly of the binary system in radians at the initial condition.
        Default is 0.0.
    wf_approximant : str, optional
        The eccentric waveform model string.
        Default is "teobresums".
    """
    eccentricity   : float
    anomaly        : float = 0.0
    wf_approximant : str = "teobresums"


@dataclass(frozen=True)
class Lensed:
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
        Default is 0.0 (equivalent to treating `m_lens` as the redshifted lens mass).
    wf_approximant : str, optional
        The waveform model string for the unlensed source signal.
        Default is "IMRPhenomXO4a".
    """
    m_lens: float
    y_lens: float
    z_lens: float = 0.0
    wf_approximant: str = "IMRPhenomXO4a"