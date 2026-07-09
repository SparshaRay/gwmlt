"""
Structs to Define Waveform Type
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Quasicircular :

    """
    Quasi-circular binary system.
    No precession (aligned L and S), no lensing effects, and no eccentricity effects are considered.

    Attributes
    ----------
    wf_approximant : str
        The LAL/PyCBC waveform model string.
        Default is "IMRPhenomXO4a".
    override_pars : dict
        A dictionary of parameters to override in waveform generation.
        Default is {"spin1x": 0.0, "spin1y": 0.0, "spin2x": 0.0, "spin2y": 0.0},
        i.e. no precession is considered.
    """

    wf_approximant : str = field(default="IMRPhenomXO4a", init=False)
    override_pars  : dict[str, float] = field(init=False, default_factory=lambda: {
        "spin1x" : 0.0,
        "spin1y" : 0.0,
        "spin2x" : 0.0,
        "spin2y" : 0.0,
    })


@dataclass(frozen=True)
class Eccentric :

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
    override_pars : dict
        A dictionary of parameters to override in waveform generation.
        Default is {"spin1x": 0.0, "spin1y": 0.0, "spin2x": 0.0, "spin2y": 0.0},
        i.e. no precession is considered.
    """

    eccentricity   : float
    anomaly        : float = 0.0
    wf_approximant : str = field(default="teobresums", init=False)
    override_pars  : dict[str, float] = field(init=False, default_factory=lambda: {
        "spin1x" : 0.0,
        "spin1y" : 0.0,
        "spin2x" : 0.0,
        "spin2y" : 0.0,
    })


@dataclass(frozen=True)
class Lensed :

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
    
    Attributes
    ----------
    wf_approximant : str
        The eccentric waveform model string.
        Default is "teobresums".
    override_pars : dict
        A dictionary of parameters to override in waveform generation.
        Default is {"spin1x": 0.0, "spin1y": 0.0, "spin2x": 0.0, "spin2y": 0.0},
        i.e. no precession is considered.
    """
    
    m_lens: float
    y_lens: float
    z_lens: float = 0.0
    wf_approximant: str = field(default="IMRPhenomXO4a", init=False)
    override_pars: dict[str, float] = field(init=False, default_factory=lambda: {
        "spin1x": 0.0,
        "spin1y": 0.0,
        "spin2x": 0.0,
        "spin2y": 0.0,
    })