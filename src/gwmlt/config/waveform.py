"""
Waveform Generation Configurations
"""

from dataclasses import dataclass, field

@dataclass(frozen=True)
class WaveformConfig :
    """
    Configurations for waveform generation.
    """

    imrphenomx_modes : list[tuple[int, int]] | None = None
    """
    Modes to use for IMRPhenomXO4a waveform approximant.
    Set to None to use all available modes.
    Otherwise, a list of tuples of (l, m) modes is expected.
    Default global value is None.
    """

    teobresums_modes : list[int] = field(default_factory=lambda : [1])
    """
    Modes to use for TEOBResumS waveform approximant.
    A list of linear integer indices k is expected, 
    where k = (l*(l-1)/2 + m-2) for the (l, m) mode. 
    Default global value is [1], i.e. only the (2, 2) mode is used.
    """

    ecc_freq : int = 3
    """
    The ecc_freq parameter for TEOBResumS waveform approximant. 
    This parameter is used to define the starting frequency.
    Can be set to 0, 1, 2, or 3. These correspond to the periastron frequency, 
    average of periastron & apastron frequencies, apastron frequency, 
    and orbit-average frequency respectively. Refer to arxiv:2302.11257 for implications.
    Default global value is 3, i.e. orbit-average frequency.
    """

    f_ref : float = 20.0
    """
    The reference frequency (in Hz) for spin values for IMRPhenomXO4a.
    This parameter is automatically set to the f_lower value for TEOBResumS Dali.
    Changing this has no effect on non-precessing systems beyond phase shifts.
    Default global value is 20.0 Hz, following Bilby values.
    """