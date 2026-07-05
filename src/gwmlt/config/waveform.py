from dataclasses import dataclass, field

@dataclass
class WaveformConfig :
    """
    Configuration for waveform generation.
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
    The reference frequency (in Hz) for spin values.
    Changing this has no effect on non-precessing systems.
    Default global value is 20.0 Hz, following Bilby values.
    """

    td_wf_gen_srate : int = 4096
    """
    The frequency (in Hz) for time-domain waveform generation.
    All time-domain waveforms are generated at this sample rate 
    and then downsampled to the requested sample rate.
    Default global value is 4096 Hz.
    """