"""
Configurations for the modules. 
"""

import os
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", default=Path(__file__).resolve().parents[2]))


# Global waveform approximant configs -------------------------------------------------------------

_DEFAULT_IMRPHENOMX_MODES : None | list[tuple[int, int]] = None    # Use all available modes in IMRPhenomXO4a
_DEFAULT_TEOBRESUMS_MODES : list[int] = [1]                        # Use only the (2, 2) mode in TEOBResumS
_DEFAULT_ECC_FREQ         : int = 3                                # Use orbit-average frequency for TEOBResumS
_DEFAULT_F_REF            : float = 20.0                           # Reference frequency (Hz) for spin values
_DEFAULT_TD_WF_GEN_SRATE  : int = 4096                             # Internal frequency (Hz) for waveform generation 

def set_imrphenomx_modes(modes : None | list[tuple[int, int]]) :
    """
    Set the default modes for IMRPhenomXO4a waveform approximant.
    To use all available modes, set to None. Otherwise, a list of tuples of (l, m) modes is expected.
    Default global value is None.
    """
    global _DEFAULT_IMRPHENOMX_MODES
    print(f"Setting default IMRPhenomXO4a modes from {_DEFAULT_IMRPHENOMX_MODES} to {modes}")
    _DEFAULT_IMRPHENOMX_MODES = modes

def set_teobresums_modes(modes : list[int]) :
    """
    Set the default modes for TEOBResumS waveform approximant.
    A list of linear indices k is expected, where k = (l*(l-1)/2 + m-2) for the (l, m) mode. 
    Default global value is [1], i.e. only the (2, 2) mode is used.
    """
    global _DEFAULT_TEOBRESUMS_MODES
    print(f"Setting default TEOBResumS modes from {_DEFAULT_TEOBRESUMS_MODES} to {modes}")
    _DEFAULT_TEOBRESUMS_MODES = modes

def set_teobresums_ecc_freq(freq : int) :
    """
    Set the default ecc_freq parameter for TEOBResumS waveform approximant. 
    This parameter is used to define the starting frequency.
    Can be set to 0, 1, 2, or 3. These correspond to the periastron frequency, 
    average of periastron & apastron frequencies, apastron frequency, 
    and orbit-average frequency respectively. Refer to arxiv:2302.11257 for implications.
    Default global value is 3, i.e. orbit-average frequency.
    """
    global _DEFAULT_ECC_FREQ
    print(f"Setting default TEOBResumS ecc_freq from {_DEFAULT_ECC_FREQ} to {freq}")
    _DEFAULT_ECC_FREQ = freq

def set_f_ref(freq : float) :
    """
    Set the default reference frequency (in Hz) for spin values.
    Default global value is 20.0 Hz, following Bilby values.
    Changing this has no effect on non-precessing systems.
    """
    global _DEFAULT_F_REF
    print(f"Setting default reference frequency from {_DEFAULT_F_REF} to {freq}")
    _DEFAULT_F_REF = freq

def set_td_wf_gen_srate(srate : int) :
    """
    Set the default internal frequency (in Hz) for time-domain waveform generation.
    All time-domain waveforms are generated at this sample rate and then downsampled to the requested sample rate.
    Default global value is 4096 Hz.
    """
    global _DEFAULT_TD_WF_GEN_SRATE
    print(f"Setting default time-domain waveform generation sample rate from {_DEFAULT_TD_WF_GEN_SRATE} to {srate}")
    _DEFAULT_TD_WF_GEN_SRATE = srate


# Noise configs -----------------------------------------------------------------------------------

NOISE_DATA_DIR = PROJECT_ROOT / "data" / "noise"

DEFAULT_NOISE_PROFILES = {

    # Real glitch-free event-free detector noise timeseries 
    "O4_real": {
        "H1": NOISE_DATA_DIR / "real_noise" / "H1_noise.hdf5",
        "L1": NOISE_DATA_DIR / "real_noise" / "L1_noise.hdf5",
    },

    # O4 High
    "O4_gaus": {
        "H1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O4_high_psd.npz",
        "L1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O4_high_psd.npz",
        "V1": NOISE_DATA_DIR / "detector_PSDs" / "aVirgo_O4_high_psd.npz",
    },

    # O5b
    "O5b_gaus": {
        "H1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O5b_psd.npz",
        "L1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O5b_psd.npz",
        "V1": NOISE_DATA_DIR / "detector_PSDs" / "aVirgo_O5_high_psd.npz",
    },

    # IndIGO-D High
    "S1_gaus": {
        "IndIGO-D": NOISE_DATA_DIR / "detector_PSDs" / "IndIGO-D_S2_psd.npz",
    },

}