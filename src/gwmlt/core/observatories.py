"""
Structs Defining Observatory Networks and Configurations.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from gwmlt.config import config


# Detector Networks -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True, kw_only=True)
class Observatory :

    """
    Base class for observatory network configurations.
    
    Parameters
    ----------
    noise_psds : list of str or Path or None
        Target noise psd (e.g., [None, 'O4_gaus', Path('path/to/psd.npz')] etc.)
        If Path, it should point to a npz file with 'freq' and 'psd' arrays,
        or an hdf5 file with noise timeseries.
    detector_list : list of str, or None, optional
        Target detectors. Default is None,
        i.e., select all available detectors for the specified PSD tag.
    f_low : float
        Lower bound of frequency range for all analyses (in Hz).
    sample_rate : float
        Sampling rate of the timeseries data (in Hz). 
        Sets upper bound of frequency range for all analyses at the Nyquist frequency.
    wf_gen_srate : float
        Internal sampling rate for the TD waveform generation (in Hz).
    
    Attributes
    ----------
    resolved_paths : dict
        A map of detector names to their respective noise file paths
    active_detectors : list of str
        A list of the active detectors in the observatory network.
    
    Note
    ----
    If `detector_list` is None, it defaults to all detectors available for the specified noise psd.
    Otherwise, `noise_psds` and `detector_list` must be broadcastable to a common shape.
    The values which must be specified in the children have no defaults here.
    """

    noise_psds    : list[str | Path | None]
    detector_list : list[str] | None = None
    f_low         : float
    sample_rate   : float
    wf_gen_srate  : float

    resolved_paths : dict[str, Path] = field(default_factory=dict, init=False, repr=True)

    def __post_init__(self) -> None:
        paths_map = _resolve_paths(self.noise_psds, self.detector_list)
        object.__setattr__(self, "resolved_paths", paths_map)

    @property
    def active_detectors(self) -> list[str]:
        return list(self.resolved_paths.keys())


@dataclass(frozen=True)
class GroundBased(Observatory) :

    """
    Ground-based detector network configurations.

    Parameters
    ----------
    noise_psds : list of str or Path or None
        Target noise psd (e.g., [None, 'O4_gaus', Path('path/to/psd.npz')] etc.)
        If Path, it should point to a npz file with 'freq' and 'psd' arrays,
        or an hdf5 file with noise timeseries.
    detector_list : list of str, or None, optional
        Target detectors (e.g., ['H1', 'L1'] etc.). Default is None, 
        i.e., select all available detectors for the specified PSD tag.
    f_low : float, optional
        Lower bound of frequency range for all analyses (in Hz).
        Also sets reference frequency and starting frequency for waveform generation.
        Default is 20.0
    sample_rate : float, optional
        Sampling rate of the timeseries data (in Hz).
        Sets upper bound of frequency range for all analyses at the Nyquist frequency.
        Default is 4096.0
    wf_gen_srate : float, optional
        Internal sampling rate for the TD waveform generation (in Hz).
        Should be kept at sample_rate for ground-based detectors.
        Default is 4096.0

    Attributes
    ----------
    resolved_paths : dict
        A map of detector names to their respective noise file paths: 
        { 'H1': Path(...), 'L1': Path(...) }
    active_detectors : list of str
        A list of the active detectors in the observatory network.

    Note
    ----
    If `detector_list` is None, it defaults to all detectors available for the specified noise psd.
    Otherwise, `noise_psds` and `detector_list` must be broadcastable to a common shape.
    """

    f_low         : float = 20.000
    sample_rate   : float = 4096.0
    wf_gen_srate  : float = 4096.0


@dataclass(frozen=True)
class Decihertz(Observatory) :

    """
    Space-based decihertz detector configurations.

    Parameters
    ----------
    wf_duration : float
        Desired duration of the waveform (in seconds).
    noise_psds : list of str or Path or None
        Target noise psd (e.g., ['S1_gaus'] etc.)
        If Path, it should point to a npz file with 'freq' and 'psd' arrays.
    detector_list : list of str, or None, optional
        Target detectors (e.g., ['IndIGO-D'] etc.). Default is None, i.e., 
        select all available detectors for the specified PSD tag.
    ecc_f_ref : float, optional
        Reference frequency for eccentricity. Relevant for eccentric waveforms only.
        This is different from the f_ref which set globally according to bilby PE settings.
        Default is 20.0
    f_low : float, optional
        Lower bound of frequency range for all analyses (in Hz).
        This is not the same f_start used for waveform generation, 
        which is automatically determined from the desired waveform duration.
        Default is 1.0
    sample_rate : float, optional
        Sampling rate of the timeseries data (in Hz).
        Sets upper bound of frequency range for all analyses at the Nyquist frequency.
        Default is 20.0
    wf_gen_srate : float, optional
        Internal sampling rate for the TD waveform generation (in Hz).
        Should be a power-of-two multiple of the sample_rate for butterworth downsampling.
        This is necessary to prevent aliasing artifacts in the projected strain.
        Default is 2560.0

    Attributes
    ----------
    resolved_paths : dict
        A map of detector names to their respective noise file paths: 
        { 'IndIGO-D': Path(...) }
    active_detectors : list of str
        A list of the active detectors in the observatory network.

    Note
    ----
    If `detector_list` is None, it defaults to all detectors available for the specified noise psd.
    Otherwise, `noise_psds` and `detector_list` must be broadcastable to a common shape.
    """

    wf_duration  : float = field(kw_only=True)
    ecc_f_ref    : float = 20.0
    f_low        : float = 1.00
    sample_rate  : float = 20.0
    wf_gen_srate : float = 2560.0


# Helper function to resolve noise psd paths for given detectors ----------------------------------
# -------------------------------------------------------------------------------------------------

def _resolve_paths(
    noise_psds: list[str | Path | None], 
    detector_list: list[str] | None
) -> dict[str, Path] :

    """
    Resolve the noise psd paths for the given detectors.

    Parameters
    ----------
    noise_psds : list of str or Path or None
    detector_list : list of str, or None

    Returns
    -------
    dict
        A map of detector names to their resolved noise file paths.
    """

    # Wrap single values in lists for user convenience, because python does not care about type hints anyways
    # Bad practice, but oh well, saves a second
    noise_psds          = noise_psds    if isinstance(noise_psds,    list) else [noise_psds]
    available_detectors = detector_list if isinstance(detector_list, list) else [detector_list]

    default_noise_psds = config.noise.noise_psds

    # Handle the case when detector_list is None, defaults to all detectors available in the specified noise psd
    if available_detectors == [None] :
        if (len(noise_psds) == 1) and (noise_psds[0] in default_noise_psds) :
            available_detectors = list(default_noise_psds[noise_psds[0]].keys())
        else :
            raise ValueError(
                "When `detector_list` is None, `noise_psds` must be exactly one single built-in tag.\n"
                f"Provided: {noise_psds}\n"
                f"Available tags: {list(default_noise_psds.keys())}"
            )
    
    # Get all path maps
    try :
        kv_pairs = np.array(np.broadcast_arrays(available_detectors, noise_psds)).T
    except ValueError as e :
        raise ValueError(
            f"Failed to match noise psds to detectors due to a shape mismatch.\n"
            f"Could not broadcast shapes: psds shape {np.shape(noise_psds)} "
            f"with detector list shape {np.shape(available_detectors)}.\n"
            f"Ensure lengths match or one sequence has a length of 1."
        ) from e

    paths_map = {}
    for detector, noise_psd in kv_pairs :

        detector = str(detector)
        if detector in paths_map :
            raise ValueError(f"Ambiguous configuration : multiple noise psds assigned to detector '{detector}'.")

        if noise_psd is None :
            paths_map[detector] = None
        elif noise_psd in default_noise_psds :
            if detector in default_noise_psds[noise_psd] :
                paths_map[detector] = default_noise_psds[noise_psd][detector]
            else :
                raise ValueError(
                    f"detector '{detector}' does not have a noise psd defined under tag '{noise_psd}'.\n"
                    f"Available detectors for '{noise_psd}': {list(default_noise_psds[noise_psd].keys())}"
                )
        elif Path(noise_psd).is_file() :
            paths_map[detector] = Path(noise_psd)
        else :
            db_summary = "\n".join(
                f"  - '{tag}': {list(detectors.keys())}" for tag, detectors in default_noise_psds.items()
            )
            raise ValueError(
                f"Invalid noise psd: '{noise_psd}'. Must be an existing file path, None, or a valid built-in tag.\n"
                f"Available built-in options (tag: [supported detectors]):\n"
                f"{db_summary}"
            )
    
    return paths_map