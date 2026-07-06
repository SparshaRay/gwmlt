"""
Structs Defining Observatory Networks and Configurations.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..config import config


# Interferometer Network Configurations -----------------------------------------------------------

@dataclass(frozen=True)
class LVK :

    """
    Ground-based LIGO-Virgo-KAGRA detector network configurations.

    Parameters
    ----------
    noise_profiles : list of str or Path or None
        Target noise profile tag (e.g., [None, 'O4_gaus'] etc.)
    ifo_list : list of str, or None, optional
        Target interferometers (e.g., ['H1', 'L1'] etc.).
        Default is None.
    f_low : float, optional
        Lower bound of frequency range for all analysis (in Hz).
        Sets reference frequency for waveform generation.
        Default is 20.0
    sampling_rate : float, optional
        Sampling rate of the data (in Hz).
        Sets upper bound of frequency range for all analysis at the Nyquist frequency.
        Default is 4096.0

    Attributes
    ----------
    resolved_paths : dict
        A map of detector names to their noise file Path: 
        { 'H1': Path(...), 'L1': Path(...) }

    Note
    ----
    If `ifo_list` is None, it defaults to all interferometers available in the noise profile.
    Otherwise, `noise_profiles` and `ifo_list` must be broadcastable to a common shape.
    """

    noise_profiles : list[str | Path | None]
    ifo_list       : list[str] | None = None
    f_low          : float = 20.0
    sampling_rate  : float = 4096.0

    resolved_paths : dict[str, Path] = field(default_factory=dict, init=False, repr=True)

    def __post_init__(self) :
        paths_map = _resolve_paths(self.noise_profiles, self.ifo_list)
        object.__setattr__(self, "resolved_paths", paths_map)

    @property
    def active_ifos(self) -> list[str] :
        return list(self.resolved_paths.keys())


@dataclass(frozen=True)
class Decihertz :

    """
    Space-based decihertz detector configurations.

    Parameters
    ----------
    wf_duration : float
        Duration of the waveform (in seconds).
    noise_profiles : list of str or Path or None
        Target noise profile tag (e.g., ['S1_gaus'] etc.)
    ifo_list : list of str, or None, optional
        Target interferometers (e.g., ['IndIGO-D'] etc.).
        Default is None.
    ecc_f_ref : float, optional
        Reference frequency for eccentricity.
        Default is 20.0
    f_low : float, optional
        Lower bound of frequency range for all analysis (in Hz).
        Default is 1.0
    sampling_rate : float, optional
        Sampling rate of the data (in Hz).
        Sets upper bound of frequency range for all analysis at the Nyquist frequency.
        Default is 20.0

    Attributes
    ----------
    resolved_paths : dict
        A map of detector names to their noise file Path: 
        { 'IndIGO-D': Path(...) }

    Note
    ----
    If `ifo_list` is None, it defaults to all interferometers available in the noise profile.
    Otherwise, `noise_profiles` and `ifo_list` must be broadcastable to a common shape.
    """

    wf_duration    : float
    noise_profiles : list[str | Path | None]
    ifo_list       : list[str] | None = None
    ecc_f_ref      : float = 20.0
    f_low          : float = 1.00
    sampling_rate  : float = 20.0

    resolved_paths : dict[str, Path] = field(default_factory=dict, init=False, repr=True)

    def __post_init__(self) :
        paths_map = _resolve_paths(self.noise_profiles, self.ifo_list)
        object.__setattr__(self, "resolved_paths", paths_map)

    @property
    def active_ifos(self) -> list[str] :
        return list(self.resolved_paths.keys())


# Helper function to resolve noise profile paths for given interferometers ------------------------

def _resolve_paths(noise_profiles: list[str | Path | None], ifo_list: list[str]) -> dict[str, Path] :

    """
    Resolve the noise profile paths for the given interferometers.

    Parameters
    ----------
    noise_profiles : list of str or Path or None
    ifo_list : list of str, or None

    Returns
    -------
    dict
        A map of interferometer names to their resolved noise file paths.
    """

    # Wrap single values in lists for user convenience, because python does not care about type hints anyways
    # Bad practice, but oh well, saves a second
    noise_profiles = noise_profiles if isinstance(noise_profiles, list) else [noise_profiles]
    available_ifos = ifo_list       if isinstance(ifo_list,       list) else [ifo_list      ]

    default_noise_profiles = config.noise.noise_profiles

    # Handle the case when ifo_list is None, defaults to all ifos available in the noise profile
    if available_ifos == [None] :
        if (len(noise_profiles) == 1) and (noise_profiles[0] in default_noise_profiles) :
            available_ifos = list(default_noise_profiles[noise_profiles[0]].keys())
        else :
            raise ValueError(
                "When `ifo_list` is None, `noise_profiles` must be exactly one single built-in tag.\n"
                f"Provided: {noise_profiles}\n"
                f"Available tags: {list(default_noise_profiles.keys())}"
            )
    
    # Get all path maps
    try :
        kv_pairs = np.array(np.broadcast_arrays(available_ifos, noise_profiles)).T
    except ValueError as e :
        raise ValueError(
            f"Failed to match noise profiles to interferometers due to a shape mismatch.\n"
            f"Could not broadcast shapes: profiles shape {np.shape(noise_profiles)} "
            f"with IFO list shape {np.shape(available_ifos)}.\n"
            f"Ensure lengths match or one sequence has a length of 1."
        ) from e

    paths_map = {}
    for ifo, noise_profile in kv_pairs :

        ifo = str(ifo)
        if ifo in paths_map :
            raise ValueError(f"Ambiguous configuration : multiple noise profiles assigned to interferometer '{ifo}'.")

        if noise_profile is None :
            paths_map[ifo] = None
        elif noise_profile in default_noise_profiles :
            if ifo in default_noise_profiles[noise_profile] :
                paths_map[ifo] = default_noise_profiles[noise_profile][ifo]
            else :
                raise ValueError(
                    f"Interferometer '{ifo}' does not have a profile defined under tag '{noise_profile}'.\n"
                    f"Available interferometers for '{noise_profile}': {list(default_noise_profiles[noise_profile].keys())}"
                )
        elif Path(noise_profile).is_file() :
            paths_map[ifo] = Path(noise_profile)
        else :
            db_summary = "\n".join(
                f"  - '{tag}': {list(ifos.keys())}" for tag, ifos in default_noise_profiles.items()
            )
            raise ValueError(
                f"Invalid noise profile: '{noise_profile}'. Must be an existing file path, None, or a valid built-in tag.\n"
                f"Available built-in options (tag: [supported ifos]):\n"
                f"{db_summary}"
            )
    
    return paths_map