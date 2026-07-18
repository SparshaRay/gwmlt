"""
State Objects for GW Signals
"""

from dataclasses import dataclass
from pycbc.types.timeseries import TimeSeries


@dataclass(frozen=True)
class Polarizations :
    """
    Wrapper for generated gravitational wave polarizations,
    its geocenter trigger time and the generating parameters.

    Attributes
    ----------
    hp : TimeSeries
        The h_plus timeseries.
    hc : TimeSeries
        The h_cross timeseries.
    geocent_time : float
        The geocentric trigger time (coalescence time at the center of the Earth).
    generation_pars : dict
        A dictionary containing the parameters which were
        passed to the  PyCBC `get_td_waveform` function.
    lensing_time_delay : float | None
        The time delay due to lensing, if lensed. None if not lensed.
    """
    hp: TimeSeries
    hc: TimeSeries
    geocent_time: float
    generation_pars: dict
    lensing_time_delay: float | None = None