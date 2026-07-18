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
    """
    hp: TimeSeries
    hc: TimeSeries
    geocent_time: float
    generation_pars: dict