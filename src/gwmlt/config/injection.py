"""
Injection Duration Configurations
"""

from dataclasses import dataclass, field
from pycbc.detector import get_available_detectors


@dataclass(frozen=True)
class InjectionConfig :
    """
    Configurations for Injections.
    """

    pre_merger_datapoints   : int = int(4096 * 2.5)
    """
    The minimum number of datapoints of data before the merger that must be available after bandpass and whitening.
    Default global value is 10240 (2.5 seconds in 4096 Hz sampling).
    """

    post_merger_datapoints  : int = int(4096 * 1.0)
    """
    The minimum number of datapoints of data after the merger that must be available after bandpass and whitening.
    Default global value is 4096 (1.0 second in 4096 Hz sampling).
    """

    correlation_grace_order : int = 2048 * 2
    """
    Total number of datapoints around the merger that must be available for cross-correlation after bandpass and whitening. 
    This is used to generate SNR time series and to compute cross-correlation maps.
    Must be even integer. Default global value is 4096.
    """

    highpass_fir_order  : int = 512 * 2
    """
    Total number of corrupted datapoints to be removed in highpass filtering. Must be even integer.
    Default global value is 1024.
    """

    lowpass_fir_order   : int = 512 * 2
    """
    Total number of corrupted datapoints to be removed in lowpass filtering. Must be even integer.
    Default global value is 1024.
    """

    whiten_filter_order : int = 1024 * 2
    """
    Total number of corrupted datapoints to be removed after whitening. Must be even integer.
    Default global value is 2048.
    """

    max_geocent_to_det_td : dict[str, float] = field(default_factory=lambda : {
        # All ground based detectors have a maximum geocenter to interferometer time delay of 21.3 ms.
        **{detector: 0.0213 for detector in get_available_detectors()},
        # Heliocentric IndIGO-D has a maximum geocenter to interferometer time delay of around 175 seconds.
        **{'IndIGO-D' : 175},
        # If None then no padding will be applied.
        **{None : 0.0},
    })
    """
    Dictionary containing maximum possible time delays from geocenter to each interferometer (in seconds).
    """

    safety_padding : int = 4 * 2
    """
    Total number of extra datapoints to be added to the start and end of the final timeseries.
    This is to accommodate for any precision or rounding or off-by-one errors in the downstream pipeline.
    Usually not needed, but its a good idea to have a few extra datapoints just in case. Must be even integer.
    """