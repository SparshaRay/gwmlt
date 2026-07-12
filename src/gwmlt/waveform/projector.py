"""
Project GW Polarizations onto Detectors
"""

from warnings import warn

import numpy as np
from fractions import Fraction
from scipy.signal import resample_poly

from pycbc.types.timeseries import TimeSeries
from pycbc.detector import Detector, get_available_detectors
from pycbc.filter import resample_to_delta_t

from ..core.sources import BBHSystem
from ..core.observatories import LVK, Decihertz
from ..decihertz.detectors import IndIGO_D


def project_waveform(
    hp: TimeSeries,
    hc: TimeSeries,
    source: BBHSystem,
    observatory: LVK | Decihertz,
) -> dict[str, tuple[TimeSeries, float]] :

    """
    Project the polarization timeseries onto the observatory detectors.

    Parameters
    ----------
    hp : TimeSeries
        The h_plus polarization timeseries.
    hc : TimeSeries
        The h_cross polarization timeseries.
    source : BBHSystem
        The source parameters of the binary black hole system.
    observatory : LVK | Decihertz
        The observatory parameters.

    Returns
    -------
    dict[str, tuple[TimeSeries, float]]
        A dictionary containing the projected timeseries and their corresponding 
        time delays from geocenter for each detector in the observatory.
    """

    hp = hp.copy()
    hc = hc.copy()

    projected_strains = {}
    pycbc_detectors = get_available_detectors()
    
    for detector_str in observatory.active_detectors :
        
        # Check if the detector is available in PyCBC or is a custom Decihertz detector
        if detector_str in pycbc_detectors :
            det = Detector(detector_str)
        elif detector_str == "IndIGO-D" :
            det = IndIGO_D()
        else :
            raise NotImplementedError(
                f"Detector {detector_str} is not implemented in PyCBC or the custom Decihertz detectors."
            )
        
        # Project the waveform onto the detector
        strain = det.project_wave(hp, hc, source.ra, source.dec, source.psi)
        geocenter_td = det.time_delay_from_earth_center(source.ra, source.dec, float(source.geocent_time))
        
        # If the strain is already at the desired sampling rate, no resampling is needed
        if strain.sample_rate == observatory.sample_rate : projected_strains[detector_str] = (strain, geocenter_td)

        elif abs(strain.sample_rate - observatory.sample_rate) < 1.0 :
            warn(
                f"Waveform sample rate {strain.sample_rate} is very close to target observatory sample rate {observatory.sample_rate}.\n"
                "If not intentional, this usually indicates OBOE either in the waveform generation or post-processing steps."
            )
            recast_strain = TimeSeries(strain, delta_t=1.0/observatory.sample_rate, epoch=strain.start_time)
            projected_strains[detector_str] = (recast_strain, geocenter_td)
        
        elif strain.sample_rate < observatory.sample_rate :
            raise RuntimeError(
                f"Waveform sample rate {strain.sample_rate} is lower than target observatory sample rate {observatory.sample_rate}.\n"
                "Please ensure that the waveform generation sample rate is more than or equal to the target observatory sample rate."
            )

        else :

            # ![NOTE]
            # The following section is intended to downsample decihertz band signals with a sampling rate far below the kHz range.
            # This assumes long wavelength approximation which is not valid above ~333 Hz for IndIGO-D.
            # The signal is projected first and then downsampled because downsampling is not an linear operation.
            # If the waveform itself is generated at the observatory sample rate of ~20Hz, then projection creates aliasing artifacts. 
            # To downsample exclusively by butterworth, the waveform sample rate must be a power-of-two multiple of the observatory sample rate.

            padding_datapoints  = 128 # Put 128 datapoints of padding wrt observatory sample rate
            truncate_datapoints = 96  # Would remove 96 of the padding datapoints after resampling
            zero_pad_seconds    = padding_datapoints / observatory.sample_rate
            strain.append_zeros (int(zero_pad_seconds * strain.sample_rate))
            strain.prepend_zeros(int(zero_pad_seconds * strain.sample_rate))

            try :
                butterworth_strain = resample_to_delta_t(strain, 1.0/observatory.sample_rate, method='butterworth')
                butterworth_strain = butterworth_strain[truncate_datapoints:-truncate_datapoints]
                projected_strains[detector_str] = (butterworth_strain, geocenter_td)

            except Exception as e :
                raise RuntimeError(
                    f"Failed to resample waveform from {strain.sample_rate} Hz to {observatory.sample_rate} Hz for detector {detector_str}.\n"
                    "This usually indicates that the waveform sample rate is not a power-of-two multiple of the target observatory sample rate.\n"
                    "Please ensure that the waveform generation sample rate is set appropriately in the observatory configuration."
                ) from e

    return projected_strains