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
    
    for detector_str in observatory.active_ifos :
        
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
        geocenter_td = det.time_delay_from_earth_center(source.ra, source.dec, source.geocent_time)

        # ![WARNING] : This assumes long wavelength approximation which is not valid above ~333 Hz for IndIGO-D.
        # The signal is first projected and then downsampled because the downsampling is not an linear operation.
        # If the waveform itself is generated at IndIGO-D band, then projection creates aliasing artifacts.
        # Usually, for 4096Hz to 20Hz downsampling, the waveform is projected down to 32Hz by Butterworth first,
        # followed by polyphase resampling down to 20Hz. This is done to minimize aliasing artifacts.
        
        # If the strain is already at the desired sampling rate, no resampling is needed
        if strain.sample_rate == observatory.sample_rate : projected_strains[detector_str] = (strain, geocenter_td)
        else : # Resample the projected strain to the desired sampling rate

            if abs(strain.sample_rate - observatory.sample_rate) < 1.0 :
                warn(
                    f"Strain sample rate {strain.sample_rate} is very close to observatory sample rate {observatory.sample_rate}.\n"
                    "If not intentional, this usually indicates OBOE either in the waveform generation or post-processing steps."
                )

            sig_start_time = strain.start_time

            # Padding
            padding_datapoints  = 256 # Do 256 datapoints of padding wrt observatory sample rate
            truncate_datapoints = 192 # Would remove 192 datapoints after resampling
            zero_pad_seconds    = padding_datapoints / observatory.sample_rate
            strain.append_zeros (int(zero_pad_seconds * strain.sample_rate))
            strain.prepend_zeros(int(zero_pad_seconds * strain.sample_rate))

            # Do as much with Butterworth as possible
            resample_ratio    = strain.sample_rate / observatory.sample_rate
            butterworth_ratio = max(1, int(2 ** np.floor(np.log2(resample_ratio))))

            strain_butterworth = resample_to_delta_t(
                strain, 
                1.0/(strain.sample_rate/butterworth_ratio), 
                method='butterworth'
            )

            # Polyphase resampling
            polyphase_ratio = Fraction(
                observatory.sample_rate / strain_butterworth.sample_rate
            )
            up   = polyphase_ratio.numerator
            down = polyphase_ratio.denominator

            if (up>8) and (polyphase_ratio<=1) :
                # Try ldas downsampling if the scale up factor for polyphase is too high 
                warn(
                    f"Polyphase upsampling factor {up} is too high. Switching to ldas.\n"
                    f"Consider changing internal waveform generation sample rate."
                )
                ldas_strain = resample_to_delta_t(
                    strain_butterworth, 
                    1.0/observatory.sample_rate, 
                    method='ldas'
                )
                ldas_strain = ldas_strain[truncate_datapoints:-truncate_datapoints]
                ldas_strain.corrupted_samples = 0
                projected_strains[detector_str] = (ldas_strain, geocenter_td)
                    
            else :
                # Otherwise, use polyphase resampling
                polyphase_strain = resample_poly(strain_butterworth, up, down)
                polyphase_strain = polyphase_strain[truncate_datapoints:-truncate_datapoints]
                polyphase_strain = TimeSeries(
                    polyphase_strain, 
                    delta_t=1.0/observatory.sample_rate, 
                    epoch=sig_start_time
                )
                projected_strains[detector_str] = (polyphase_strain, geocenter_td)

    return projected_strains