"""
Inject Timeseries into Noise and calculate PSD
"""

from typing import Literal
from pathlib import Path

import numpy as np

from pycbc.types.timeseries import TimeSeries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.filter import sigma, matched_filter
from lal import LIGOTimeGPS

from ..config import config
from .sampler import sample_gaussian_noise, sample_real_noise


def inject_with_snr(
    sig_ts : TimeSeries,
    noise_file_path : Path | None,
    geocent_time : LIGOTimeGPS,
    geocent_td_pad_for_det : str | None,
    f_low : float,
    seed : int,
    snr_type : Literal['optimal', 'matched'] = 'optimal'
) -> tuple[TimeSeries, TimeSeries, float] :
    
    """
    Inject strain into noise and calculate SNR.

    Parameters
    ----------
    sig_ts : TimeSeries
        The strain timeseries to be injected.
    noise_file_path : Path | None
        Path to the noise PSD or timeseries file. If None, no noise will be added.
    geocent_time : LIGOTimeGPS
        The geocenter trigger time of the signal.
    geocent_td_pad_for_det : str or None
        The name of the detector (e.g., 'H1', 'IndIGO-D', etc.) for which 
        the geocenter to detector time delay padding duration is to be applied.
        If None, no padding duration will be applied.
    f_low : float
        The low frequency cutoff for PSD, noise and SNR calculation.
    seed : int
        The seed for noise sampling.
    snr_type : 'optimal' or 'matched', optional
        The type of SNR to calculate. 'optimal' calculates the optimal SNR (sigma) of the
        injected signal, while 'matched' calculates the matched filter SNR of the padded/cropped
        recovered signal with respect to the injected signal. Default is 'optimal'.
    
    Returns
    -------
    tuple[TimeSeries, TimeSeries, float]
        A tuple containing the injected timeseries, padded/cropped noise-free strain, 
        and the specified type of SNR.
    """

    sample_rate = sig_ts.sample_rate
    max_geocent_to_det_td = config.injection.max_geocent_to_det_td[geocent_td_pad_for_det]

    # Calculate the start and end times for the final timeseries

    start_time = (
        geocent_time 
        - max_geocent_to_det_td
        - (
            + config.injection.pre_merger_datapoints
            + config.injection.correlation_grace_order / 2
            + config.injection.highpass_fir_order      / 2
            + config.injection.lowpass_fir_order       / 2
            + config.injection.whiten_filter_order     / 2
            + config.injection.safety_padding          / 2
        ) / sample_rate
    )

    end_time = (
        geocent_time
        + max_geocent_to_det_td
        + (
            max(
                config.injection.post_merger_datapoints,
                config.injection.correlation_grace_order / 2
            )
            + config.injection.highpass_fir_order        / 2
            + config.injection.lowpass_fir_order         / 2
            + config.injection.whiten_filter_order       / 2
            + config.injection.safety_padding            / 2
        ) / sample_rate
    )


    # If no noise file is provided, just return padded and/or cropped copy of the signal
    if noise_file_path is None :

        noise_ts = TimeSeries(
            np.zeros(int((end_time - start_time) * sample_rate)),
            delta_t = 1.0 / sample_rate,
            epoch = start_time
        )

        injected_ts = noise_ts.inject(sig_ts)
        recovered_ts = (injected_ts - noise_ts).copy()
        snr = np.nan

        return injected_ts, recovered_ts, snr
    
    # If npz file is provided, generate colored Gaussian noise
    elif noise_file_path.suffix == '.npz' :
        noise_ts, psd = sample_gaussian_noise(
            file_path = noise_file_path,
            start_time = start_time,
            end_time = end_time,
            sample_rate = sample_rate,
            f_low = f_low,
            seed = seed
        )
    
    # If hdf5 file is provided, sample real noise timeseries
    elif noise_file_path.suffix == '.hdf5' :
        noise_ts, psd = sample_real_noise(
            file_path = noise_file_path,
            start_time = start_time,
            end_time = end_time,
            sample_rate = sample_rate,
            f_low = f_low,
            seed = seed
        )
    
    # My Disappointment Is Immeasurable And My Day Is Ruined 
    else :
        raise ValueError(
            f"Unsupported noise file format: '{noise_file_path.suffix}'. "
            "Supported formats are .npz and .hdf5."
        )


    # Inject the signal into the noise and extract back the padded/cropped signal
    # Automatically ensures all the timeseries are of the same length and delta_f
    injected_ts = noise_ts.inject(sig_ts) 
    recovered_ts = (injected_ts - noise_ts).copy()

    # Calculate the optimal SNR of the sig_ts
    if snr_type == 'optimal' :
        
        psd = interpolate(psd, 1.0/sig_ts.duration, length=(len(sig_ts)//2+1))
        psd = inverse_spectrum_truncation(
            psd,
            max_filter_len=len(sig_ts),
            low_frequency_cutoff=f_low,
            trunc_method='hann'
        )

        snr = sigma(
            sig_ts,
            psd=psd,
            low_frequency_cutoff=f_low,
            high_frequency_cutoff=sample_rate/2.0
        )

        return injected_ts, recovered_ts, snr
    
    # Calculate the matched filter SNR of the recovered_ts with respect to injected_ts
    elif snr_type == 'matched' :

        # We need to highpass it to get reliable SNR values
        hpass_recovered_ts = recovered_ts.highpass_fir(f_low, int(config.injection.highpass_fir_order/2.0))
        hpass_injected_ts = injected_ts.highpass_fir(f_low, int(config.injection.highpass_fir_order/2.0))

        psd = interpolate(psd, 1.0/hpass_recovered_ts.duration, length=(len(hpass_recovered_ts)//2+1))
        psd = inverse_spectrum_truncation(
            psd, 
            max_filter_len=len(hpass_recovered_ts),
            low_frequency_cutoff=f_low,
            trunc_method='hann'
        )

        snr = max(np.abs(matched_filter(
            template=hpass_recovered_ts,
            data=hpass_injected_ts,
            psd=psd,
            low_frequency_cutoff=f_low,
            high_frequency_cutoff=sample_rate/2.0
        )))

        return injected_ts, recovered_ts, snr

    else :
        raise ValueError(
            f"Unsupported SNR type: '{snr_type}'. "
            "Supported types are 'optimal' and 'matched'."
        )