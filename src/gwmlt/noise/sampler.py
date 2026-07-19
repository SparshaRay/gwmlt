"""
Generate Noise Time Series and Corresponding PSDs 
"""

import h5py
import numpy as np
import pandas as pd
from pathlib import Path

from pycbc.types.timeseries import TimeSeries
from pycbc.types.frequencyseries import FrequencySeries
from pycbc.psd import read, interpolate, inverse_spectrum_truncation
from pycbc.noise.reproduceable import colored_noise
from lal import LIGOTimeGPS


def sample_gaussian_noise(
    file_path: Path, 
    start_time: LIGOTimeGPS,
    end_time: LIGOTimeGPS, 
    sample_rate: float, 
    f_low: float,
    seed: int,
    psd_len: int = 65537
) -> tuple[TimeSeries, FrequencySeries] :
    
    """
    Sample colored Gaussian noise time series from a given PSD.

    Parameters
    ----------
    file_path : Path
        Path to the PSD file (in .npz format, containing 'freq' and 'psd' arrays).
    start_time : LIGOTimeGPS
        Start time of the noise time series.
    end_time : LIGOTimeGPS
        End time of the noise time series.
    sample_rate : float
        Sampling rate of the noise time series.
    f_low : float
        Lower frequency cutoff for the noise and PSD.
    seed : int
        Seed for reproducibility.
    psd_len : int, optional
        Length of the PSD to be returned. Default is 65537.
    
    Returns
    -------
    tuple[TimeSeries, FrequencySeries]
        A tuple containing the time series of the sampled noise and the corresponding PSD.
    """

    with np.load(file_path) as f:
        freq = f['freq']
        psd  = f['psd']

    psd_duration = (psd_len - 1) / (sample_rate / 2)
    delta_f = 1.0 / psd_duration

    length = int(sample_rate / (2 * delta_f)) + 1

    psd = read.from_numpy_arrays(
        freq_data=freq,
        noise_data=psd,
        length=length,
        delta_f=delta_f,
        low_freq_cutoff=f_low
    )

    # Flat out below f_low to avoid interpolation artifacts downstream
    cutoff_idx = (psd!=0.0).argmax()
    psd[:cutoff_idx] = psd[cutoff_idx]

    # We have to generate a longer time segment and crop it down because PyCBC applies
    # an integer cast in colored_noise that sometimes throws off-by-one errors.
    ts = colored_noise(
        psd,
        start_time=np.floor(start_time),
        end_time=np.ceil(end_time),
        seed=seed,
        sample_rate=sample_rate,
        low_frequency_cutoff=f_low,
        filter_duration=4.0
    ).time_slice(start_time, end_time)

    return ts, psd


def sample_real_noise(
    file_path: Path, 
    start_time: LIGOTimeGPS,
    end_time: LIGOTimeGPS, 
    sample_rate: float, 
    f_low: float,
    seed: int,
    psd_len: int = 65537
) -> tuple[TimeSeries, FrequencySeries] :
    
    """
    Sample a segment of real noise from a HDF5 file containing timeseries data.

    Parameters
    ----------
    file_path : Path
        Path to the HDF5 file containing the noise timeseries data.
    start_time : LIGOTimeGPS
        The start time of the desired noise timeseries.
    end_time : LIGOTimeGPS
        The end time of the desired noise timeseries.
    sample_rate : float
        The sample rate of the noise data.
    f_low : float
        The low frequency cutoff for the PSD.
    seed : int
        Seed for reproducibility.
    psd_len : int, optional
        Length of the PSD to be generated. Default is 65537.

    Returns
    -------
    tuple[TimeSeries, FrequencySeries]
        A tuple containing the sampled noise TimeSeries and its corresponding Power Spectral Density (PSD).
    """

    index_df = pd.read_hdf(file_path, key='index')
    ts_len = int(float(end_time - start_time) * sample_rate)

    with h5py.File(file_path, 'r') as f :

        ts_group = f['timeseries']

        if sample_rate != ts_group.attrs['sample_rate'] : raise ValueError(
            f"Requested sample rate {sample_rate} does not match noise data sample rate {ts_group.attrs['sample_rate']}."
        )

        df_valid = index_df[index_df['length'] >= ts_len].copy()

        if df_valid.empty : raise ValueError(
            f"No timeseries segments exceed the requested length threshold of {ts_len}"
        )

        df_valid['sampling_weight'] = df_valid['length'] - ts_len + 1
        selected_row = df_valid.sample(n=1, weights='sampling_weight', random_state=seed).iloc[0]

        ts_id, total_length = selected_row['ts_id'], selected_row['length']

        max_start_index = total_length - ts_len
        start_index = np.random.default_rng(seed).integers(0, max_start_index+1)

        selected_ts = ts_group[f'ts_{ts_id:04d}'][start_index:start_index+ts_len]
    
    ts = TimeSeries(selected_ts, delta_t=1.0/sample_rate, epoch=start_time)

    psd = ts.psd(segment_duration=ts.duration/5) # Average PSD over 5*2+1 = 11 segments
    psd = interpolate(psd, delta_f=(sample_rate/2)/(psd_len-1), length=psd_len)
    psd = inverse_spectrum_truncation(
        psd, 
        max_filter_len=int(sample_rate*4), # Time domain filter is 4 seconds long
        low_frequency_cutoff=f_low,
        trunc_method='hann'
    )

    return ts, psd