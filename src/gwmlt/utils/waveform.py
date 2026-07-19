"""
Utilities Related to Waveform Generation and Analysis
"""   

from numpy import np
from scipy.signal import find_peaks
from scipy.optimize import root_scalar

from pycbc.types.timeseries import TimeSeries
from pycbc.waveform import get_waveform_filter_length_in_time

from gwmlt.eccentric.analytic import ecc_from_envelop_freqs


# Frequency and phase evolution analysis ----------------------------------------------------------

def pols_to_freq_features(
    hp : TimeSeries,
    hc : TimeSeries,
) -> dict[str, np.ndarray] :
    
    """
    Get various frequency and phase evolution trends over time from the plus and cross 
    polarizations of the waveform (assumes zero inclination and coalescence time of 0.0).
    The phase and frequency evolution, points of frequency peaks and troughs 
    (corresponding to periastron and apastron points for eccentric binaries),
    points of calculating the orbit average frequency, and a function for obtaining
    interpolated eccentricity values is returned. This is primarily useful to work with 
    eccentric waveforms, but can also be repurposed for microlensed or precessing waveforms. 
    Based on : https://arxiv.org/pdf/2302.11257

    Parameters
    ----------
    hp : TimeSeries
        Plus polarization of the waveform.
    hc : TimeSeries
        Cross polarization of the waveform.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing the frequency and phase evolutions over time, 
        as well as the peak and through locations of the frequency evolution.
    """

    hp = hp.time_slice(hp.start_time, 0.0)
    hc = hc.time_slice(hc.start_time, 0.0)

    h = hp + 1j * hc
    time = np.arange(len(h)) * (hp.delta_t)

    phase  = np.unwrap(np.angle(h))
    f_inst = np.gradient(phase, time) / (2.0 * np.pi)

    peaks, _   = find_peaks(+f_inst)
    troughs, _ = find_peaks(-f_inst)

    t_peri = -time[::-1][peaks]
    f_peri = f_inst[peaks]

    phase_peri = phase[peaks]

    t_apos = -time[::-1][troughs]
    f_apos = f_inst[troughs]

    omg_22_dir = []
    omg_22_tms = []

    for i in range(len(peaks)-1) :
        dphase = phase_peri[i+1] - phase_peri[i]
        dt     = t_peri[i+1] - t_peri[i]
        omg_22_dir.append((dphase / dt) / (2 * np.pi))
        omg_22_tms.append(t_peri[i]+dt/2)

    omg_22_dir = np.array(omg_22_dir)
    omg_22_tms = np.array(omg_22_tms)

    peri_fit_poly = np.polynomial.Polynomial.fit(np.log10(-t_peri), np.log10(f_peri), deg=6)
    apos_fit_poly = np.polynomial.Polynomial.fit(np.log10(-t_apos), np.log10(f_apos), deg=6)

    def interp_ecc(t) :
        interp_apos = 10**apos_fit_poly(np.log10(-t))
        interp_peri = 10**peri_fit_poly(np.log10(-t))
        return ecc_from_envelop_freqs(interp_apos, interp_peri)

    return {
        "time"          : -time[::-1],   # Full time array
        "phase"         : phase,         # Phase
        "f_inst"        : f_inst,        # Instantaneous frequency
        "t_peri"        : t_peri,        # Periastron times
        "f_peri"        : f_peri,        # Periastron frequencies
        "peri_fit_poly" : peri_fit_poly, # Periastron frequency fit polynomial
        "phase_peri"    : phase_peri,    # Phase at periastron
        "t_apos"        : t_apos,        # Apastron times
        "f_apos"        : f_apos,        # Apastron frequencies
        "apos_fit_poly" : apos_fit_poly, # Apastron frequency fit polynomial
        "omg_22_tms"    : omg_22_tms,    # omega_22 times
        "omg_22_dir"    : omg_22_dir,    # omega_22 from periastron phase
        "interp_ecc"    : interp_ecc,    # Interpolated eccentricity function
    }


# IMRPhenom family initial conditions -------------------------------------------------------------

def imrphenom_fillen_initconds(
    mass_1 : float, mass_2 : float,
    chi1z  : float, chi2z  : float,
    waveform_duration : float
) -> float :
    
    """
    Get the starting frequency (f_start) for the 
    IMRPhenom family of waveform approximants.

    Parameters
    ----------
    mass_1 : float
        Mass of the primary in solar masses.
    mass_2 : float
        Mass of the secondary in solar masses.
    chi1z : float
        Dimensionless spin of the primary along the z-axis.
    chi2z : float
        Dimensionless spin of the secondary along the z-axis.
    waveform_duration : float
        Duration of the waveform in seconds.
    
    Returns
    -------
    float
        Starting frequency (f_start) in Hz.

    Note
    ----
    This function internally uses IMRPhenomD to estimate the filter length, 
    but it is applicable for IMRPhenomXO4a as well.
    """

    # Calibrate requested waveform duration for IMRPhenomXO4a
    calibrated_waveform_duration = waveform_duration / 0.9

    def waveform_duration_error(f_start) :
        sig_len = get_waveform_filter_length_in_time(
            # Since IMRPhenomXO4a does not support the 
            # `get_waveform_filter_length_in_time` function, 
            # we use IMRPhenomD as a proxy.
            approximant = "IMRPhenomD", 
            mass1   = mass_1, 
            mass2   = mass_2, 
            spin1z  = chi1z, 
            spin2z  = chi2z, 
            f_lower = f_start,
        )
        return sig_len - calibrated_waveform_duration
    
    f_start = root_scalar(
        waveform_duration_error,
        bracket = [1e-2, 256],
        method  = 'brentq',
    ).root

    return f_start