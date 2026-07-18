"""
Apply Lensing Effects on TD Waveforms
"""

import numpy as np
from plum import dispatch
from pycbc.types.timeseries import TimeSeries

from gwmlt.core.morphologies import LensedProtocol
from gwmlt.core.timeseries import Polarizations
from gwmlt.lensing.pointlens import Ff_hybrid, time_delay


@dispatch
def apply_lensing(
    polarizations : Polarizations,
    lens : LensedProtocol
) -> Polarizations :

    """
    Apply lensing effects on Polarizations state object.

    Parameters
    ----------
    polarizations : Polarizations
        The Polarizations state object containing the h_plus and h_cross timeseries.
    lens : LensedProtocol or child thereof
        The lensing struct containing the lensing parameters.

    Returns
    -------
    Polarizations
        The lensed Polarizations state object.
    """

    lensed_hp, lensed_hc, td = apply_lensing(
        polarizations.hp,
        polarizations.hc,
        lens.m_lens,
        lens.y_lens,
        lens.z_lens,
    )

    return Polarizations(
        hp = lensed_hp,
        hc = lensed_hc,
        geocent_time    = polarizations.geocent_time,
        generation_pars = polarizations.generation_pars,
        lensing_time_delay = td
    )


@dispatch
def apply_lensing(
    hp : TimeSeries,
    hc : TimeSeries,
    ml : float,
    y  : float,
    zl : float = 0.0
) -> tuple[TimeSeries, TimeSeries, float] :

    """
    Apply lensing effects on the input TD waveforms.

    Parameters
    ----------
    hp : TimeSeries
        The h_plus timeseries.
    hc : TimeSeries
        The h_cross timeseries.
    ml : float
        The lens mass in solar masses.
    y : float
        The impact parameter (dimensionless).
    zl : float, optional
        The lens redshift. Default is 0.

    Returns
    -------
    tuple[TimeSeries, TimeSeries, float]
        A tuple containing the lensed h_plus and lensed h_cross timeseries 
        and the lensing time delay.
    """

    hp = hp.copy()
    hc = hc.copy()

    td = time_delay(ml, y, zl)

    td_pad = int(td * hp.sample_rate)
    wf_len = len(hp)
    ev_pad = (wf_len + td_pad) % 2

    hp.append_zeros(td_pad + ev_pad)
    hc.append_zeros(td_pad + ev_pad)

    hp_tilde = hp.to_frequencyseries()
    hc_tilde = hc.to_frequencyseries()

    fs = hp_tilde.sample_frequencies
    Ff = Ff_hybrid(fs, ml, y, zl)

    lensed_hp_tilde = hp_tilde * np.conj(Ff)
    lensed_hc_tilde = hc_tilde * np.conj(Ff)

    lensed_hp = lensed_hp_tilde.to_timeseries()
    lensed_hc = lensed_hc_tilde.to_timeseries()

    return lensed_hp, lensed_hc, td