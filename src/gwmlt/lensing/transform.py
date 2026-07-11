"""
Apply Lensing Effects on TD Waveforms
"""

import numpy as np
from pycbc.types.timeseries import TimeSeries

from ..core.morphologies import Lensed
from .pointlens import Ff_hybrid, time_delay


def apply_lensing(
    hp: TimeSeries,
    hc: TimeSeries,
    morphology: Lensed
) -> tuple[TimeSeries, TimeSeries] :

    """
    Apply lensing effects on the input TD waveforms.

    Parameters
    ----------
    hp : TimeSeries
        The h_plus timeseries.
    hc : TimeSeries
        The h_cross timeseries.
    morphology : Lensed
        The lensing morphology containing lensing parameters.

    Returns
    -------
    tuple[TimeSeries, TimeSeries]
        A tuple containing the lensed h_plus and h_cross timeseries.
    """

    hp = hp.copy()
    hc = hc.copy()

    ml = morphology.m_lens
    y  = morphology.y_lens
    zl = morphology.z_lens

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

    return lensed_hp, lensed_hc