"""
PyCBC TD Waveform Generator
"""

from pycbc.waveform import get_td_waveform

from ..core.sources import BBHSystem
from ..core.profiles import Quasicircular, Eccentric, Lensed
from ..core.observatories import LVK, Decihertz
from ._pars import get_pars


def generate_waveform(
    source: BBHSystem,
    profile: Quasicircular | Eccentric | Lensed,
    observatory: LVK | Decihertz
) :
    """
    Generate time-domain polarization (h_plus, h_cross) waveforms
    """

    pars = get_pars(source, profile, observatory)
    hp, hc = get_td_waveform(**pars)

    hp = hp.taper_timeseries(location='TAPER_STARTEND')
    hc = hc.taper_timeseries(location='TAPER_STARTEND')
    
    return hp, hc