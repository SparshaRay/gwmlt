"""
PyCBC TD Waveform Generator
"""

from warnings import warn

from pycbc.waveform import get_td_waveform
from pycbc.types.timeseries import TimeSeries

from ..core.sources import BBHSystem
from ..core.morphologies import Quasicircular, Eccentric, Lensed
from ..core.observatories import LVK, Decihertz
from ..config import config


def generate_waveform(
    source: BBHSystem,
    morphology: Quasicircular | Eccentric | Lensed,
    observatory: LVK | Decihertz | None,
    **user_override_kwargs
) -> tuple[TimeSeries, TimeSeries, dict] :
    
    """
    Generate polarization (h_plus, h_cross) timeseries

    Parameters
    ----------
    source : BBHSystem
        The binary black hole system parameters.
    morphology : Quasicircular | Eccentric | Lensed
        The waveform kind.
    observatory : LVK | Decihertz | None
        The observatory parameters.
    **user_override_kwargs
        Additional keyword arguments to override any waveform generation parameter.

    Returns
    -------
    tuple[TimeSeries, TimeSeries, dict]
        A tuple containing the generated h_plus timeseries, the h_cross timeseries,
        and a dictionary containing the parameters for PyCBC `get_td_waveform` function.
    """

    pars = _get_pars(source, morphology, observatory, **user_override_kwargs)
    hp, hc = get_td_waveform(**pars)

    hp = hp.taper_timeseries(location='TAPER_STARTEND')
    hc = hc.taper_timeseries(location='TAPER_STARTEND')
    
    return hp, hc, pars


def _get_pars(
    source: BBHSystem,
    morphology: Quasicircular | Eccentric | Lensed,
    observatory: LVK | Decihertz | None,
    **user_override_kwargs
) -> dict :

    """
    Generate the parameters dictionary for PyCBC TD waveform generation.

    Parameters
    ----------
    source : BBHSystem
        The binary black hole system parameters.
    morphology : Quasicircular | Eccentric | Lensed
        The waveform kind.
    observatory : LVK | Decihertz | None
        The observatory parameters.
    **user_override_kwargs
        Additional keyword arguments to override any waveform generation parameter.

    Returns
    -------
    dict
        A dictionary containing the parameters for PyCBC `get_td_waveform` function.
    """


    # Binary Parameters ---------------------------------------------------------------------------

    binary_params = {
        "mass1"  : source.mass_1,
        "mass2"  : source.mass_2,

        "spin1x" : source.l_frame.spin1x,
        "spin1y" : source.l_frame.spin1y,
        "spin1z" : source.l_frame.spin1z,

        "spin2x" : source.l_frame.spin2x,
        "spin2y" : source.l_frame.spin2y,
        "spin2z" : source.l_frame.spin2z,

        "inclination" : source.l_frame.inclination,
        "coa_phase"   : source.phase,
        "distance"    : source.luminosity_distance,
    }

    # Since teobresums expects lambda1 and lambda2 values to be provided as well
    if morphology.wf_approximant == "teobresums" :
        # No tidal effects for black holes
        binary_params["lambda1"] = 0.0
        binary_params["lambda2"] = 0.0


    # Initial Conditions for Waveform Generation --------------------------------------------------
    
    initial_conditions = {}

    if observatory is None :
        warn('Observatory is set to None. Setting f_lower to 20 Hz for waveform generation. '
             'You can override this and other parameters by passing them as keyword arguments')
        initial_conditions["f_lower"] = 20.0
        if morphology.wf_approximant == "teobresums" :
            initial_conditions["ecc"]     = morphology.eccentricity
            initial_conditions["anomaly"] = morphology.anomaly

    elif isinstance(observatory, LVK) :

        if morphology.wf_approximant == "IMRPhenomXO4a" :
            initial_conditions["f_lower"] = observatory.f_low

        elif morphology.wf_approximant == "teobresums" :
            initial_conditions["f_lower"] = observatory.f_low
            initial_conditions["ecc"]     = morphology.eccentricity
            initial_conditions["anomaly"] = morphology.anomaly
        
        else : raise NotImplementedError(f"Waveform approximant {morphology.wf_approximant} is not supported")

    elif isinstance(observatory, Decihertz) :

        raise NotImplementedError("Decihertz band initial conditions are not implemented yet")

        # if morphology.wf_approximant == "IMRPhenomXO4a" :
        #     f_start = get_imrphenomx_fstart(source, observatory.wf_duration)
        #     initial_conditions["f_lower"] = f_start
        
        # elif morphology.wf_approximant == "teobresums" :
        #     f_start, ecc_start = get_teobresums_fstart(
        #         source, observatory.wf_duration, morphology.eccentricity, observatory.ecc_f_ref
        #     )
        #     initial_conditions["f_lower"] = f_start
        #     initial_conditions["ecc"]     = ecc_start
        #     initial_conditions["anomaly"] = morphology.anomaly
        
        # else : raise NotImplementedError(f"Waveform approximant {morphology.wf_approximant} is not supported")
    
    else : raise NotImplementedError(f"Observatory {observatory} is not supported")


    # Waveform Generation Parameters --------------------------------------------------------------

    waveform_params = {
        "approximant" : morphology.wf_approximant,
        "delta_t"     : 1.0 / config.waveform.td_wf_gen_srate,
    }

    if morphology.wf_approximant == "IMRPhenomXO4a" :
        waveform_params["mode_array"] = config.waveform.imrphenomx_modes
        if initial_conditions["f_lower"] <= config.waveform.f_ref :
            waveform_params["f_ref"] = config.waveform.f_ref
        else :
            warn(
                f"f_lower ({initial_conditions['f_lower']} Hz) is greater than default f_ref ({config.waveform.f_ref} Hz). "
                "Setting f_ref = f_lower for IMRPhenomXO4a waveform generation."
            )
            waveform_params["f_ref"] = initial_conditions["f_lower"]
    
    if morphology.wf_approximant == "teobresums" :
        waveform_params["use_mode_lm"] = config.waveform.teobresums_modes
        waveform_params["ecc_freq"]    = config.waveform.ecc_freq


    # Merge all parameters ------------------------------------------------------------------------

    pars = {**binary_params, **initial_conditions, **waveform_params}
    pars.update(morphology.override_pars)
    pars.update(user_override_kwargs)
    return pars