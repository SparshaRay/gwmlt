"""
PyCBC TD Waveform Generator
"""

from warnings import warn

from pycbc.waveform import get_td_waveform
from pycbc.types.timeseries import TimeSeries

from gwmlt.core.sources import BinarySystem
from gwmlt.core.morphologies import Morphology, EccentricProtocol
from gwmlt.core.observatories import Observatory, GroundBased, Decihertz
from gwmlt.utils import imrphenom_initconds, teobresums_initconds
from gwmlt.config import config


def generate_waveform(
    source: BinarySystem,
    morphology: Morphology,
    observatory: Observatory | None,
    **override_kwargs
) -> tuple[TimeSeries, TimeSeries, dict] :
    
    """
    Generate the polarization (h_plus, h_cross) timeseries.

    Parameters
    ----------
    source : BinarySystem
        The compact binary system parameters.
    morphology : Morphology
        The waveform kind.
    observatory : Observatory | None
        The observatory parameters.
    **override_kwargs
        Additional keyword arguments to override any waveform generation parameter.

    Returns
    -------
    tuple[TimeSeries, TimeSeries, dict]
        A tuple containing the generated h_plus timeseries, the h_cross timeseries, 
        and a dictionary containing the parameters passed to the PyCBC `get_td_waveform` function.
    """

    pars = _get_pars(source, morphology, observatory, **override_kwargs)
    hp, hc = get_td_waveform(**pars)

    hp = hp.taper_timeseries(location='TAPER_START')
    hc = hc.taper_timeseries(location='TAPER_START')

    # Pad the waveform with zeros
    # This is not necessary but makes it clearer to see in plots
    hp.append_zeros(64)
    hc.append_zeros(64)
    hp.prepend_zeros(64)
    hc.prepend_zeros(64)

    hp.start_time += source.geocent_time
    hc.start_time += source.geocent_time
    
    # return hp, hc, pars
    return hp, hc, pars


def _get_pars(
    source: BinarySystem,
    morphology: Morphology,
    observatory: Observatory | None,
    **override_kwargs
) -> dict :

    """
    Generate the parameters dictionary for PyCBC TD waveform generation.

    Parameters
    ----------
    source : BinarySystem
        The compact binary system parameters.
    morphology : Morphology
        The waveform kind.
    observatory : Observatory | None
        The observatory parameters.
    **override_kwargs
        Additional keyword arguments to override any waveform generation parameter.

    Returns
    -------
    dict
        A dictionary containing the parameters for PyCBC `get_td_waveform` function.
    """

    if observatory is None :
        warn("Observatory is set to None. Setting f_lower to 20 Hz "
             "and sample_rate to 4096 Hz for waveform generation.\n"
             "You can override this and other parameters by passing them as keyword arguments")

    # Binary Parameters ---------------------------------------------------------------------------

    binary_params = {
        "mass1"  : source.mass_1,
        "mass2"  : source.mass_2,

        "lambda1" : source.lambda_1,
        "lambda2" : source.lambda_2,

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

    # Set the in-plane spins to zero if precession is not allowed
    if morphology.allow_precession == False :
        binary_params.update({
            "spin1x": 0.0, "spin1y": 0.0,
            "spin2x": 0.0, "spin2y": 0.0
        })

    # Initial Conditions for Waveform Generation --------------------------------------------------
    
    initial_conditions = {}

    if observatory is None :
        initial_conditions["f_lower"] = 20.0
        if isinstance(morphology, EccentricProtocol) :
            initial_conditions["ecc"]     = morphology.eccentricity
            initial_conditions["anomaly"] = morphology.anomaly

    elif isinstance(observatory, GroundBased) :
        initial_conditions["f_lower"] = observatory.f_low
        if isinstance(morphology, EccentricProtocol) :
            initial_conditions["ecc"]     = morphology.eccentricity
            initial_conditions["anomaly"] = morphology.anomaly

    elif isinstance(observatory, Decihertz) :

        if morphology.wf_approximant == "IMRPhenomXO4a" :
            f_start = imrphenom_initconds(
                mass_1 = source.mass_1,
                mass_2 = source.mass_2,
                spin1z = source.l_frame.spin1z,
                spin2z = source.l_frame.spin2z,
                waveform_duration = observatory.wf_duration,
            )
            initial_conditions["f_lower"] = f_start
        
        elif morphology.wf_approximant == "teobresums" and isinstance(morphology, EccentricProtocol) :
            ecc_start, f_start = teobresums_initconds(
                mass_1 = source.mass_1,
                mass_2 = source.mass_2,
                spin1z = source.l_frame.spin1z,
                spin2z = source.l_frame.spin2z,
                ecc_ref = morphology.eccentricity,
                f_ref   = observatory.ecc_f_ref,
                waveform_duration = observatory.wf_duration,
            )
            initial_conditions["f_lower"] = f_start
            initial_conditions["ecc"]     = ecc_start
            initial_conditions["anomaly"] = morphology.anomaly
        
        else : raise NotImplementedError(
            f"Waveform approximant '{morphology.wf_approximant}' with morphology "
            f"'{type(morphology).__name__}' is not supported for Decihertz band.")
          
    else : raise NotImplementedError(f"Observatory {observatory} is not supported")

    # Waveform Generation Parameters --------------------------------------------------------------

    sample_rate = observatory.wf_gen_srate if observatory is not None else 4096.0
    waveform_params = {
        "approximant" : morphology.wf_approximant,
        "delta_t"     : 1.0 / sample_rate,
    }

    if morphology.wf_approximant == "IMRPhenomXO4a" :
        waveform_params["mode_array"] = config.waveform.imrphenomx_modes
        if initial_conditions["f_lower"] <= config.waveform.f_ref :
            waveform_params["f_ref"] = config.waveform.f_ref
        else :
            warn(
                f"f_lower ({initial_conditions['f_lower']} Hz) is greater than default f_ref ({config.waveform.f_ref} Hz).\n"
                "Setting f_ref = f_lower for IMRPhenomXO4a waveform generation."
            )
            waveform_params["f_ref"] = initial_conditions["f_lower"]
    
    elif morphology.wf_approximant == "teobresums" :
        # TEOBResumS does not support f_ref, it uses the f_lower as the reference frequency for spin evolution.
        waveform_params["use_mode_lm"] = config.waveform.teobresums_modes
        waveform_params["ecc_freq"]    = config.waveform.ecc_freq
    
    else : warn(f"No special parameters defined for waveform approximant {morphology.wf_approximant}. Using default settings.")

    # Merge all parameters ------------------------------------------------------------------------

    # Unpack and repack
    pars = {**binary_params, **initial_conditions, **waveform_params}
    # Let user override any of the parameters (e.g., changing f_lower, delta_t, etc.)
    pars.update(override_kwargs)
    # Return the final parameters dictionary
    return pars