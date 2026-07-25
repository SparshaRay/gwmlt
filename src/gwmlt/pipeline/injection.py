"""
Generate GW Injections
"""

from dataclasses import dataclass
import warnings
import numpy as np

from pycbc.types.timeseries import TimeSeries
from lal import LIGOTimeGPS

from gwmlt.core.sources import BinarySystem
from gwmlt.core.morphologies import Morphology, LensedProtocol
from gwmlt.core.observatories import Observatory
from gwmlt.waveform.generator import generate_waveform
from gwmlt.lensing.transform import apply_lensing
from gwmlt.waveform.projector import project_waveform
from gwmlt.noise.inject import inject_noise
from gwmlt.config import config


# Data classes for injection results --------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class DetectorInjection :
    """
    Injection results for a single detector.

    Attributes
    ----------
    noisy_strain : TimeSeries
        The noisy strain time series for the detector.
    clean_strain : TimeSeries
        The clean cropped/padded strain time series for the detector.
    geocenter_td : float
        The signal time delay from the geocenter to the detector.
    snr : float
        The optimal SNR of the injected signal in the detector.
        Will be NaN if no noise was specified.
    _projected_strain : TimeSeries
        The projected strain time series for the detector (before cropping/padding).
        Intended for debugging purposes.
    """
    noisy_strain      : TimeSeries
    clean_strain      : TimeSeries
    geocenter_td      : float
    snr               : float
    _projected_strain : TimeSeries


@dataclass(frozen=True)
class InjectionResult :
    """
    Complete output of the GW injection pipeline.

    Attributes
    ----------
    injections : dict[str, DetectorInjection]
        Dictionary mapping detector names to their corresponding injection results.
    network_snr : float
        The optimal network SNR of the injected signal across all detectors.
        This will be NaN if even one of the detectors have no noise specified.
    lensing_td : float | None
        The time delay introduced by lensing, if applicable. None if no lensing was applied.
    geocent_time : LIGOTimeGPS
        The geocenter trigger time of the injection.
    gen_pars : dict
        Dictionary containing the parameters used for waveform generation.
    _pre_lens_pols : tuple[TimeSeries, TimeSeries]
        The plus and cross polarizations of the waveform before lensing, intended for debugging purposes.
    _post_lens_pols : tuple[TimeSeries, TimeSeries]
        The plus and cross polarizations of the waveform after lensing, intended for debugging purposes.
    """
    injections      : dict[str, DetectorInjection]
    network_snr     : float
    lensing_td      : float | None
    geocent_time    : LIGOTimeGPS
    gen_pars        : dict
    _pre_lens_pols  : tuple[TimeSeries, TimeSeries]
    _post_lens_pols : tuple[TimeSeries, TimeSeries]


# Injection Pipeline ------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def generate_injections(
    source      : BinarySystem,
    morphology  : Morphology,
    observatory : Observatory,
    seed : int = 0,
    **kwargs
) -> InjectionResult :
    
    """
    Full GW injection pipeline for a given source, morphology, and observatory.

    Parameters
    ----------
    source : BinarySystem
        The binary system source parameters.
    morphology : Morphology
        The lensing morphology parameters.
    observatory : Observatory
        The observatory parameters.
    seed : int, optional
        The random seed for noise injection. Default is 0.
    **kwargs
        Additional keyword arguments for waveform generation.

    Returns
    -------
    InjectionResult
        The complete results of the injection pipeline, including detector injections, 
        network SNR, lensing time delay, and waveform generation parameters.
    """

    # Part 0. Sanity checks and warnings

    if hasattr(observatory, 'wf_duration') :

        max_geocent_tds = []
        for detector_str in observatory.active_detectors :
            max_geocent_tds.append(config.injection.max_geocent_to_det_td[detector_str])

        min_req_wf_dur = (
            + max(max_geocent_tds) * 2
            + (
                + config.injection.pre_merger_datapoints
                + config.injection.correlation_grace_order / 2
                + config.injection.highpass_fir_order      / 2
                + config.injection.lowpass_fir_order       / 2
                + config.injection.whiten_filter_order     / 2
                + config.injection.safety_padding          / 2
            ) / observatory.sample_rate
        )

        if int(observatory.wf_duration) < int(min_req_wf_dur) :
            warnings.warn(
                f"Observatory waveform duration ({int(observatory.wf_duration)} s) is less than "
                f"the minimum required duration ({int(min_req_wf_dur)} s) for the current injection configuration.\n"
                "This may lead to unintended truncations. Consider increasing the observatory's `wf_duration`."
            )


    # Part 1. Generate the unlensed waveform in the source frame

    pre_lens_hp, pre_lens_hc, gen_pars = generate_waveform(source, morphology, observatory, **kwargs)


    # Part 2. Apply lensing effects if the morphology is child of lensed protocol

    if isinstance(morphology, LensedProtocol) :
        post_lens_hp, post_lens_hc, lensing_td = apply_lensing(
            pre_lens_hp, pre_lens_hc, morphology.m_lens, morphology.y_lens, morphology.z_lens)
    else :
        post_lens_hp, post_lens_hc = pre_lens_hp.copy(), pre_lens_hc.copy()
        lensing_td = None


    # Part 3. Project the waveforms to get the detector strains

    projected_strains = project_waveform(
        post_lens_hp, post_lens_hc,
        source.ra, source.dec, source.psi,
        observatory
    )


    # Part 4. Inject noise into the projected strains for each detector

    detector_snrs = []
    detector_injections = {}

    for detector_str, (projected_strain, geocenter_td) in projected_strains.items() :

        sample_rate = projected_strain.sample_rate
        max_geocent_to_det_td = config.injection.max_geocent_to_det_td[detector_str]

        injection_start_time = (
            source.geocent_time 
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

        injection_end_time = (
            source.geocent_time
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

        noisy_strain, clean_strain, snr = inject_noise(
            sig_ts          = projected_strain,
            noise_file_path = observatory.resolved_paths[detector_str],
            f_low           = observatory.f_low,
            injection_start_time = injection_start_time,
            injection_end_time   = injection_end_time,
            snr_type = 'optimal',
            seed = int.from_bytes(f"{seed}{detector_str}".encode("utf-8"), "big") % 2147483647,
        )

        detector_snrs.append(snr)
        detector_injections[detector_str] = DetectorInjection(
            noisy_strain      = noisy_strain,
            clean_strain      = clean_strain,
            geocenter_td      = geocenter_td,
            snr               = snr,
            _projected_strain = projected_strain
        )
    

    # Part 5. Compute the network SNR and return the complete injection results

    network_snr = np.linalg.norm(np.array(detector_snrs))

    return InjectionResult(
        injections      = detector_injections,
        network_snr     = network_snr,
        lensing_td      = lensing_td,
        geocent_time    = source.geocent_time,
        gen_pars        = gen_pars,
        _pre_lens_pols  = (pre_lens_hp,  pre_lens_hc),
        _post_lens_pols = (post_lens_hp, post_lens_hc)
    )