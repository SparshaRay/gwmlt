import pytest
import numpy as np
from pycbc.types.timeseries import TimeSeries

from gwmlt.core.sources import BBHSystem, LFrameSpins
from gwmlt.core.morphologies import Quasicircular, Eccentric, Lensed
from gwmlt.core.observatories import GroundBased, Decihertz
from gwmlt.pipeline.injection import generate_injections


def test_observatory_broadcasting() :
    
    # 1. Single PSD, Multiple Detectors
    obs1 = GroundBased(
        noise_psds=['O4_gaus'], 
        detector_list=['H1', 'L1']
    )
    assert 'H1' in obs1.active_detectors
    assert 'L1' in obs1.active_detectors
    assert obs1.resolved_paths['H1'].name == "aLIGO_O4_high_psd.npz"

    # 2. Multiple PSDs, Multiple Detectors
    obs2 = GroundBased(
        noise_psds=['O4_gaus', 'O5b_gaus'], 
        detector_list=['H1', 'V1']
    )
    assert 'H1' in obs2.active_detectors
    assert 'V1' in obs2.active_detectors
    # V1 should be assigned the O5b PSD
    assert obs2.resolved_paths['V1'].name == "aVirgo_O5_high_psd.npz"

    # 3. Default inference (detector_list=None)
    obs3 = GroundBased(
        noise_psds=['O4_real'],
        detector_list=None
    )
    assert set(obs3.active_detectors) == {'H1', 'L1'}
    assert obs3.resolved_paths['H1'].name == "H1_noise.hdf5"


def test_observatory_invalid_configurations() :
    
    # Shape mismatch (2 PSDs for 3 detectors)
    with pytest.raises(ValueError, match="Failed to match noise psds to detectors"):
        GroundBased(
            noise_psds=['O4_gaus', 'O5b_gaus'], 
            detector_list=['H1', 'L1', 'V1']
        )
        
    # Invalid built-in tag
    with pytest.raises(ValueError, match="Invalid noise psd: 'NonExistentTag'"):
        GroundBased(
            noise_psds=['NonExistentTag'], 
            detector_list=['H1']
        )


@pytest.fixture
def base_source():
    return BBHSystem(
        mass_1=45.0,
        mass_2=30.0,
        luminosity_distance=500.0,
        spins=LFrameSpins(spin1z=0.1, spin2z=-0.5)
    )

@pytest.mark.parametrize(
    "morphology, observatory, expected_detectors",
    [
        # 1. GroundBased + Quasicircular (Standard LIGO/Virgo run with O4 real noise)
        (
            Quasicircular(),
            GroundBased(noise_psds=['O4_real'], detector_list=['H1', 'L1']),
            ['H1', 'L1']
        ),
        # 2. GroundBased + Lensed (Microlensing test with O5b Gaussian noise)
        (
            Lensed(m_lens=200.0, y_lens=1.5),
            GroundBased(noise_psds=['O5b_gaus'], detector_list=['V1']),
            ['V1']
        ),
        # 4. Decihertz + Eccentric (TEOBResumS with IndIGO-D S2 noise)
        (
            Eccentric(eccentricity=0.075),
            Decihertz(wf_duration=1200.0, noise_psds=['S2_gaus'], detector_list=['IndIGO-D']),
            ['IndIGO-D']
        ),
    ]
)
def test_full_injection_pipeline(base_source, morphology, observatory, expected_detectors):

    result = generate_injections(
        source=base_source, 
        morphology=morphology, 
        observatory=observatory, 
        seed=42
    )

    # Validate active detectors
    assert list(result.injections.keys()) == expected_detectors

    # Validate the Network SNR 
    assert not np.isnan(result.network_snr)
    assert result.network_snr > 0.0

    # Validate Lensing outputs
    if isinstance(morphology, Lensed):
        assert result.lensing_td is not None
        assert result.lensing_td > 0.0
    else:
        assert result.lensing_td is None

    # Validate individual detector injections
    for det in expected_detectors:
        inj = result.injections[det]
        
        # Ensure timeseries objects are generated
        assert isinstance(inj.noisy_strain, TimeSeries)
        assert isinstance(inj.clean_strain, TimeSeries)
        
        # Ensure lengths match
        assert len(inj.noisy_strain) == len(inj.clean_strain)
        
        # Verify valid time delay evaluations
        assert isinstance(inj.geocenter_td, float)
        
        # Verify that actual noise was added (noisy should NOT exactly equal clean)
        assert not np.array_equal(inj.noisy_strain.data, inj.clean_strain.data)
        
        # Validate that the detector SNR is calculated successfully
        assert not np.isnan(inj.snr)
        assert inj.snr > 0.0

        # Length check for Decihertz
        if isinstance(observatory, Decihertz):
            # Check if the generated (pre-padding) waveform duration is within 5% of requested wf_duration
            actual_wf_duration = inj._projected_strain.duration
            assert actual_wf_duration == pytest.approx(observatory.wf_duration, rel=0.05), \
                f"Generated waveform duration {actual_wf_duration}s is not within 5% of requested {observatory.wf_duration}s."