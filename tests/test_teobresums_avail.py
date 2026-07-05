import pytest
from pycbc.waveform import td_approximants, get_td_waveform

def test_teobresums_dali_availability() :
    """
    Verify that the TEOBResumS (Dalí) approximant is available via PyCBC.
    """
    
    assert 'teobresums' in td_approximants(), 'TEOBResumS (Dali) approximant is not available in PyCBC time-domain approximants.'
    try : 
        get_td_waveform(approximant="teobresums", mass1=30, mass2=30, ecc=0.5, f_lower=20, lambda1=0, lambda2=0, delta_t=1.0/4096)
    except Exception as e :
        pytest.fail(f"TEOBResumS (Dali) approximant is not functioning correctly: {e}")