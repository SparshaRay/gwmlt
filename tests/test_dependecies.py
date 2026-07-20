import pytest
import numpy as np


def test_julia_hypergeometric_functions() :
    try:
        from gwmlt.lensing.pointlens import hyp1f1
    except Exception as e:
        pytest.fail(f"Failed to import Julia hypergeometric binding: {e}")

    # $_1F_1(1, 1, z) = e^z$ identity test
    z_vals = np.array([0.0, 1.0, 2.5, -1.0])
    a_vec = np.ones_like(z_vals)
    b_vec = np.ones_like(z_vals)

    try:
        results = hyp1f1(a_vec, b_vec, z_vals)
        expected = np.exp(z_vals)
        results_np = np.array(results)

        np.testing.assert_allclose(
            results_np, 
            expected, 
            rtol=1e-7, 
            err_msg="HypergeometricFunctions.jl output does not match expected values."
        )
    except Exception as e:
        pytest.fail(f"Julia HypergeometricFunctions evaluation failed: {e}")


def test_teobresums_dali_availability() :
    from pycbc.waveform import td_approximants, get_td_waveform
    assert 'teobresums' in td_approximants(), 'TEOBResumS (Dali) approximant is not available in PyCBC time-domain approximants.'
    try : 
        get_td_waveform(approximant="teobresums", mass1=30, mass2=30, ecc=0.5, f_lower=20, lambda1=0, lambda2=0, delta_t=1.0/4096)
    except Exception as e :
        pytest.fail(f"TEOBResumS (Dali) approximant is not functioning correctly: {e}")