import numpy as np
from gwmlt.lensing.pointlens import Ff_hybrid, time_delay

def test_unlensed_limit() :
    freqs = np.linspace(10, 512, 100)
    Ff = Ff_hybrid(freqs, ml=0.0, y=1.5)
    np.testing.assert_allclose(np.abs(Ff), 1.0, rtol=1e-3)

def test_lensing_time_delay() :
    td = time_delay(ml=100.0, y=1.0, zl=0.5)
    assert td > 0.0