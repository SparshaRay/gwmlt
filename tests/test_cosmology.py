import pytest

def test_cosmology_z_dl_roundtrip() :
    from gwmlt.utils.cosmology import luminosity_distance, z_from_luminosity_distance

    z_true = 0.5
    dl = luminosity_distance(z_true)
    z_calc = z_from_luminosity_distance(dl)
    
    assert z_calc == pytest.approx(z_true, rel=1e-4)