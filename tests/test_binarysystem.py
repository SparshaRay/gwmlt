import pytest
import numpy as np
from dataclasses import asdict


def test_Lframe_Jframe_equivalence() :
    from gwmlt.core.sources import BBHSystem, JFrameSpins, LFrameSpins

    jframe_init = BBHSystem(mass_1=10, mass_2=10, spins=JFrameSpins(a_1=0.5, a_2=0.75, tilt_1=0.0, tilt_2=np.pi))
    lframe_init = BBHSystem(mass_1=10, mass_2=10, spins=LFrameSpins(inclination=0.0, spin1z=0.5, spin2z=-0.75))

    l_frame_eqiv = asdict(jframe_init.l_frame) == pytest.approx(asdict(lframe_init.l_frame))
    j_frame_eqiv = asdict(jframe_init.j_frame) == pytest.approx(asdict(lframe_init.j_frame))

    if not l_frame_eqiv or not j_frame_eqiv :
        pytest.fail(f"L-frame and J-frame conversions are not equivalent :\n"
                    f"J-frame (as initialized):       {asdict(jframe_init.j_frame)}\n"
                    f"J-frame (derived from L-frame): {asdict(lframe_init.j_frame)}\n"
                    f"L-frame (as initialized):       {asdict(lframe_init.l_frame)}\n"
                    f"L-frame (derived from J-frame): {asdict(jframe_init.l_frame)}")


def test_binary_system_factory_roundtrip() :
    from gwmlt.core.sources import BBHSystem
    sys1 = BBHSystem(mass_1=30.0, mass_2=20.0)
    sys2 = BBHSystem.from_mtot(mtot=sys1.mtot, q=1.0/sys1.q)
    assert sys1.mass_1 == pytest.approx(sys2.mass_2)
    assert sys1.mass_2 == pytest.approx(sys2.mass_1)
    assert sys1.mchirp == pytest.approx(sys2.mchirp)


def test_chi_p_mass_ordering_invariance() :
    from gwmlt.core.sources import LFrameSpins, BBHSystem
    spins1 = LFrameSpins(spin1x=0.4, spin1y=0.3, spin2x=0.1, spin2y=0.2)
    spins2 = LFrameSpins(spin2x=0.4, spin2y=0.3, spin1x=0.1, spin1y=0.2)
    sys1 = BBHSystem(mass_1=40.0, mass_2=20.0, spins=spins1)
    sys2 = BBHSystem(mass_1=20.0, mass_2=40.0, spins=spins2)
    assert sys1.chi_p == pytest.approx(sys2.chi_p)