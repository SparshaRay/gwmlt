import pytest
import numpy as np
from dataclasses import asdict

def test_Lframe_Jframe_equivalence() :
    """
    Test lalsimulation's L-frame and J-frame conversions are working as expected.
    """
    
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