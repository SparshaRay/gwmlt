"""
PSD and Real Noise File Paths
"""

from dataclasses import dataclass, field
from pathlib import Path

from gwmlt.config.base import NOISE_DATA_DIR

@dataclass(frozen=True)
class NoiseConfig :
    """
    Configurations for noise data and PSDs.
    """

    data_dir : Path = NOISE_DATA_DIR
    """Root directory holding noise files and PSDs."""

    noise_psds : dict[str, dict[str, Path]] = field(default_factory=lambda : {

        # Real glitch-free event-free detector noise timeseries 
        "O4_real": {
            "H1": NOISE_DATA_DIR / "real_noise" / "H1_noise.hdf5",
            "L1": NOISE_DATA_DIR / "real_noise" / "L1_noise.hdf5",
        },

        # O4 High
        "O4_gaus": {
            "H1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O4_high_psd.npz",
            "L1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O4_high_psd.npz",
            "V1": NOISE_DATA_DIR / "detector_PSDs" / "aVirgo_O4_high_psd.npz",
        },

        # O5b
        "O5b_gaus": {
            "H1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O5b_psd.npz",
            "L1": NOISE_DATA_DIR / "detector_PSDs" / "aLIGO_O5b_psd.npz",
            "V1": NOISE_DATA_DIR / "detector_PSDs" / "aVirgo_O5_high_psd.npz",
        },

        # IndIGO-D S2
        "S2_gaus": {
            "IndIGO-D": NOISE_DATA_DIR / "detector_PSDs" / "IndIGO-D_S2_psd.npz",
        },

    })
    """Dictionary mapping detector network setups to their respective PSD or timeseries paths."""