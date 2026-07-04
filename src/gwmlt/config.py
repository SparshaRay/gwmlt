"""
Configurations for the modules. 
"""

import os
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", default=Path(__file__).resolve().parents[2]))
NOISE_DATA_DIR = PROJECT_ROOT / "data" / "noise"

DEFAULT_NOISE_PROFILES = {

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

    # IndIGO-D High
    "S1_gaus": {
        "IndIGO-D": NOISE_DATA_DIR / "detector_PSDs" / "IndIGO-D_S2_psd.npz",
    },

}