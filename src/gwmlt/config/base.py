"""
Paths
"""

import os
from pathlib import Path
from importlib import resources

NOISE_DATA_DIR  = resources.files("gwmlt").joinpath("data", "noise")
TEOBRESUMS_FITS = resources.files("gwmlt").joinpath("data", "teobresums_fits") / "teobresums_dimless_evol_fits.npz"

try :
    PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT"))
except TypeError :
    raise TypeError(
        "PROJECT_ROOT environment variable is not set.\n"
        "Please set it to the root directory of your project."
    )