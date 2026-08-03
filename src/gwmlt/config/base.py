"""
Resources and Paths
"""

import os
import warnings
from pathlib import Path
from importlib import resources

# Resources
NOISE_DATA_DIR  = resources.files("gwmlt").joinpath("data", "noise")
POPULATION_DIR  = resources.files("gwmlt").joinpath("data", "population")
TEOBRESUMS_FITS = resources.files("gwmlt").joinpath("data", "teobresums_fits") / "teobresums_dimless_evol_fits.npz"

# Database root
try :
    DATABASE_ROOT = Path(os.environ.get("DATABASE_ROOT"))
except TypeError :
    warnings.warn("`DATABASE_ROOT` environment variable is not set.\n"
    "Defaulting to the `database` sibling directory of the GWMLT package root directory.")
    DATABASE_ROOT = Path(__file__).resolve().parents[4] / "database"