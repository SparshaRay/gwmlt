"""
Paths
"""

import os
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", default=Path(__file__).resolve().parents[3]))
NOISE_DATA_DIR = PROJECT_ROOT / "data" / "noise"