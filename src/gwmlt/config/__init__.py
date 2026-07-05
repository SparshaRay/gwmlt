"""
Global Settings and Configurations
"""

from dataclasses import dataclass, field

from .base import PROJECT_ROOT
from .noise import NoiseConfig
from .waveform import WaveformConfig


@dataclass
class GlobalConfig :
    waveform: WaveformConfig = field(default_factory=WaveformConfig)
    """Configuration for waveform generation"""
    
    noise: NoiseConfig = field(default_factory=NoiseConfig)
    """Configuration for noise data and PSDs"""
    
    project_root = PROJECT_ROOT
    """Root directory of the project
    Determined from the PROJECT_ROOT environment variable.
    If not set, it defaults to the parent root directory of this module.
    """

# Instantiate the global configuration instance
config = GlobalConfig()


__all__ = ["config", "GlobalConfig"]