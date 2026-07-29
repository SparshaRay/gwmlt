"""
Global Settings and Configurations
"""

from dataclasses import dataclass, field
from contextlib import contextmanager
from typing import Any, Generator
from pathlib import Path

from gwmlt.config.base import PROJECT_ROOT, TEOBRESUMS_FITS
from gwmlt.config.noise import NoiseConfig
from gwmlt.config.waveform import WaveformConfig
from gwmlt.config.injection import InjectionConfig
from gwmlt.config.cosmology import CosmologyConfig
from gwmlt.config.population import PopulationConfig


@dataclass(frozen=True)
class GlobalConfig :
    """
    Struct for global configuration settings for the GWMLT library.
    """

    waveform: WaveformConfig = field(default_factory=WaveformConfig)
    """Configuration for waveform generation"""
    
    noise: NoiseConfig = field(default_factory=NoiseConfig)
    """Configuration for noise data and PSDs"""

    injection: InjectionConfig = field(default_factory=InjectionConfig)
    """Configuration for injection timeseries"""

    cosmology: CosmologyConfig = field(default_factory=CosmologyConfig)
    """Configuration for cosmology constants"""

    population: PopulationConfig = field(default_factory=PopulationConfig)
    """Configuration for population sampling"""

    teobresums_fits_path: Path = TEOBRESUMS_FITS
    """Path to the TEOBResumS evolution fits file."""
    
    project_root: Path = PROJECT_ROOT
    """
    Root directory of the project
    Determined from the PROJECT_ROOT environment variable.
    If not set, it defaults to the parent root directory of this module.
    """

# Instantiate the global configuration instance
config = GlobalConfig()
"""Global configuration instance for the GWMLT library."""


# Context manager to temporarily override global configuration variables
@contextmanager
def config_override(overrides: dict[str, Any]) -> Generator[GlobalConfig, None, None] :

    """
    Temporarily override global configuration variables.
    Accepts a dictionary of keys as strings with the desired values.

    Example usage:
    ```python
    print(config.waveform.f_ref)  # Prints: 20.0
    with config_override({
        "waveform.f_ref": 40.0,
        "project_root": Path("/new/root")
    }):
        print(config.waveform.f_ref)  # Prints: 40.0
    ```    
    """

    old_values = []
    
    try :
        for dot_path, new_value in overrides.items() :
            parts = dot_path.split(".")
            
            target_obj = config
            for part in parts[:-1] : target_obj = getattr(target_obj, part)
                
            attr_name = parts[-1]
            if not hasattr(target_obj, attr_name) :
                raise AttributeError(f"GlobalConfig has no attribute '{dot_path}'")
                
            old_values.append((target_obj, attr_name, getattr(target_obj, attr_name)))
            object.__setattr__(target_obj, attr_name, new_value)
            
        yield config
        
    finally :
        for target_obj, attr_name, old_value in reversed(old_values) :
            object.__setattr__(target_obj, attr_name, old_value)


# Expose
__all__ = [
    "config", 
    "GlobalConfig", 
    "config_override"
]