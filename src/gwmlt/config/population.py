"""
Population Sampling Configurations 
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np

from gwmlt.config.base import POPULATION_DIR


@dataclass(frozen=True)
class PopulationConfig :
    """
    Configurations for population sampling.
    """

    population_file_path : Path = POPULATION_DIR / "GWTC3_pop.json"
    """
    Path to the population file containing the MAP parameters for the population model.
    Default file is "data/population/GWTC3_pop.json", which can be obtained as the
    `o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json` under `analyses/PowerLawPeak` of 
    `GWTC-3-population-data.tar.gz` (https://zenodo.org/records/5655785/files/GWTC-3-population-data.tar.gz)
    """

    m1_src_range : tuple[float, float] = (5.0, 125.0)
    """
    Low and high limits of the primary mass in the source frame.
    Defaults to (5.0, 125.0).
    """

    m_det_cutoff_range : tuple[float, float] = (5.0, np.inf)
    """
    Low and high cutoff limits of binary masses in the detector frame.
    Defaults to (5.0, np.inf).
    """

    source_redshift_range : tuple[float, float] = (1e-3, 10.0)
    """
    Low and high limits of the source redshift.
    Defaults to (1e-3, 10.0).
    """

    lens_mass_range : tuple[float, float] = (1e2, 1e5)
    """
    Low and high limits of the lens mass.
    Defaults to (1e2, 1e5).
    """

    impact_parameter_range : tuple[float, float] = (1e-2, 1.25)
    """
    Low and high limits of the impact parameter.
    Defaults to (1e-2, 1.25).
    """

    geocent_time_range : tuple[float, float] = (1261872018.0, 2208643218.0)
    """
    Geocenter trigger time range (GPS time in seconds).
    Defaults to (1261872018.0, 2208643218.0), which is 01/01/2020 to 01/01/2050.
    """