"""
Orchestrator Configurations
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class OrchestratorConfig :
    """
    Configurations for various orchestrators.
    """

    max_concurrent_jobs : int = 5000
    """
    Maximum number of concurrent HTCondor jobs. Default is 5000. 
    """

    populations_dir_name : str = "population"
    """Directory prefix for population databases. Default is 'population'."""
    populations_table_name : str = "params"
    """Table name for storing population samples. Default is 'params'."""
    populations_samples_per_worker : int = 100_000
    """Number of samples to be created by each population sampling worker. Default is 100,000."""