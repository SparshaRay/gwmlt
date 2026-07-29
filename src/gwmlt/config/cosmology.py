"""
Cosmology Constants Configuration
"""

from dataclasses import dataclass

@dataclass(frozen=True)
class CosmologyConfig :
    """
    Cosmology constants configuration.
    Based on Planck 2018 results (TT,TE,EE+lowE+lensing+BAO constraints).

    References :
        1. Planck 2018 results : https://arxiv.org/pdf/1807.06209
        2. Hogg (2000) : https://arxiv.org/abs/astro-ph/9905116
    """

    OMEGA_M : float = 0.3111
    """
    Matter (baryonic plus CDM) density parameter.
    Default is 0.3111
    """

    OMEGA_LAMBDA : float = 0.6889
    """
    Dark energy density parameter.
    Default is 0.6889
    """

    OMEGA_K : float = 0.0
    """
    Curvature density parameter.
    Default is 0.0 (flat universe)
    """

    HUBBLE_CONSTANT : float = 67.66
    """
    Hubble constant in (km/s/Mpc).
    Default is 67.66
    """