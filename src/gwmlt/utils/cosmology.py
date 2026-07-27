"""
Cosmology Utilities
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d

from gwmlt.constants import C_SI
from gwmlt.config import config


# Basic functions ---------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def inverse_ez(z: float) -> float :

    """
    Returns 1/E(z), where E(z) = H(z)/H0 is the dimensionless Hubble parameter.
    This serves as the integration kernel for distance calculations.

    Parameters:
        z (float): Redshift.
    
    Reference: 
        - Eq. 14, Hogg (2000): https://arxiv.org/pdf/astro-ph/9905116.pdf
    """

    OMEGA_M      = config.cosmology.OMEGA_M
    OMEGA_LAMBDA = config.cosmology.OMEGA_LAMBDA
    OMEGA_K      = config.cosmology.OMEGA_K
    OMEGA_R = 1.0 - OMEGA_M - OMEGA_LAMBDA - OMEGA_K 

    ez = np.sqrt(
        OMEGA_R * (1 + z)**4 + 
        OMEGA_M * (1 + z)**3 + 
        OMEGA_K * (1 + z)**2 + 
        OMEGA_LAMBDA
    )

    return 1.0 / ez


def comoving_distance(z: float, z_start: float = 0.0) -> float :

    """
    Returns the line-of-sight comoving distance to redshift z, 
    or between two redshifts if z_start is provided.

    Reference: 
        - Eq. 15, Hogg (2000)

    Parameters:
        z (float): Target redshift.
        z_start (float, optional): Starting redshift. Defaults to 0.0 (observer).

    Returns:
        float: Comoving distance in Mpc.
    """

    HUBBLE_CONSTANT = config.cosmology.HUBBLE_CONSTANT
    HUBBLE_DISTANCE_MPC = C_SI / (HUBBLE_CONSTANT * 1000.0) 

    integral, _ = quad(inverse_ez, z_start, z)
    return HUBBLE_DISTANCE_MPC * integral


def transverse_comoving_distance(z: float) -> float :

    """
    Returns the transverse comoving distance, 
    accounting for spatial curvature (Omega_k).
    
    Reference: 
        - Eq. 16, Hogg (2000)

    Parameters:
        z (float): Redshift.

    Returns:
        float: Transverse comoving distance in Mpc.
    """

    OMEGA_K         = config.cosmology.OMEGA_K
    HUBBLE_CONSTANT = config.cosmology.HUBBLE_CONSTANT
    HUBBLE_DISTANCE_MPC = C_SI / (HUBBLE_CONSTANT * 1000.0) 

    dc = comoving_distance(z)
    sqrt_ok = np.sqrt(abs(OMEGA_K))

    if OMEGA_K == 0.0 : 
        return dc
    if OMEGA_K >  0.0 :
        return (HUBBLE_DISTANCE_MPC / sqrt_ok) * np.sinh(sqrt_ok * dc / HUBBLE_DISTANCE_MPC)
    if OMEGA_K <  0.0 :
        return (HUBBLE_DISTANCE_MPC / sqrt_ok) * np.sin(sqrt_ok * dc / HUBBLE_DISTANCE_MPC)


def angular_diameter_distance_between(z1: float, z2: float) -> float :

    """
    Returns the angular diameter distance between two redshifts (z1 < z2).
    
    Reference: 
        - Eq. 19, Hogg (2000) 
        (Valid for all Omega_k, not just for positive Omega_k as stated in Hogg's text).
    
    Parameters:
        z1 (float): Lower redshift.
        z2 (float): Higher redshift.
    
    Returns:
        float: Angular diameter distance between z1 and z2 in Mpc.
    """

    OMEGA_K         = config.cosmology.OMEGA_K
    HUBBLE_CONSTANT = config.cosmology.HUBBLE_CONSTANT
    HUBBLE_DISTANCE_MPC = C_SI / (HUBBLE_CONSTANT * 1000.0) 

    dm1 = transverse_comoving_distance(z1)
    dm2 = transverse_comoving_distance(z2)
    dh = HUBBLE_DISTANCE_MPC
    
    term1 = dm2 * np.sqrt(1.0 + OMEGA_K * (dm1 / dh)**2)
    term2 = dm1 * np.sqrt(1.0 + OMEGA_K * (dm2 / dh)**2)
    
    return (term1 - term2) / (1.0 + z2)


def comoving_volume(z: float) -> float :

    """
    Returns the total comoving volume enclosed out to redshift z, 
    accounting for spatial curvature.
    
    Reference: 
        - Eq. 29, Hogg (2000)

    Parameters:
        z (float): Redshift.
    
    Returns:
        float: Total comoving volume in Mpc^3.
    """

    OMEGA_K         = config.cosmology.OMEGA_K
    HUBBLE_CONSTANT = config.cosmology.HUBBLE_CONSTANT
    HUBBLE_DISTANCE_MPC = C_SI / (HUBBLE_CONSTANT * 1000.0) 

    dm = transverse_comoving_distance(z)
    dh = HUBBLE_DISTANCE_MPC
    
    if OMEGA_K == 0.0 : return (4.0 * np.pi / 3.0) * dm**3
        
    term1 = 4.0 * np.pi * dh**3 / (2.0 * OMEGA_K)
    term2 = (dm / dh) * np.sqrt(1.0 + OMEGA_K * (dm / dh)**2)
    sqrt_ok = np.sqrt(abs(OMEGA_K))
    
    if OMEGA_K > 0.0 :
        return term1 * (term2 - (1.0 / sqrt_ok) * np.arcsinh(sqrt_ok * dm / dh))
    if OMEGA_K < 0.0 :
        return term1 * (term2 - (1.0 / sqrt_ok) * np.arcsin(sqrt_ok * dm / dh))


def differential_comoving_volume(z: float) -> float :

    """
    Computes the differential comoving volume element dVc/dz in Gpc^3.

    Parameters:
        z (float): Redshift.

    Returns:
        float: Differential comoving volume element (Gpc^3 per unit redshift).
    """

    HUBBLE_CONSTANT = config.cosmology.HUBBLE_CONSTANT
    HUBBLE_DISTANCE_MPC = C_SI / (HUBBLE_CONSTANT * 1000.0) 

    D_A = angular_diameter_distance_between(0.0, z)
    inv_Ez = inverse_ez(z)

    # 4 * pi * (c/H0) * (1+z)^2 * D_A^2 / E(z)
    # 1e-9 converts Mpc^3 to Gpc^3
    d_Omega = 4.0 * np.pi
    dVc_dz = (
        d_Omega
        * HUBBLE_DISTANCE_MPC
        * ((1.0 + z) ** 2)
        * (D_A ** 2)
        * inv_Ez
    ) * 1e-9

    return dVc_dz


# SFR and merger rates ----------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def madau_dickinson_fishbach_sfr(
    z: float | np.typing.NDArray[np.float64]
) -> float | np.typing.NDArray[np.float64] :

    """
    Star formation rate following Fishbach et al. (2018) [https://arxiv.org/abs/1805.10270],
    based on Madau & Dickinson (2014) [https://arxiv.org/abs/1403.0007].

    Parameters:
        z (float or numpy.ndarray): Redshift.

    Returns:
        float or numpy.ndarray: Star formation rate in Gpc^-3 yr^-1
    """

    return 0.015 * ((1.0 + z) ** 2.7) / (1.0 + ((1.0 + z) / 2.9) ** 5.6)


def madau_fragos_sfr(
    z: float | np.typing.NDArray[np.float64]
) -> float | np.typing.NDArray[np.float64] :

    """
    Star formation rate following Madau & Fragos (2016) [https://arxiv.org/abs/1606.07887].

    Parameters:
        z (float or numpy.ndarray): Redshift.

    Returns:
        float or numpy.ndarray: Star formation rate in Gpc^-3 yr^-1
    """

    return 0.010 * ((1.0 + z) ** 2.6) / (1.0 + ((1.0 + z) / 3.2) ** 6.2)


def merger_rate_density(
    z: float, 
    sfr_model: callable = madau_dickinson_fishbach_sfr
) -> float :

    """
    The merger rate density (per unit redshift) in the detector frame :
    dR_det/dz = [R_src(z) / (1 + z)] * (dVc/dz).

    Parameters:
        z (float): Redshift.
        sfr_model (callable, optional): Function returning the source-frame star formation rate. 
                                        Defaults to `madau_dickinson_fishbach_sfr`.

    Returns:
        float: Merger rate density in the detector frame.
    """

    rate_src = sfr_model(z)
    dVc_dz = differential_comoving_volume(z)
    return (rate_src * dVc_dz) / (1.0 + z)  # dt_det = dt_src * (1+z)


def merger_rate_pdf_interpolator(
    z_min: float, 
    z_max: float, 
    num_points: int = 10000,
    sfr_model: callable = madau_dickinson_fishbach_sfr
) -> interp1d :

    """
    Normalized probability density function (PDF) for merger rate in the detector frame.

    Parameters:
        z_min (float): Minimum redshift.
        z_max (float): Maximum redshift.
        num_points (int, optional): Number of grid points for interpolation. Defaults to 10000.
        sfr_model (callable, optional): The source-frame star formation rate function. 
                                        Defaults to `madau_dickinson_fishbach_sfr`.

    Returns:
        interp1d: A SciPy interpolation object representing the normalized PDF.
    """

    norm, _ = quad(merger_rate_density, z_min, z_max, args=(sfr_model,))
    z_grid = np.linspace(z_min, z_max, num=num_points)
    pdf_grid = np.array([merger_rate_density(z, sfr_model)/norm for z in z_grid])
    return interp1d(z_grid, pdf_grid, kind='linear', bounds_error=False, fill_value=0.0)


# Utility functions -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def luminosity_distance(z: float) -> float :

    """
    Returns the luminosity distance to redshift z.
    
    Reference: 
        - Eq. 21, Hogg (2000)
    
    Parameters:
        z (float): Redshift.
    
    Returns:
        float: Luminosity distance in Mpc.
    """

    return (1.0 + z) * transverse_comoving_distance(z)


def z_from_comoving_distance(target_dc: float, z_max: float = 20.0) -> float :

    """
    Finds the redshift corresponding to a given comoving distance.

    Parameters:
        target_dc (float): Target comoving distance in Mpc.
        z_max (float, optional): Maximum redshift to search for the root. Defaults to 20.0.
    
    Returns:
        float: Redshift corresponding to the target comoving distance.
    """

    return brentq(lambda z: comoving_distance(z) - target_dc, 0.0, z_max)


def z_from_luminosity_distance(target_dl: float, z_max: float = 20.0) -> float :

    """
    Finds the redshift corresponding to a given luminosity distance.

    Parameters:
        target_dl (float): Target luminosity distance in Mpc.
        z_max (float, optional): Maximum redshift to search for the root. Defaults to 20.0.
    
    Returns:
        float: Redshift corresponding to the target luminosity distance.
    """

    return brentq(lambda z: luminosity_distance(z) - target_dl, 0.0, z_max)