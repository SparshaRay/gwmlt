"""
Lensing by Isolated Point-mass Lenses
"""

import numpy as np
from scipy.special import loggamma

from ..core.morphologies import Lensed
from ..constants import C_SI, G_SI, MSUN_SI


# Kummer’s function -----------------------------------------------------------

from juliacall import Main as julia

julia.seval("""
using HypergeometricFunctions
function hyp1f1(a_vec, b_vec, z_vec)
    return HypergeometricFunctions.pFq.(Tuple.(a_vec), Tuple.(b_vec), z_vec)
end
""")

hyp1f1 = julia.hyp1f1
"""
Robust, fast, and vectorized implementation of the
confluent hypergeometric function of the first kind.
"""


# Images ----------------------------------------------------------------------

def _x_minima(y):
    """Returns the image position for the minima (type I) image."""
    return (y + np.sqrt(y**2 + 4)) / 2

def _x_saddle(y):
    """Returns the image position for the saddle-point (type II) image."""
    return (y - np.sqrt(y**2 + 4)) / 2

def _magnification_minima(y):
    """Returns the image magnification for the minima (type I) image."""
    return 1 / 2 + (y**2 + 2) / (2 * y * np.sqrt(y**2 + 4))

def _magnification_saddle(y):
    """Returns the image magnification for the saddle-point (type II) image."""
    return 1 / 2 - (y**2 + 2) / (2 * y * np.sqrt(y**2 + 4))


# Time delay ------------------------------------------------------------------

def _time_delay_dimensionless(y):
    """Returns the dimensionless time-delay between the two micro-images."""
    sqrt_term = np.sqrt(y**2 + 4)
    return (y * sqrt_term) / 2.0 + np.log((sqrt_term + y) / (sqrt_term - y))

# Expose if necessary
def time_delay(
    ml : float,
    y  : float,
    zl : float = 0
) -> float :

    """
    Returns the dimensionful time-delay between the two micro-images in seconds.

    Parameters
    ----------
    ml : float
        Lens mass in solar masses.
    y : float
        Dimensionless impact parameter (in units of the Einstein radius).
    zl : float, optional
        Lens redshift. Default is 0.

    Returns
    -------
    float
        Time-delay between the two micro-images in seconds.
    
    """

    return (4 * G_SI * MSUN_SI * ml * (1 + zl) / C_SI**3) * _time_delay_dimensionless(y)

# Primary time delay user entry point function
def get_time_delay(morphology: Lensed) -> float :
    """
    Returns the time-delay between the two micro-images for a given lensing morphology.

    Parameters
    ----------
    morphology : Lensed
        The lensing morphology containing lensing parameters.

    Returns
    -------
    float
        Time-delay between the two micro-images in seconds.
    """
    return time_delay(morphology.m_lens, morphology.y_lens, morphology.z_lens)


# Frequency conversions -------------------------------------------------------

F_CONST = 8.0 * np.pi * G_SI * MSUN_SI / C_SI**3

def _w_of_f(f, ml, zl=0):
    """Converts a dimensionful frequency (f in Hz) to dimensionless frequency (w)."""
    return f * F_CONST * ml * (1 + zl)

def _f_of_w(w, ml, zl=0):
    """Converts a dimensionless frequency (w) to dimensionful frequency (f in Hz)."""
    return w / (F_CONST * ml * (1 + zl))


# Regime transition cutoffs ---------------------------------------------------
# Only valid for y in (0.01, 5.00)

def _w_cutoff_geometric(y):
    """Dimensionless cutoff frequency where Geometric optics error < 0.1%."""
    if   y <= 0.12:
        return 15112.5 - 52563.5 * y
    elif y <= 1.50:
        return -34.08 - 12.84 * (y**-1.0) + 114.33 * (y**-2.0) + 0.89 * (y**-3.0)
    else :
        return -15.02 + 18.25 * y - 2.66 * y**2

# Only kept for backwards compatibility.
def _w_cutoff_geometric_tolerance_1p0(y):
    """Dimensionless cutoff frequency where Geometric optics error < 1.0%."""
    if y <= 0.071 : return 16604 - 202686 * y
    return 0.64 + 0.97 * (y**-1.0) + 6 * (y**-2.0) + 0.38 * (y**-3.0)

def _w_cutoff_quasigeometric(y):
    """Dimensionless cutoff frequency where Quasi-geometric optics error < 0.1%."""
    return 9 * (y**-1.0) + 0.04 * (y**-2.0)

# Only kept for backwards compatibility.
def _w_cutoff_quasigeometric_tolerance_1p0(y):
    """Dimensionless cutoff frequency where Quasi-geometric optics error < 1.0%."""
    return 4 * (y**-1.0) - np.log(y) / 5.0


# Lensing amplification in different regimes ----------------------------------

def _Fw_geometric_optics(w, y):
    """Amplification factor F(w) under the Geometric Optics approximation."""
    return (
        np.sqrt(np.abs(_magnification_minima(y))) 
        - 1j * np.sqrt(np.abs(_magnification_saddle(y)))
        * np.exp(1j * w * _time_delay_dimensionless(y))
    )

def _Fw_quasigeometric_optics(w, y):
    """Amplification factor F(w) under the Quasi-Geometric Optics approximation."""
    xm = _x_minima(y)
    xs = _x_saddle(y)
    
    term_minima = (4 * xm**2 - 1) / (((xm**2 + 1)**3) * (xm**2 - 1))
    term_saddle = (4 * xs**2 - 1) / (((xs**2 + 1)**3) * (xs**2 - 1))
    
    return (
        _Fw_geometric_optics(w, y)
        + (1j / (3 * w)) * term_minima * np.sqrt(np.abs(_magnification_minima(y)))
        + (1  / (3 * w)) * term_saddle * np.sqrt(np.abs(_magnification_saddle(y)))
        * np.exp(1j * w * _time_delay_dimensionless(y))
    )

def _Fw_exact(w, y):
    """Returns the exact lensing amplification factor using the analytical formula.
    As the system approaches geometric optics limit, numerical instabilities may arise.
    Under such conditions, it is recommended to use the _Fw_hybrid method."""

    w_arr = np.atleast_1d(w)
    
    result = np.ones(w_arr.shape, dtype=complex)
    nz = w_arr != 0
    if not np.any(nz) : return result if np.ndim(w)>0 else result[0]

    w_nz = w_arr[nz]

    xm = (y + np.sqrt(y * y + 4.0)) / 2.0
    pm = ((xm - y) ** 2) / 2.0 - np.log(xm)
    hp = np.log(w_nz / 2.0) - (2.0 * pm)

    log_h  = (np.pi * w_nz / 4.0) + 1j * (hp * w_nz / 2.0)
    log_gm = loggamma(1.0 - (1j * w_nz / 2.0))
    log_hf = np.log(hyp1f1(
        1j * w_nz / 2.0, 
        np.ones_like(w_nz), 
        1j * y * y * w_nz / 2.0
    ))
    
    result[nz] = np.exp(log_h + log_gm + log_hf)
    return result if np.ndim(w)>0 else result[0]


# Dimensionful variants of the above functions

def Ff_geometric_optics(f, ml, y, zl=0):
    """Dimensionful variant of _Fw_geometric_optics."""
    return _Fw_geometric_optics(_w_of_f(f, ml, zl), y)

def Ff_quasigeometric_optics(f, ml, y, zl=0):
    """Dimensionful variant of _Fw_quasigeometric_optics."""
    return _Fw_quasigeometric_optics(_w_of_f(f, ml, zl), y)

def Ff_exact(f, ml, y, zl=0):
    """Dimensionful variant of _Fw_exact."""
    return _Fw_exact(_w_of_f(f, ml, zl), y)


# Lensing Amplification Factor with Hybrid Regime Switching -------------------

def _Fw_hybrid(w, y):
    """
    Point-lens amplification factor F(w) that automatically switches 
    between exact, quasi-geometric, and geometric optics regimes.
    """

    w_arr = np.atleast_1d(w)
    
    if y<0.01 or y>5.00 :
        return _Fw_exact(w_arr, y) if np.ndim(w)>0 else _Fw_exact(w_arr, y)[0]
        
    wc_geo  = _w_cutoff_geometric(y)
    wc_qgeo = _w_cutoff_quasigeometric(y)
    
    result = np.empty(w_arr.shape, dtype=complex)
    
    if wc_qgeo < wc_geo :
        cond_exact = (w_arr < wc_qgeo)
        cond_geo   = (w_arr >= wc_geo)
        cond_qgeo  = ~(cond_exact | cond_geo)
        
        if np.any(cond_exact) :
            result[cond_exact] = _Fw_exact(w_arr[cond_exact], y)
        if np.any(cond_qgeo) :
            result[cond_qgeo] = _Fw_quasigeometric_optics(w_arr[cond_qgeo], y)
        if np.any(cond_geo) :
            result[cond_geo] = _Fw_geometric_optics(w_arr[cond_geo], y)

    else : result[:] = _Fw_quasigeometric_optics(w_arr, y)
        
    return result if np.ndim(w)>0 else result[0]


# Primary Ff user entry point function
def Ff_hybrid(
    f  : float | np.typing.NDArray[np.float64],
    ml : float,
    y  : float,
    zl : float = 0
) -> float | np.typing.NDArray[np.complex128] :
    
    """
    Numerically stable evaluation of the point-lens amplification factor F(f)
    for given dimensionful frequencies f, lens mass, lens redshift, and impact parameter.

    Parameters
    ----------
    f : float or array_like
        Dimensionful frequency in Hz.
    ml : float
        Lens mass (in solar masses).
    y : float
        Dimensionless impact parameter (in units of the Einstein radius).
    zl : float, optional
        Lens redshift. Default is 0.

    Returns
    -------
    float or array_like
        Complex amplification factor F(f) for the given parameters.
    """

    return _Fw_hybrid(_w_of_f(f, ml, zl), y)