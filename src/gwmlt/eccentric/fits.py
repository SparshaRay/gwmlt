"""
Tools and Functions Related to Phenomenological Fits of TEOBResumS Evolutions
"""

# [NOTE] :
# All the teobresumsfits functions except the `teobresumsfits_generalized_initconds` function
# operate with the dimensionless time and frequency in logarithmic (base 10) scale.


import warnings

import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.optimize import root_scalar

from gwmlt.config import config
from gwmlt.utils.physics import DimfulToDimless, DimlessToDimful
from gwmlt.eccentric.analytic import ecc_from_envelop_freqs


# Setting Up the Interpolators --------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

class AnisotropicTEFPhenom :
    
    def __init__(
        self, 
        x_vals : np.ndarray, 
        y_vals : np.ndarray, 
        z_vals : np.ndarray, 
        fit_vals : np.ndarray,
        smoothing : float = 0.025, 
        neighbors : int = 256,
        poly_interp_low : float = np.log10(5e2),
        poly_interp_high : float | None = None
    ) -> None :

        """
        Initialize an anisotropic R3 -> R1 RBF interpolator for TEF.

        Parameters
        ----------
        x_vals : array-like
            The 1st dimension coordinates of the data points.
        y_vals : array-like
            The 2nd dimension coordinates of the data points.
        z_vals : array-like
            The 3rd dimension coordinates of the data points.
        fit_vals : array-like
            The values to be interpolated at the training data points.
        smoothing : float, optional
            Smoothing parameter for the RBF interpolator. Default is 0.025.
        neighbors : int, optional
            Number of nearest neighbors to use for interpolation. Default is 256.
        poly_interp_low : float, optional
            The lower bound for datapoint polynomial interpolation. Default is np.log10(5e2).
        poly_interp_high : float | None, optional
            The upper bound for polynomial interpolation. If None, the query x_val is used.
            Default is None.
        """

        self.poly_interp_low  = poly_interp_low
        self.poly_interp_high = poly_interp_high

        assert x_vals.shape[:3] == y_vals.shape[:3] == z_vals.shape[:3] == fit_vals.shape[:3], \
            "Input arrays must have the same shape in the first three dimensions."
        
        self.anchor_points = fit_vals.shape[3] if len(fit_vals.shape)==4 else 1

        points = np.vstack((
            np.array(x_vals).flatten(), 
            np.array(y_vals).flatten(), 
            np.array(z_vals).flatten()
        )).T

        values = np.array(fit_vals).reshape(-1, self.anchor_points)

        self.p_min = points.min(axis=0)
        self.p_max = points.max(axis=0)
        self.p_range = self.p_max - self.p_min
        
        scaled_points = (points - self.p_min) / self.p_range
        
        self.rbf = RBFInterpolator(
            y = scaled_points,
            d = values,
            kernel    = 'thin_plate_spline',
            smoothing = smoothing,
            neighbors = neighbors
        )
        
    def __call__(self, query_points: np.ndarray) -> np.ndarray :

        """
        Queries the interpolator at given points.
        
        Parameters
        ----------
        query_points: array-like of shape (M, 3) or (3,)
            Coordinates of the points where the interpolation is to be evaluated.
        
        Returns
        -------
        array-like or scalar
            Interpolated values at the query points.
        
        Notes
        -----
        As the data is setup, the three coordinates, in order, are: 
        1 .log_waveform_duration_dimensionless or log_requested_f_start_dimensionless, 
        2. ecc_start, 
        3. chi_eff.
        """

        # Required for brentq
        q_points = [np.array(i).flatten() for i in query_points]
        q_points = [i.item() if len(i)==1 else i for i in q_points]

        q_points = np.atleast_2d(q_points)

        q_min = q_points.min(axis=0)
        q_max = q_points.max(axis=0)

        if np.any(q_min < self.p_min) or np.any(q_max > self.p_max):
            warnings.warn("Query points are outside the range of the training data.\n" 
                 "Extrapolation may lead to unreliable results.")

        q_scaled = (q_points - self.p_min) / self.p_range
        res = self.rbf(q_scaled)

        if res.shape[0]==1 or res.shape[1]==1 : res = res.flatten()
        if len(res) == 1 : res = res.item()
        return res
        
    def poly(self, x_val: float, y_val: float, z_val: float, deg:int=6) -> np.polynomial.Polynomial :

        """
        Return polynomial fit at the given point.
        Only valid for non-scalar data.
        
        Parameters
        ----------
        x_val : float
            The x-coordinate of the point.
        y_val : float
            The y-coordinate of the point.
        z_val : float
            The z-coordinate of the point.
        deg : int
            The degree of the polynomial fit.
        
        Returns
        -------
        np.polynomial.Polynomial
            Polynomial fit at the given point.
        
        Notes
        -----
        As the data is setup, the three coordinates, in order, are:
        1 .log_waveform_duration_dimensionless or log_requested_f_start_dimensionless, 
        2. ecc_start, 
        3. chi_eff.
        """

        if self.anchor_points == 1 :
            raise ValueError("Polynomial fit is not available for scalar data.")

        anchor_ys = self.__call__([x_val, y_val, z_val])
        anchor_xs = np.linspace(self.poly_interp_low, self.poly_interp_high or x_val, self.anchor_points)

        return np.polynomial.Polynomial.fit(anchor_xs, anchor_ys, deg=deg)


# Initializing the Interpolators 

_interpolators = {}
_current_teobresums_fits_path = None

# The whole interpolator initialization is wrapped in a function for 
# integrating nicely with the config_override context manager

def get_interpolators() :

    global _interpolators, _current_teobresums_fits_path

    if _interpolators=={} or _current_teobresums_fits_path!=config.teobresums_fits_path :

        with np.load(config.teobresums_fits_path) as f : data = dict(f)

        # min_log_tau_interp sets the lower bound for datapoints and polynomial interpolation.
        # All logarithms are base 10.

        _interpolators = {

            "log_f_peri_bar_interpolator" : AnisotropicTEFPhenom(
                x_vals          = data['log_waveform_duration_dimensionless'],
                y_vals          = data['ecc_start'],
                z_vals          = data['chi_eff'],
                fit_vals        = data['log_f_peri_anchors_dimensionless'],
                poly_interp_low = data['min_log_tau_interp'],
            ),
            # Interpolator for log of dimensionless periastron frequencies.
            # The three coordinates are (respectively) : 
            # log_waveform_duration_dimensionless, ecc_start, chi_eff.

            "log_f_apos_bar_interpolator" : AnisotropicTEFPhenom(
                x_vals          = data['log_waveform_duration_dimensionless'],
                y_vals          = data['ecc_start'],
                z_vals          = data['chi_eff'],
                fit_vals        = data['log_f_apos_anchors_dimensionless'],
                poly_interp_low = data['min_log_tau_interp'],
            ),
            # Interpolator for log of dimensionless apastron frequencies.
            # The three coordinates are (respectively) : 
            # log_waveform_duration_dimensionless, ecc_start, chi_eff.

            "log_f_orbavg_bar_interpolator" : AnisotropicTEFPhenom(
                x_vals          = data['log_waveform_duration_dimensionless'],
                y_vals          = data['ecc_start'],
                z_vals          = data['chi_eff'],
                fit_vals        = data['log_f_orbavg_anchors_dimensionless'],
                poly_interp_low = data['min_log_tau_interp'],
            ),
            # Interpolator for log of dimensionless orbital average frequencies.
            # The three coordinates are (respectively) : 
            # log_waveform_duration_dimensionless, ecc_start, chi_eff.

            "log_f_calib_bar_interpolator" : AnisotropicTEFPhenom(
                x_vals          = data['log_requested_f_start_dimensionless'],
                y_vals          = data['ecc_start'],
                z_vals          = data['chi_eff'],
                fit_vals        = data['log_measured_f_start_dimensionless'],
                poly_interp_low = data['min_log_tau_interp'],
            ),
            # Interpolator for calibrating log of dimensionless starting frequencies.
            # This is used to correct f_start offsets in TEOBResums.
            # The three coordinates are (respectively) : 
            # log_requested_f_start_dimensionless, ecc_start, chi_eff.
            # The output is the log_measured_f_start_dimensionless.
        }

        _current_teobresums_fits_path = config.teobresums_fits_path

    return _interpolators


# TEOBResumS Phenomenological Fits ----------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def teobresumsfits_tau_at_fbar(
    log_sig_tau : float,
    ecc_start : float,
    chi_eff : float,
    log_fbar : float,
    extrapolate : bool = False,
    suppress_warnings : bool = False
) -> float :

    """
    Returns the log of dimensionless time (tau) at which 
    the log of dimensionless frequency (log_fbar) is reached; 
    given the starting eccentricity (ecc_start), effective spin (chi_eff), 
    and log of dimensionless total signal duration (log_sig_tau).
    Uses the TEOBResums phenomenological fits. 
    This is more accurate than the Gergely ODEs, especially at high eccentricities.

    Parameters
    ----------
    log_sig_tau : float
        Logarithm of the dimensionless total signal duration.
    ecc_start : float
        Starting eccentricity.
    chi_eff : float
        Effective spin parameter.
    log_fbar : float
        Logarithm of the dimensionless frequency f_bar for which to get the time.
    extrapolate : bool, optional
        Whether to allow extrapolation outside the training data range.
        Default is False.
    suppress_warnings : bool, optional
        Whether to suppress warnings about no valid roots found within the time bracket.
        Default is False.

    Returns
    -------
    float
        Log of the dimensionless time (tau) at which 
        the specified log dimensionless frequency (log_fbar) is reached.

    Note
    ----
    All logarithms are base 10.
    """

    orbavg_poly = get_interpolators()['log_f_orbavg_bar_interpolator'].poly(log_sig_tau, ecc_start, chi_eff)
    all_roots   = (orbavg_poly - log_fbar).roots()
    real_roots  = all_roots[np.isreal(all_roots)].real
    if len(real_roots) == 0 : raise ValueError("No real roots found for the given parameters.")
    interp_low  = get_interpolators()['log_f_orbavg_bar_interpolator'].poly_interp_low
    valid_root  = real_roots[np.all([real_roots>=interp_low, real_roots<=log_sig_tau], axis=0)]

    if len(valid_root) == 1 : return valid_root[0]
    if len(valid_root) > 1  : raise ValueError("Multiple roots found in time bracket")
    if len(valid_root) == 0 :
        if not extrapolate : 
            if not suppress_warnings :
                warnings.warn("No valid roots found in time bracket, fallback to closest bound.")
            return log_sig_tau if orbavg_poly(log_sig_tau)>log_fbar else interp_low
        else :
            if not suppress_warnings :
                warnings.warn("No valid roots found in time bracket, extrapolating outside bound.")
            if orbavg_poly(log_sig_tau)>log_fbar and np.any(real_roots>log_sig_tau) :
                return real_roots[real_roots>log_sig_tau].min()
            else : return real_roots[real_roots<interp_low].max()


def teobresumsfits_ecc_at_fbar(
    log_sig_tau : float,
    ecc_start : float,
    chi_eff : float,
    log_fbar : float,
    extrapolate : bool = False,
    suppress_warnings : bool = False
) -> float :

    """
    Returns the eccentricity (ecc) when 
    the log of dimensionless frequency (log_fbar) is reached; 
    given the starting eccentricity (ecc_start), effective spin (chi_eff), 
    and log of dimensionless total signal duration (log_sig_tau).
    Uses the TEOBResums phenomenological fits. 
    This is more accurate than the Gergely ODEs, especially at high eccentricities.

    Parameters
    ----------
    log_sig_tau : float
        Logarithm of the dimensionless total signal duration.
    ecc_start : float
        Starting eccentricity.
    chi_eff : float
        Effective spin parameter.
    log_fbar : float
        Logarithm of the dimensionless frequency f_bar for which to get the eccentricity.
    extrapolate : bool, optional
        Whether to allow extrapolation outside the training data range.
        Default is False.
    suppress_warnings : bool, optional
        Whether to suppress warnings about no valid roots found within the time bracket.
        Default is False.

    Returns
    -------
    float
        The eccentricity (ecc) when the specified 
        log of dimensionless frequency (log_fbar) is reached.

    Note
    ----
    All logarithms are base 10.
    """
    
    tau = teobresumsfits_tau_at_fbar(
        log_sig_tau, ecc_start, chi_eff, log_fbar, 
        extrapolate=extrapolate, suppress_warnings=suppress_warnings
    )
    f_peri_bar = 10 ** get_interpolators()['log_f_peri_bar_interpolator'].poly(log_sig_tau, ecc_start, chi_eff)(tau)
    f_apos_bar = 10 ** get_interpolators()['log_f_apos_bar_interpolator'].poly(log_sig_tau, ecc_start, chi_eff)(tau)
    return ecc_from_envelop_freqs(f_apos_bar, f_peri_bar)


def teobresumsfits_initial_conditions(
    log_sig_tau : float,
    chi_eff : float,
    ecc_ref : float,
    log_f_ref_bar : float,
    extrapolate : bool = False
) -> tuple[float, float] :
    
    """
    Get the starting eccentricity (ecc_start) and 
    log of the dimensionless starting frequency (log_f_start_bar);
    given the log of dimensionless total signal duration (log_sig_tau),
    effective spin (chi_eff), and reference eccentricity (ecc_ref) at 
    the log of dimensionless reference frequency (log_f_ref_bar).

    Based on a database of highly faithful phenomenological fits to TEOBResumS.
    This method holds better accuracies than the Gergely ODEs at high eccentricities.

    Parameters
    ----------
    log_sig_tau : float
        Logarithm of the dimensionless total signal duration.
    chi_eff : float
        Effective spin parameter.
    ecc_ref : float
        Reference eccentricity at the reference frequency.
    log_f_ref_bar : float
        Logarithm of the dimensionless reference frequency f_ref_bar.
    extrapolate : bool, optional
        Whether to allow extrapolation outside the training data range.
        If you get "initial conditions failed" errors, try setting this to True.
        Default is False.
    
    Returns
    -------
    tuple[float, float]
        Tuple containing the starting eccentricity (ecc_start) and 
        log of the dimensionless starting frequency (log_f_start_bar).
    
    Note
    ----
    All logarithms are base 10.
    If extrapolate is set to False, log_f_ref_bar must be reached after 
    the start time. Set extrapolate to True to enable backwards interpolation.
    """
    
    # Find the starting eccentricity for which we will get the 
    # reference eccentricity at the reference frequency.
    def ecc_error(ecc_start) :
        ecc_guess = teobresumsfits_ecc_at_fbar(
            log_sig_tau, ecc_start, chi_eff, log_f_ref_bar, 
            extrapolate=extrapolate, suppress_warnings=True
            # Just let brentq query whatever point it wants, dont throw warnings here.
        )
        return ecc_guess - ecc_ref
    
    # Since ecc is the 2nd argument (or coordinate) of the interpolators, 
    # we can use the min and max values of the 2nd argument as the bracket for root finding.
    ecc_bracket = [get_interpolators()['log_f_orbavg_bar_interpolator'].p_min[1], 
                   get_interpolators()['log_f_orbavg_bar_interpolator'].p_max[1]]
    
    try : ecc_start = root_scalar(ecc_error, bracket=ecc_bracket, method='brentq').root
    except ValueError as e : raise ValueError(f"Initial conditions failed: {e}")
    
    required_log_fbar_start = get_interpolators()['log_f_orbavg_bar_interpolator'](
        [log_sig_tau, ecc_start, chi_eff])[-1] # Fetch the last index, i.e. starting freq anchor point
    
    # Verify the results
    with warnings.catch_warnings() :
        if not extrapolate : warnings.filterwarnings("error")
        # Catch the warning from AnisotropicTEFPhenom.__call__ range check
        # when the query points are outside the training data range.
        try :
            ecc_obt = teobresumsfits_ecc_at_fbar(
                log_sig_tau, ecc_start, chi_eff, log_f_ref_bar, extrapolate=extrapolate)
        except Warning as w : 
            if not extrapolate : raise Exception(f"Initial conditions failed: {w}")
        if np.abs(ecc_obt - ecc_ref) > 1e-2 : raise Exception("Initial conditions failed")

    # Since TEOBResumS don't start the waveform at the requested f_start, 
    # we need to find an starting frequency that will give us our desired f_start.
    def f_calib_error(requested_log_fbar_start) :
        measured_log_fbar_start = get_interpolators()['log_f_calib_bar_interpolator'](
            [requested_log_fbar_start, ecc_start, chi_eff])
        return measured_log_fbar_start - required_log_fbar_start
    
    try :
        log_f_start_bar = root_scalar(
            f_calib_error, 
            # This bracket is usually far more than enough. Decrease if necessary.
            bracket=[required_log_fbar_start-0.25, required_log_fbar_start+0.25],
            method='brentq'
        ).root
    except ValueError as e : raise ValueError(f"Initial conditions failed: {e}")

    return ecc_start, log_f_start_bar


def teobresumsfits_generalized_initconds(
    mass_1 : float,
    mass_2 : float,
    chi1z  : float,
    chi2z  : float,
    ecc_ref : float,
    f_ref   : float,
    waveform_duration : float,
    extrapolate : bool = False
) -> tuple[float, float] :
    
    """
    Get the starting eccentricity (ecc_start) and starting frequency (f_start)
    for TEOBResumS waveform generation.

    Parameters
    ----------
    mass_1 : float
        Mass of the primary in solar masses.
    mass_2 : float
        Mass of the secondary in solar masses.
    chi1z : float
        Dimensionless spin of the primary along the orbital angular momentum.
    chi2z : float
        Dimensionless spin of the secondary along the orbital angular momentum.
    ecc_ref : float
        Reference eccentricity at the reference frequency.
    f_ref : float
        Reference frequency in Hz at which the reference eccentricity is defined.
    waveform_duration : float
        Duration of the waveform in seconds.
    extrapolate : bool, optional
        Whether to allow extrapolation outside the training data range.
        If you get "initial conditions failed" errors, try setting this to True.
        Setting this to True may lead to inaccurate results. Default is False.
    
    Returns
    -------
    tuple[float, float]
        Tuple containing the starting eccentricity (ecc_start) and 
        starting frequency (f_start) in Hz for TEOBResumS.
    """

    eta     = (mass_1 * mass_2) / ((mass_1 + mass_2)**2)
    chi_eff = (mass_1*chi1z + mass_2*chi2z) / (mass_1 + mass_2)

    to_dimensionless = DimfulToDimless(mass_1, mass_2)
    f_ref_bar = to_dimensionless.get_dimless_frequency(f_ref)
    sig_tau   = to_dimensionless.get_dimless_time(waveform_duration)

    # Account for signal duration scaling with symmetric mass ratio
    calibrated_sig_tau = sig_tau * (4 * eta) 

    log_sig_tau   = np.log10(calibrated_sig_tau)
    log_f_ref_bar = np.log10(f_ref_bar)

    ecc_start, log_f_start_bar = teobresumsfits_initial_conditions(
        log_sig_tau=log_sig_tau,
        chi_eff=chi_eff,
        ecc_ref=ecc_ref,
        log_f_ref_bar=log_f_ref_bar,
        extrapolate=extrapolate
    )

    to_physical = DimlessToDimful(mass_1 + mass_2)
    f_start_bar = 10 ** log_f_start_bar
    f_start     = to_physical.get_dimful_frequency(f_start_bar)

    return ecc_start, f_start