"""
Various Functions Related to Eccentric Systems
"""

# [NOTE] :
#
# Throughout this file, the "bar" suffix denotes dimensionless quantities.
#
# Except for the `ecc_from_envelop_freqs` and `teobresumsfits_generalized_initconds` functions
# the rest of the functions are exclusively dimensionless. To get dimensionless quantities from
# physical quantities, use the conversion classes in `.general.py`.
#
# All the teobresumsfits functions except the `teobresumsfits_generalized_initconds` function
# operate with the dimensionless time and frequency in logarithmic (base 10) scale.


import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

from .general import DimfulToDimless, DimlessToDimful
from .teobresums import (log_f_peri_bar_interpolator, 
                         log_f_apos_bar_interpolator, 
                         log_f_orbavg_bar_interpolator, 
                         log_f_calib_bar_interpolator)


# General Functions -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def ecc_from_envelop_freqs(
    f_apos : float | np.ndarray, 
    f_peri : float | np.ndarray
) -> float | np.ndarray :
    
    """
    Get eccentricity from apocenter and pericenter frequencies.\n
    ref: https://arxiv.org/pdf/2302.11257 (Eq. 4 to 9)

    Parameters
    ----------
    f_apos : float | np.ndarray
        Apocenter frequencies.
    f_peri : float | np.ndarray
        Pericenter frequencies.
    
    Returns
    -------
    float | np.ndarray
        Eccentricity corresponding to the given apocenter and pericenter frequencies.

    Note
    ----
    The function can handle both scalar and array inputs for f_apos and f_peri.
    The output will be a scalar if the inputs are scalars, and an array if the 
    inputs are arrays of the same shape. The frequencies may both be dimensionless or physical.
    """

    f_apos = np.asarray(f_apos)
    f_peri = np.asarray(f_peri)

    assert f_apos.shape == f_peri.shape, "Input arrays must have the same shape."
    
    ecc = np.zeros_like(f_apos, dtype=float)
    valid_mask = f_peri > f_apos
    
    f_a = f_apos[valid_mask]
    f_p = f_peri[valid_mask]
    
    if np.any(valid_mask) :

        sqrt_p = np.sqrt(f_p)
        sqrt_a = np.sqrt(f_a)
        
        ecc_omg22 = (sqrt_p - sqrt_a) / (sqrt_p + sqrt_a)
        phi       = np.arctan((1.0 - ecc_omg22**2) / (2.0 * ecc_omg22))
        
        ecc[valid_mask] = np.cos(phi / 3.0) - np.sqrt(3.0) * np.sin(phi / 3.0)
        
    return ecc.item() if ecc.ndim == 0 else ecc


# Gergely Functions -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def GergelyODEs(
    tau : float, 
    state : tuple[float, float], 
    nu : float, 
    q : float, 
    S1z_bar : float, 
    S2z_bar : float, 
    e_term : float = 1e-5
) -> tuple[float, float]:

    """
    Dimensionless da/dt and de/dt up to 1.5PN.
    Assumes L and S are aligned. If necessary add the dk/dt term for precessing systems.
    Based on Gergely (1998). [https://arxiv.org/pdf/gr-qc/9808063]

    Parameters
    ----------
    tau : float
        Dimensionless time. Not used in the equations, but required for solvers.
    state : tuple[float, float]
        Tuple containing the values of [a_bar, e],
        where a_bar is the dimensionless generalized semimajor axis and e is the eccentricity.
    nu : float
        Dimensionless symmetric mass ratio, i.e., m1*m2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., chi1 / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., chi2 / (m2 / m_tot)^2
    e_term : float, optional
        Terminal eccentricity, the threshold below which the system is considered circularized.
        Default is 1e-5.

    Returns
    -------
    tuple[float, float]
        Tuple containing the derivatives (da_bar/dtau, de/dtau).
    """

    # Setup --------------------------------------------------------------

    # Unpack the state variables
    a_bar, e = state

    # Handle terminal conditions and unphysical values
    if (
        (e > 0.99)      or
        (e < e_term)    or 
        np.isnan(e)     or
        (a_bar < 1.00)  or 
        np.isnan(a_bar)
    ) : return (0.0, 0.0)
    
    # Precompute powers of e and (1 - e^2)
    e2=e**2; e4=e**4; e6=e**6
    one_minus_e2 = 1.0 - e2
    
    # d(a_bar)/d(tau) ----------------------------------------------------

    # 0PN Term (same as Peters)
    da_dtau_0PN = - (2 * nu * (37*e4 + 292*e2 + 96)) / \
                    (15 * a_bar**3 * (one_minus_e2)**3.5)
    
    # 1.5PN Term
    da_dtau_15PN_coeff = nu / (15 * a_bar**4.5 * (one_minus_e2)**5)
    da_dtau_15PN_poly1 = (363*e6 + 3510*e4 + 7936*e2 + 2128) * (S1z_bar + S2z_bar)
    da_dtau_15PN_poly2 = (291*e6 + 4224*e4 + 7924*e2 + 1680) * (q * S1z_bar + (1/q) * S2z_bar)
    
    da_dtau = da_dtau_0PN + da_dtau_15PN_coeff * (da_dtau_15PN_poly1 + da_dtau_15PN_poly2)
    
    # d(e)/d(tau) --------------------------------------------------------

    # 0PN Term (same as Peters)
    de_dtau_0PN = - (nu * e * (121*e2 + 304)) / \
                    (15 * a_bar**4 * (one_minus_e2)**2.5)
    
    # 1.5PN Term
    de_dtau_15PN_coeff = (nu * e) / (30 * a_bar**5.5 * (one_minus_e2)**4)
    de_dtau_15PN_poly1 = (1313*e4 + 5592*e2 + 7032) * (S1z_bar + S2z_bar)
    de_dtau_15PN_poly2 = (1097*e4 + 6822*e2 + 6200) * (q * S1z_bar + (1/q) * S2z_bar)
    
    de_dtau = de_dtau_0PN + de_dtau_15PN_coeff * (de_dtau_15PN_poly1 + de_dtau_15PN_poly2)

    # Return the derivatives ---------------------------------------------

    return (da_dtau, de_dtau)


def GergelyEvolve(
    nu : float, 
    q : float,
    S1z_bar : float, 
    S2z_bar : float,
    f_start_bar : float, 
    ecc_start : float,
    max_tau : float = 1e8,
    e_term : float = 1e-5,
    **kwargs
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] :
    
    """
    Evolves the Gergely ODEs in dimensionless time.

    Parameters
    ----------
    nu : float
        Dimensionless symmetric mass ratio, i.e., m1*m2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., chi1 / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., chi2 / (m2 / m_tot)^2
    f_start_bar : float
        Starting dimensionless GW frequency.
    ecc_start : float
        Starting eccentricity.
    max_tau : float, optional
        Maximum dimensionless time to evolve to.
        Default is 1e8.
    e_term : float, optional
        Terminal eccentricity, the threshold below which the system is considered circularized.
        Default is 1e-5.
    **kwargs
        Additional keyword arguments to pass to the ODE solver.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Tuple containing the arrays of (tau, a_bar, ecc, f_bar).
    """
    
    # Termination events -------------------------------------------------
    
    def circularized(tau, state, *args):
        return state[1] - e_term
    circularized.terminal = True

    def eccentricity_minimum(*args):
        _, de_dtau = GergelyODEs(*args)
        return de_dtau
    eccentricity_minimum.terminal = True
    eccentricity_minimum.direction = 1

    # Solve the ODEs -----------------------------------------------------
    
    # Initial conditions
    omega_orb_bar = np.pi * f_start_bar 
    a_start_bar   = 1.0 / omega_orb_bar**(2/3) # Keplerian approximation
    
    # Integrate using solve_ivp
    sol = solve_ivp(
        fun    = GergelyODEs,
        t_span = [0, max_tau],
        y0     = [a_start_bar, ecc_start],
        args   = (nu, q, S1z_bar, S2z_bar, e_term),
        events = [circularized, eccentricity_minimum],
        method = 'RK45',       
        dense_output = True,
        **kwargs
    )
    
    # Extract results ----------------------------------------------------

    tau_arr   = sol.t
    a_bar_arr = sol.y[0]
    ecc_arr   = sol.y[1]

    last_valid_idx = np.argmin(ecc_arr)

    tau_arr   = tau_arr[: last_valid_idx+1]
    a_bar_arr = a_bar_arr[: last_valid_idx+1]
    ecc_arr   = ecc_arr[: last_valid_idx+1]

    # Again Keplerian approximation
    f_bar_arr = (1.0 / np.pi) * (a_bar_arr**(-3/2))

    return tau_arr, a_bar_arr, ecc_arr, f_bar_arr


def gergely_f_start_bar(
    nu : float,
    q : float, 
    S1z_bar : float, 
    S2z_bar : float, 
    ecc_start : float, 
    target_tau : float, 
    f_min_bar : float = 1e-7,
    f_max_bar : float = 1e-2,
    **kwargs
) -> float :
    
    """
    Finds the starting dimensionless frequency for a target dimensionless duration.
    Uses the Gergely ODEs. This method looses accuracy at high eccentricities.
    Use the teobresumsfits routines for more accurate phenomenological fits.

    Parameters
    ----------
    nu : float
        Dimensionless symmetric mass ratio, i.e., m1*m2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., chi1 / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., chi2 / (m2 / m_tot)^2
    ecc_start : float
        Starting eccentricity.
    target_tau : float
        Target dimensionless duration of the waveform.
    f_min_bar : float, optional
        Minimum dimensionless frequency to consider in the search.
        Default is 1e-7.
    f_max_bar : float, optional
        Maximum dimensionless frequency to consider in the search.
        Default is 1e-2.
    **kwargs
        Additional keyword arguments to pass to the ODE solver.
    
    Returns
    -------
    float
        The starting dimensionless frequency f_start_bar that achieves the target duration.
    """
    
    def duration_error(f_guess_bar) :
        tau, _, _, _ = GergelyEvolve(
            nu=nu, 
            q=q, 
            S1z_bar=S1z_bar, 
            S2z_bar=S2z_bar, 
            f_start_bar=f_guess_bar, 
            ecc_start=ecc_start,
            max_tau=target_tau*1e1,
            **kwargs
        )
        return (tau[-1] - target_tau)

    while True :

        try : 
            result = root_scalar(
                duration_error,
                bracket=[f_min_bar, f_max_bar], 
                method='brentq',
                xtol=1e-10
            )
            return result.root
        
        except ValueError :
            f_max_bar *= 2.0
            f_min_bar *= 0.5
            continue

        except Exception as e :
            print(f"Initial conditions failed : {e}")
            return 0


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

    orbavg_poly = log_f_orbavg_bar_interpolator.poly(log_sig_tau, ecc_start, chi_eff)
    all_roots   = (orbavg_poly - log_fbar).roots()
    real_roots  = all_roots[np.isreal(all_roots)].real
    if len(real_roots) == 0 : raise ValueError("No real roots found for the given parameters.")
    interp_low  = log_f_orbavg_bar_interpolator.poly_interp_low
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
    f_peri_bar = 10 ** log_f_peri_bar_interpolator.poly(log_sig_tau, ecc_start, chi_eff)(tau)
    f_apos_bar = 10 ** log_f_apos_bar_interpolator.poly(log_sig_tau, ecc_start, chi_eff)(tau)
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
    ecc_bracket = [log_f_orbavg_bar_interpolator.p_min[1], 
                   log_f_orbavg_bar_interpolator.p_max[1]]
    
    try : ecc_start = root_scalar(ecc_error, bracket=ecc_bracket, method='brentq').root
    except ValueError as e : raise ValueError(f"Initial conditions failed: {e}")
    
    required_log_fbar_start = log_f_orbavg_bar_interpolator(
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
        measured_log_fbar_start = log_f_calib_bar_interpolator(
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