"""
Analytic Functions Related to Eccentricity
"""

# [NOTE] :
# Throughout this file, the "bar" suffix denotes dimensionless quantities.
# Except for the `ecc_from_envelop_freqs` function, the rest of the functions 
# are exclusively dimensionless. To get dimensionless quantities from physical 
# quantities and vice versa, use the conversion methods under `core.sources.BinarySystem`.


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar


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
        Dimensionless symmetric mass ratio, i.e., m1*m2 / m_tot^2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., spin1z / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., spin2z / (m2 / m_tot)^2
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
        Dimensionless symmetric mass ratio, i.e., m1*m2 / m_tot^2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., spin1z / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., spin2z / (m2 / m_tot)^2
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
        Dimensionless symmetric mass ratio, i.e., m1*m2 / m_tot^2
    q : float
        Mass ratio m1/m2
    S1z_bar : float
        Dimensionless angular momentum of the primary along the orbital angular momentum.
        i.e., spin1z / (m1 / m_tot)^2
    S2z_bar : float
        Dimensionless angular momentum of the secondary along the orbital angular momentum.
        i.e., spin2z / (m2 / m_tot)^2
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