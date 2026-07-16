"""
Interpolators for Teobresums Evolution Fits
"""

from warnings import warn

import numpy as np
from scipy.interpolate import RBFInterpolator

from ..config import config


# The Interpolator Class --------------------------------------------------------------------------

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
        Initialize an R3 -> R1 RBF Interpolator.

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
        log_waveform_duration_dimensionless, ecc_start, chi_eff.
        """

        # Required for brentq
        q_points = [np.array(i).flatten() for i in query_points]
        q_points = [i.item() if len(i)==1 else i for i in q_points]

        q_points = np.atleast_2d(q_points)

        q_min = q_points.min(axis=0)
        q_max = q_points.max(axis=0)

        if np.any(q_min < self.p_min) or np.any(q_max > self.p_max):
            warn("Query points are outside the range of the training data.\n" 
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
        log_waveform_duration_dimensionless, ecc_start, chi_eff.
        """

        if self.anchor_points == 1 :
            raise ValueError("Polynomial fit is not available for scalar data.")

        anchor_ys = self.__call__([x_val, y_val, z_val])
        anchor_xs = np.linspace(self.poly_interp_low, self.poly_interp_high or x_val, self.anchor_points)

        return np.polynomial.Polynomial.fit(anchor_xs, anchor_ys, deg=deg)


# Initializing the Interpolators ------------------------------------------------------------------

with np.load(config.teobresums_fits_path) as f : data = dict(f)

# min_log_tau_interp sets the lower bound for datapoints and polynomial interpolation.
# All logarithms are base 10.

log_f_peri_bar_interpolator = AnisotropicTEFPhenom(
    x_vals          = data['log_waveform_duration_dimensionless'],
    y_vals          = data['ecc_start'],
    z_vals          = data['chi_eff'],
    fit_vals        = data['log_f_peri_anchors_dimensionless'],
    poly_interp_low = data['min_log_tau_interp'],
)
"""Interpolator for log of dimensionless periastron frequencies.
The three coordinates are (respectively) : 
log_waveform_duration_dimensionless, ecc_start, chi_eff."""

log_f_apos_bar_interpolator = AnisotropicTEFPhenom(
    x_vals          = data['log_waveform_duration_dimensionless'],
    y_vals          = data['ecc_start'],
    z_vals          = data['chi_eff'],
    fit_vals        = data['log_f_apos_anchors_dimensionless'],
    poly_interp_low = data['min_log_tau_interp'],
)
"""Interpolator for log of dimensionless apastron frequencies.
The three coordinates are (respectively) : 
log_waveform_duration_dimensionless, ecc_start, chi_eff."""

log_f_orbavg_bar_interpolator = AnisotropicTEFPhenom(
    x_vals          = data['log_waveform_duration_dimensionless'],
    y_vals          = data['ecc_start'],
    z_vals          = data['chi_eff'],
    fit_vals        = data['log_f_orbavg_anchors_dimensionless'],
    poly_interp_low = data['min_log_tau_interp'],
)
"""Interpolator for log of dimensionless orbital average frequencies.
The three coordinates are (respectively) : 
log_waveform_duration_dimensionless, ecc_start, chi_eff."""

log_f_calib_bar_interpolator = AnisotropicTEFPhenom(
    x_vals          = data['log_requested_f_start_dimensionless'],
    y_vals          = data['ecc_start'],
    z_vals          = data['chi_eff'],
    fit_vals        = data['log_measured_f_start_dimensionless'],
    poly_interp_low = data['min_log_tau_interp'],
)
"""Interpolator calibrating log of dimensionless starting frequencies.
This is used to correct f_start offsets in TEOBResums.
The three coordinates are (respectively) : 
log_requested_f_start_dimensionless, ecc_start, chi_eff.
The output is the log_measured_f_start_dimensionless."""

del data