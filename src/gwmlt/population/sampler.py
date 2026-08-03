"""
Sample Binary Parameters
"""

import random
import numpy as np
import pandas as pd

from pathlib import Path
from joblib import Memory
from scipy.interpolate import RectBivariateSpline

import bilby
import gwpopulation

from gwmlt.config import config
from gwmlt.utils.cosmology import (
    merger_rate_pdf_interpolator, 
    luminosity_distance
)

# q|m_src PDF interpolator cache
memory = Memory(str(Path(__file__).resolve().parent/"__pycache__"), verbose=0) 


def population_sampler(
    num_samples : int, 
    seed : int,
    grid_points : int = 1000,
    ecc_range : tuple[int, int] = (0.05, 0.50),
    use_q_pdf_interpolator : bool = False
) -> pd.DataFrame :

    """
    Sample binary parameters from GWTC data.

    Parameters
    ----------
    num_samples : int
        Number of samples to generate.
    seed : int
        Random seed
    grid_points
        Number of grid points to use for the interpolation of various PDFs.
        Must not be same as `num_samples`. Defaults to 1000.
    ecc_range : tuple[int, int]
        Range of eccentricities to sample from. Defaults to (0.05, 0.50).
    use_q_pdf_interpolator : bool
        Use an bivariate spline interpolator for obtaining the `q|m1_src` PDF. 
        This is much faster than constructing the distributions from scratch.
        The first time this function is executed, a cache is created, which may take a few seconds.
        Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the sampled parameters.
    """

    # Part 0. Setup -------------------------------------------------------------------------------

    assert num_samples != grid_points, \
        "Number of samples can not be equal to number of grid points, since it leads to unintended broadcasting"

    ecc_min,          ecc_max           = ecc_range
    m1_src_min,       m1_src_max        = config.population.m1_src_range
    m_det_low_cutoff, m_det_high_cutoff = config.population.m_det_cutoff_range
    z_min,            z_max             = config.population.source_redshift_range
    ml_min,           ml_max            = config.population.lens_mass_range
    yl_min,           yl_max            = config.population.impact_parameter_range
    t_start,          t_end             = config.population.geocent_time_range


    # Part 1. Loading the data --------------------------------------------------------------------

    # Read the population file path from the configuration
    population_file_path = config.population.population_file_path
    pop_result = bilby.core.result.read_in_result(filename=population_file_path)

    # Extract the MAP parameters from the population result
    log_likelihood = pop_result.posterior["log_likelihood"]
    log_prior      = pop_result.posterior["log_prior"]
    map_index      = np.argmax(log_prior + log_likelihood)
    pop_params     = pop_result.posterior.iloc[map_index].to_dict()

    # Convert spin distribution parameters to beta function parameters
    pop_params, _  = gwpopulation.conversions.convert_to_beta_parameters(pop_params)


    # Part 2. Distribution PDFs -------------------------------------------------------------------

    # 2a. Redshift distribution model ------------------------------------

    z_interp_grid = np.linspace(z_min, z_max, grid_points)
    z_pdf = merger_rate_pdf_interpolator(z_min, z_max)(z_interp_grid)

    # 2b. Mass distribution model ----------------------------------------

    mass_model = gwpopulation.models.mass.SinglePeakSmoothedMassDistribution(cache=False)

    m1_src_interp_grid = np.geomspace(m1_src_min, m1_src_max, grid_points)
    q_interp_grid      = np.linspace(1.0, 0.0, grid_points, endpoint=False)[::-1]

    # Primary mass (source frame) PDF

    m1_src_pdf = mass_model.p_m1(
        dataset = dict(mass_1 = m1_src_interp_grid),
        **{
            "alpha"   : pop_params["alpha"],
            "mmin"    : pop_params["mmin"],
            "mmax"    : pop_params["mmax"],
            "lam"     : pop_params["lam"],
            "mpp"     : pop_params["mpp"],
            "sigpp"   : pop_params["sigpp"],
            "delta_m" : pop_params["delta_m"],
        }
    )

    # Mass ratio PDF for given primary mass (source frame)

    def q_pdf_cnstrc(m1_src) :
        dataset = dict(
            mass_1     = m1_src, 
            mass_ratio = q_interp_grid
        )
        mass_model = gwpopulation.models.mass.SinglePeakSmoothedMassDistribution(cache=False)
        return mass_model.p_q(
            dataset = dataset,
            beta    = pop_params["beta"],
            mmin    = pop_params["mmin"],
            delta_m = pop_params["delta_m"],
        ).copy()

    # Fast q|m1_src pdf interpolator

    if use_q_pdf_interpolator :

        @memory.cache
        def _get_cached_q_pdf_grid(config, grid_points) :
            return np.array([q_pdf_cnstrc(m1_src) for m1_src in m1_src_interp_grid])

        _q_pdf_spline = RectBivariateSpline(m1_src_interp_grid, q_interp_grid,
            _get_cached_q_pdf_grid(config, grid_points))

        def q_pdf_spline(m1_src) : return _q_pdf_spline(m1_src, q_interp_grid)[0]
        q_pdf = q_pdf_spline

    else : q_pdf = q_pdf_cnstrc 

    # 2c. Spin tilt model ------------------------------------------------

    cos_tilt_interp_grid = np.linspace(-1.0, 1.0, grid_points)

    # Based on `gwpopulation.models.spin.iid_spin_orientation_gaussian_isotropic`
    aligned_pdf = gwpopulation.utils.truncnorm(
        cos_tilt_interp_grid, mu=1.0, sigma=pop_params["sigma_spin"], high=1.0, low=-1.0
    )
    isotropic_pdf = 0.5 * np.ones_like(cos_tilt_interp_grid)  # uniform density on [-1,1]


    # Part 3. Bilby priors ------------------------------------------------------------------------

    # 3a. Conversion function --------------------------------------------

    def evaluate_derived_params(parameters : dict) -> dict :

        super_params = parameters.copy()
        
        # Calculate source-frame secondary mass
        super_params['m2_src'] = super_params['m1_src'] * super_params['mass_ratio']

        # Calculate detector-frame masses
        super_params['mass_1'] = super_params['m1_src'] * (1.0 + super_params['z'])
        super_params['mass_2'] = super_params['m2_src'] * (1.0 + super_params['z'])
        
        # Calculate luminosity distance 
        super_params['luminosity_distance'] = np.vectorize(luminosity_distance)(super_params['z'])

        # Calculate tilts
        super_params['tilt_1'] = np.arccos(super_params['cos_tilt_1'])
        super_params['tilt_2'] = np.arccos(super_params['cos_tilt_2'])
        
        return super_params

    # 3b. Initialize the conditional prior dictionary --------------------

    priors = bilby.core.prior.ConditionalPriorDict(conversion_function=evaluate_derived_params)

    # 3c. Independent base priors ----------------------------------------

    # Source redshift
    priors['z'] = bilby.core.prior.interpolated.Interped(
        z_interp_grid, z_pdf, 
        minimum=z_min, maximum=z_max, 
        name='z'
    )

    # Source frame primary mass
    priors['m1_src'] = bilby.core.prior.interpolated.Interped(
        m1_src_interp_grid, m1_src_pdf, 
        minimum=m1_src_min, maximum=m1_src_max, 
        name='m1_src'
    )

    # Isotropic/aligned spin selection
    priors['u_spin'] = bilby.core.prior.Uniform(minimum=0.0, maximum=1.0, name='u_spin')

    # 3d. Conditional priors ---------------------------------------------

    # Lens Redshift
    def condition_func_zl(reference_params, z):
        return dict(
            minimum = 0.0,
            maximum = z
        )

    priors['zl'] = bilby.core.prior.ConditionalUniform(
        condition_func = condition_func_zl,
        minimum = 0.0,
        maximum = z_max,
        name = 'zl'
    )

    # Mass ratio
    def condition_func_q(reference_params, m1_src):
        return dict(
            xx = q_interp_grid,
            yy = q_pdf(m1_src)
        )
    
    priors['mass_ratio'] = bilby.core.prior.ConditionalInterped(
        condition_func = condition_func_q,
        xx = q_interp_grid, yy = np.ones_like(q_interp_grid),
        minimum = np.min(q_interp_grid), maximum = 1.0,
        name = 'mass_ratio'
    )

    # Spin tilts
    def condition_func_tilt(reference_params, u_spin):
        return dict(
            xx = cos_tilt_interp_grid, 
            yy = aligned_pdf if u_spin<pop_params["xi_spin"] else isotropic_pdf
        )

    priors['cos_tilt_1'] = bilby.core.prior.ConditionalInterped(
        condition_func = condition_func_tilt,
        xx = cos_tilt_interp_grid, yy = np.ones_like(cos_tilt_interp_grid),
        minimum = -1.0, maximum = 1.0, 
        name='cos_tilt_1'
    )
    priors['cos_tilt_2'] = bilby.core.prior.ConditionalInterped(
        condition_func = condition_func_tilt,
        xx = cos_tilt_interp_grid, yy = np.ones_like(cos_tilt_interp_grid),
        minimum = -1.0, maximum = 1.0, 
        name = 'cos_tilt_2'
    )

    # 3e. Spin priors ----------------------------------------------------

    # Magnitudes (beta distribution)
    priors['a_1'] = bilby.core.prior.Beta(
        alpha=pop_params["alpha_chi"], beta=pop_params["beta_chi"], name='a_1'
    )
    priors['a_2'] = bilby.core.prior.Beta(
        alpha=pop_params["alpha_chi"], beta=pop_params["beta_chi"], name='a_2'
    )

    # In-plane spin azimuthal angles
    priors['phi_12'] = bilby.core.prior.Uniform(
        minimum=0.0, maximum=2*np.pi, 
        boundary="periodic", name='phi_12'
    )
    priors['phi_jl'] = bilby.core.prior.Uniform(
        minimum=0.0, maximum=2*np.pi, 
        boundary="periodic", name='phi_jl'
    )

    # 3f. Geometric priors -----------------------------------------------

    priors['ra']       = bilby.core.prior.Uniform (name='ra',       minimum=0.0, maximum=2*np.pi, boundary="periodic")
    priors['dec']      = bilby.core.prior.Cosine  (name='dec',                                                       )
    priors['phase']    = bilby.core.prior.Uniform (name='phase',    minimum=0.0, maximum=2*np.pi, boundary="periodic")
    priors['theta_jn'] = bilby.core.prior.Sine    (name='theta_jn',                                                  )
    priors['psi']      = bilby.core.prior.Uniform (name='psi',      minimum=0.0, maximum=1*np.pi, boundary="periodic")

    # 3g. Lensing priors -------------------------------------------------

    priors['ml'] = bilby.core.prior.LogUniform (name='ml', minimum=ml_min, maximum=ml_max           )
    priors['yl'] = bilby.core.prior.PowerLaw   (name='yl', minimum=yl_min, maximum=yl_max, alpha=1.0)

    # 3h. Eccentricity priors --------------------------------------------

    priors['ecc']     = bilby.core.prior.Uniform(name='ecc',     minimum=ecc_min, maximum=ecc_max)
    priors['anomaly'] = bilby.core.prior.Uniform(name='anomaly', minimum=0.0,     maximum=2*np.pi, boundary="periodic")

    # 3i. Trigger time ---------------------------------------------------

    priors['geocent_time'] = bilby.core.prior.Uniform(name='geocent_time', minimum=t_start, maximum=t_end)

    # 3j. Constraints ----------------------------------------------------

    # Component masses in the detector frame must be above 5.0
    priors['mass_1'] = bilby.core.prior.Constraint(minimum=m_det_low_cutoff, maximum=m_det_high_cutoff, name='mass_1')
    priors['mass_2'] = bilby.core.prior.Constraint(minimum=m_det_low_cutoff, maximum=m_det_high_cutoff, name='mass_2')


    # Part 4. Sampling ----------------------------------------------------------------------------

    bilby.core.utils.random.seed(seed)

    samples = priors.sample(num_samples)
    del samples['u_spin']

    full_samples = evaluate_derived_params(samples)
    return pd.DataFrame(full_samples)