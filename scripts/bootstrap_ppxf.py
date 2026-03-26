"""
Bootstrap error estimation for ppxf velocity dispersion fits.

Implements a hybrid wild bootstrap with local residual scaling to produce
empirical error distributions for V and sigma from ppxf fits of the
AGEL0206 deflector galaxy integrated spectrum.

Following ppxf_example_population_bootstrap.py (Cappellari) with enhancements:
  - Local residual scaling: rolling-window std captures wavelength-dependent noise
  - Rademacher wild bootstrap: preserves heteroscedastic noise structure
  - Degree sweep: errors as a function of additive polynomial degree
  - Multi-template: runs for FSPS, EMILES, XSL template libraries

Usage:
    # As import:
    from scripts.bootstrap_ppxf import run_bootstrap, setup_ppxf_inputs

    # As CLI:
    python scripts/bootstrap_ppxf.py --sps_name fsps --n_bootstrap 500
    python scripts/bootstrap_ppxf.py --sps_name all
    python scripts/bootstrap_ppxf.py --degrees 4,10,16,20 --n_bootstrap 50  # quick test

Requires: numpy, astropy, ppxf, scipy, tqdm
"""

import numpy as np
from pathlib import Path
from importlib import resources
from urllib import request
from time import perf_counter as clock
from tqdm import tqdm
import argparse
import os

from astropy.io import fits
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from scipy.ndimage import uniform_filter1d


# =============================================================================
# Constants
# =============================================================================
C_KMS = 299792.458  # speed of light in km/s

# Default AGEL0206 parameters
DEFAULT_IFU = 'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits'
DEFAULT_Z = 0.67511
DEFAULT_LAM_OBS_RANGE = (6500.0, 7500.0)
DEFAULT_LAM_RANGE_TEMP = (3500, 5000)

# Deflector and noise spaxel regions
DEFLECTOR_SLICE = (slice(45, 64), slice(45, 55))  # cube[:, 45:64, 45:55]
NOISE_SLICE = (slice(28, 40), slice(45, 70))       # cube[:, 28:40, 45:70]


# =============================================================================
# Setup functions
# =============================================================================

def setup_ppxf_inputs(ifu_file, sps_name='fsps', z=DEFAULT_Z,
                      lam_obs_range=DEFAULT_LAM_OBS_RANGE,
                      lam_range_temp=DEFAULT_LAM_RANGE_TEMP):
    """
    Recreate all ppxf inputs from the raw IFU cube.

    Mirrors the setup in 01_streamlined_IFU_ppxf.ipynb (cells 8e8bae4a, a30f240b).

    Parameters
    ----------
    ifu_file : str
        Path to the KCWI FITS cube.
    sps_name : str
        Template library: 'fsps', 'emiles', or 'xsl'.
    z : float
        Deflector redshift.
    lam_obs_range : tuple
        Observed wavelength range to fit (Angstroms).
    lam_range_temp : tuple
        Rest-frame wavelength range for template loading.

    Returns
    -------
    dict with keys:
        'galaxy'       : normalized galaxy spectrum (1001 pixels)
        'noise'        : normalized noise spectrum (1001 pixels)
        'velscale'     : velocity scale in km/s per pixel
        'start'        : starting guess [V, sigma] in km/s
        'goodpixels'   : indices of pixels included in fit
        'lam_gal_rest' : rest-frame wavelength array
        'sps'          : sps_lib object with templates
        'lam_temp'     : template wavelength array
        'sps_name'     : template library name
    """
    print(f"Setting up ppxf inputs for {sps_name}...")

    # Load cube and header
    with fits.open(ifu_file) as hdul:
        hdr = hdul[0].header
        cube = np.asarray(hdul[0].data, dtype=float)

    # Build wavelength array from WCS
    crval = hdr['CRVAL3']
    cdelt = hdr['CD3_3']
    npix = hdr['NAXIS3']
    crpix = hdr.get('CRPIX3', 1.0)
    pix = np.arange(npix)
    lam = crval + cdelt * (pix + 1 - crpix)

    # Extract integrated deflector and noise spectra
    flux_int = np.average(cube[:, DEFLECTOR_SLICE[0], DEFLECTOR_SLICE[1]], axis=(1, 2))
    noise_int = np.std(cube[:, NOISE_SLICE[0], NOISE_SLICE[1]], axis=(1, 2))

    # Rebuild wavelength and trim to fitting range
    lam_int = crval + cdelt * (pix + 1 - crpix)
    mask_wl = (lam_int >= lam_obs_range[0]) & (lam_int <= lam_obs_range[1])
    lam_int = lam_int[mask_wl]
    flux_int = flux_int[mask_wl]
    noise_int = noise_int[mask_wl]

    # Rebin to logarithmic wavelength spacing
    log_lam = np.log(lam_int)
    d_log_lam = (log_lam[-1] - log_lam[0]) / (len(lam_int) - 1)
    log_lam_new = np.arange(log_lam[0], log_lam[-1] + d_log_lam, d_log_lam)
    flux_int = np.interp(log_lam_new, log_lam, flux_int)
    noise_int = np.interp(log_lam_new, log_lam, noise_int)
    lam_int = np.exp(log_lam_new)

    # Final wavelength mask after rebinning
    mask_fit = (lam_int >= 6000.0) & (lam_int <= 7500.0)
    lam_int = lam_int[mask_fit]
    flux_int = flux_int[mask_fit]
    noise_int = noise_int[mask_fit]

    # Normalize spectrum
    median_flux = np.median(flux_int)
    galaxy = flux_int / median_flux
    noise_norm = np.sqrt(noise_int**2) / median_flux

    # Convert vacuum to air wavelengths
    lam_gal = np.copy(lam_int)
    lam_gal *= np.median(util.vac_to_air(lam_gal) / lam_gal)

    # Velocity scale
    d_ln_lam_gal = (np.log(lam_gal)[-1] - np.log(lam_gal)[0]) / (np.log(lam_gal).size - 1)
    velscale = C_KMS * d_ln_lam_gal

    # Load SPS templates
    ppxf_dir = resources.files('ppxf')
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = ppxf_dir / 'sps_models' / basename
    if not filename.is_file():
        url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
        request.urlretrieve(url, filename)

    # Instrumental resolution from FITS header
    dlam_gal = np.gradient(lam_gal)
    wdisp = hdr['DISPSCAL']
    fwhm_gal = 2.355 * wdisp * dlam_gal

    # Convert to rest frame
    lam_gal_rest = lam_gal / (1 + z)
    fwhm_gal_rest = fwhm_gal / (1 + z)

    # Load templates
    fwhm_gal_dict = {"lam": lam_gal_rest, "fwhm": fwhm_gal_rest}
    sps = lib.sps_lib(filename, velscale, fwhm_gal_dict, lam_range=list(lam_range_temp))
    goodpixels = util.determine_goodpixels(np.log(lam_gal_rest), list(lam_range_temp))

    print(f"  Templates: {sps_name} ({len(sps.templates.T)} templates)")
    print(f"  Galaxy: {len(galaxy)} pixels, Good pixels: {len(goodpixels)}")
    print(f"  Velscale: {velscale:.2f} km/s/pix")

    return {
        'galaxy': galaxy,
        'noise': noise_norm,
        'velscale': velscale,
        'start': [0., 300.],
        'goodpixels': goodpixels,
        'lam_gal_rest': lam_gal_rest,
        'sps': sps,
        'lam_temp': sps.lam_temp,
        'sps_name': sps_name,
        'z': z,
    }


# =============================================================================
# Bootstrap functions
# =============================================================================

def compute_local_residual_scaling(residuals, window=75):
    """
    Compute local-to-global residual scatter ratio for wild bootstrap scaling.

    Uses a rolling window to estimate wavelength-dependent noise structure,
    then scales each pixel's residual by the local/global std ratio.

    Parameters
    ----------
    residuals : array, shape (n_pixels,)
        Residuals = galaxy - best_fit for a given degree.
    window : int
        Rolling window size in pixels for local std computation.

    Returns
    -------
    scale_factor : array, shape (n_pixels,)
        Ratio of local_std / global_std at each pixel.
        Clipped to [0.2, 5.0] to avoid extreme outliers.
    """
    global_std = np.std(residuals)
    if global_std == 0:
        return np.ones_like(residuals)

    # Local variance via rolling window: Var = E[x^2] - E[x]^2
    local_mean = uniform_filter1d(residuals, size=window, mode='reflect')
    local_mean_sq = uniform_filter1d(residuals**2, size=window, mode='reflect')
    local_var = local_mean_sq - local_mean**2
    local_var = np.maximum(local_var, 0)  # numerical safety
    local_std = np.sqrt(local_var)

    scale_factor = local_std / global_std
    scale_factor = np.clip(scale_factor, 0.2, 5.0)

    return scale_factor


def run_bootstrap_single_degree(ppxf_inputs, degree, best_fit_spectrum,
                                 n_bootstrap=500, window=75, seed=None):
    """
    Run hybrid wild bootstrap for a single polynomial degree.

    Parameters
    ----------
    ppxf_inputs : dict
        Output of setup_ppxf_inputs().
    degree : int
        Additive polynomial degree.
    best_fit_spectrum : array, shape (n_pixels,)
        Best-fit spectrum from original ppxf run at this degree.
    n_bootstrap : int
        Number of bootstrap iterations.
    window : int
        Rolling window size for local residual scaling.
    seed : int or None
        RNG seed for reproducibility.

    Returns
    -------
    dict with keys:
        'V_samples'     : array, shape (n_bootstrap,)
        'sigma_samples' : array, shape (n_bootstrap,)
        'chi2_samples'  : array, shape (n_bootstrap,)
        'n_failed'      : int, number of failed ppxf fits
    """
    galaxy = ppxf_inputs['galaxy']
    noise = ppxf_inputs['noise']
    velscale = ppxf_inputs['velscale']
    start = ppxf_inputs['start']
    goodpixels = ppxf_inputs['goodpixels']
    lam_gal_rest = ppxf_inputs['lam_gal_rest']
    sps = ppxf_inputs['sps']
    lam_temp = ppxf_inputs['lam_temp']

    # Compute residuals and local scaling
    residuals = galaxy - best_fit_spectrum
    scale_factor = compute_local_residual_scaling(residuals, window=window)
    scaled_residuals = residuals * scale_factor

    # Input redshift for converting V back to z
    z0 = ppxf_inputs.get('z', DEFAULT_Z)

    # Initialize storage
    V_samples = np.full(n_bootstrap, np.nan)
    sigma_samples = np.full(n_bootstrap, np.nan)
    chi2_samples = np.full(n_bootstrap, np.nan)
    z_samples = np.full(n_bootstrap, np.nan)
    n_failed = 0

    rng = np.random.default_rng(seed)

    for i in range(n_bootstrap):
        # Rademacher wild bootstrap: flip sign of each scaled residual
        signs = rng.choice([-1, 1], size=len(galaxy))
        galaxy_boot = best_fit_spectrum + scaled_residuals * signs

        try:
            pp = ppxf(sps.templates, galaxy_boot, noise, velscale, start,
                      goodpixels=goodpixels, plot=False, moments=2, trig=False,
                      degree=degree, lam=lam_gal_rest, lam_temp=lam_temp, mdegree=0)
            V_samples[i] = pp.sol[0]
            sigma_samples[i] = pp.sol[1]
            chi2_samples[i] = pp.chi2
            # Convert fitted velocity to redshift (eq. 5c, Cappellari 2023)
            z_samples[i] = (1 + z0) * np.exp(pp.sol[0] / C_KMS) - 1
        except Exception:
            n_failed += 1

    return {
        'V_samples': V_samples,
        'sigma_samples': sigma_samples,
        'chi2_samples': chi2_samples,
        'z_samples': z_samples,
        'n_failed': n_failed,
    }


def run_bootstrap(ifu_file, sps_name='fsps', results_dir='results',
                  degrees=None, n_bootstrap=500, window=75, seed=42,
                  save=True, z=None, save_suffix=None):
    """
    Full bootstrap pipeline for all degrees and one template library.

    Parameters
    ----------
    ifu_file : str
        Path to the KCWI FITS cube.
    sps_name : str
        Template library name ('fsps', 'emiles', or 'xsl').
    results_dir : str
        Directory containing saved ppxf results and for output.
    degrees : array or None
        Polynomial degrees to bootstrap. None = all from saved results.
    n_bootstrap : int
        Number of bootstrap iterations per degree.
    window : int
        Rolling window for local residual scaling.
    seed : int
        Base RNG seed (each degree gets seed + degree_index).
    save : bool
        Whether to save results to .npz.
    z : float or None
        Deflector redshift. None = use default (0.67511).
    save_suffix : str or None
        Extra suffix for output filename, e.g. 'z067564'.
        Output: ppxf_bootstrap_errors_{sps_name}_{save_suffix}.npz

    Returns
    -------
    dict with keys:
        'sigma_bootstrap'    : array, shape (n_degrees, n_bootstrap) — full distribution
        'V_bootstrap'        : array, shape (n_degrees, n_bootstrap) — full distribution
        'chi2_bootstrap'     : array, shape (n_degrees, n_bootstrap)
        'sigma_boot_err_lo'  : array, shape (n_degrees,) — median - 16th percentile
        'sigma_boot_err_hi'  : array, shape (n_degrees,) — 84th percentile - median
        'sigma_p16/p50/p84'  : arrays, shape (n_degrees,) — percentiles
        'V_boot_err_lo/hi'   : arrays, shape (n_degrees,) — asymmetric V errors
        'V_p16/p50/p84'      : arrays, shape (n_degrees,)
        'sigma_original'     : array, shape (n_degrees,) — original fit values
        'V_original'         : array, shape (n_degrees,)
        'chi2_original'      : array, shape (n_degrees,)
        'degrees'            : array
        'sps_name'           : str
        'n_bootstrap'        : int
        'n_failed'           : array, shape (n_degrees,)
    """
    # Setup ppxf inputs
    setup_kwargs = {'ifu_file': ifu_file, 'sps_name': sps_name}
    if z is not None:
        setup_kwargs['z'] = z
    ppxf_inputs = setup_ppxf_inputs(**setup_kwargs)

    # Load saved best-fit results
    # Try template-specific file first, then generic
    results_path = os.path.join(results_dir, f'ppxf_integrated_spectrum_results_{sps_name}.npz')
    if not os.path.exists(results_path):
        results_path = os.path.join(results_dir, 'ppxf_integrated_spectrum_results.npz')
    saved = np.load(results_path, allow_pickle=True)

    all_degrees = saved['degrees']
    all_best_fit = saved['best_fit']
    all_vel_dis = saved['vel_dis']
    all_mean_vel = saved['mean_vel']
    all_fit_chi2 = saved['fit_chi2']

    if degrees is None:
        degrees = all_degrees
    degree_indices = [np.where(all_degrees == d)[0][0] for d in degrees]

    print(f"\n{'=' * 60}")
    print(f"Bootstrap error estimation: {sps_name}")
    print(f"  Redshift: {ppxf_inputs['z']:.5f}")
    print(f"  Degrees: {degrees}")
    print(f"  N_bootstrap: {n_bootstrap}")
    print(f"  Window: {window} pixels")
    print(f"  Seed: {seed}")
    print(f"{'=' * 60}")

    n_deg = len(degrees)
    V_bootstrap = np.full((n_deg, n_bootstrap), np.nan)
    sigma_bootstrap = np.full((n_deg, n_bootstrap), np.nan)
    chi2_bootstrap = np.full((n_deg, n_bootstrap), np.nan)
    z_bootstrap = np.full((n_deg, n_bootstrap), np.nan)
    n_failed = np.zeros(n_deg, dtype=int)

    t0 = clock()

    for j, (deg, idx) in enumerate(tqdm(zip(degrees, degree_indices),
                                         total=n_deg,
                                         desc=f"Bootstrap {sps_name}")):
        result = run_bootstrap_single_degree(
            ppxf_inputs, degree=deg,
            best_fit_spectrum=all_best_fit[idx],
            n_bootstrap=n_bootstrap,
            window=window,
            seed=seed + j,
        )
        V_bootstrap[j] = result['V_samples']
        sigma_bootstrap[j] = result['sigma_samples']
        chi2_bootstrap[j] = result['chi2_samples']
        z_bootstrap[j] = result['z_samples']
        n_failed[j] = result['n_failed']

    elapsed = clock() - t0
    print(f"\nCompleted in {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Compute bootstrap errors from 16th/84th percentiles (asymmetric)
    sigma_p16 = np.nanpercentile(sigma_bootstrap, 16, axis=1)
    sigma_p50 = np.nanpercentile(sigma_bootstrap, 50, axis=1)
    sigma_p84 = np.nanpercentile(sigma_bootstrap, 84, axis=1)
    sigma_boot_err_lo = sigma_p50 - sigma_p16   # lower 1-sigma
    sigma_boot_err_hi = sigma_p84 - sigma_p50   # upper 1-sigma

    V_p16 = np.nanpercentile(V_bootstrap, 16, axis=1)
    V_p50 = np.nanpercentile(V_bootstrap, 50, axis=1)
    V_p84 = np.nanpercentile(V_bootstrap, 84, axis=1)
    V_boot_err_lo = V_p50 - V_p16
    V_boot_err_hi = V_p84 - V_p50

    z_p16 = np.nanpercentile(z_bootstrap, 16, axis=1)
    z_p50 = np.nanpercentile(z_bootstrap, 50, axis=1)
    z_p84 = np.nanpercentile(z_bootstrap, 84, axis=1)

    # Original values at the selected degrees
    V_original = all_mean_vel[degree_indices]
    sigma_original = all_vel_dis[degree_indices]
    chi2_original = all_fit_chi2[degree_indices]

    # Report
    print(f"\nResults for {sps_name}:")
    print(f"{'Deg':>4} {'sigma':>7} {'-err':>6} {'+err':>6} {'V':>7} {'-err':>6} {'+err':>6} {'z_med':>8} {'-dz':>7} {'+dz':>7} {'fail':>5}")
    print("-" * 80)
    formal_err = saved['error_vdis'][degree_indices]
    for j, deg in enumerate(degrees):
        print(f"{deg:4d} {sigma_original[j]:7.1f} {sigma_boot_err_lo[j]:6.1f} "
              f"{sigma_boot_err_hi[j]:6.1f} "
              f"{V_original[j]:7.1f} {V_boot_err_lo[j]:6.1f} {V_boot_err_hi[j]:6.1f} "
              f"{z_p50[j]:8.5f} {z_p50[j]-z_p16[j]:7.5f} {z_p84[j]-z_p50[j]:7.5f} "
              f"{n_failed[j]:5d}")

    if np.any(n_failed > 0):
        print(f"\nWARNING: {np.sum(n_failed)} total failed fits across all degrees")

    output = {
        # Full bootstrap distributions
        'V_bootstrap': V_bootstrap,           # (n_deg, n_bootstrap)
        'sigma_bootstrap': sigma_bootstrap,   # (n_deg, n_bootstrap)
        'chi2_bootstrap': chi2_bootstrap,     # (n_deg, n_bootstrap)
        'z_bootstrap': z_bootstrap,           # (n_deg, n_bootstrap)
        # Percentile-based errors (asymmetric)
        'sigma_boot_err_lo': sigma_boot_err_lo,  # median - 16th
        'sigma_boot_err_hi': sigma_boot_err_hi,  # 84th - median
        'sigma_p16': sigma_p16,
        'sigma_p50': sigma_p50,
        'sigma_p84': sigma_p84,
        'V_boot_err_lo': V_boot_err_lo,
        'V_boot_err_hi': V_boot_err_hi,
        'V_p16': V_p16,
        'V_p50': V_p50,
        'V_p84': V_p84,
        'z_bootstrap': z_bootstrap,
        'z_p16': z_p16,
        'z_p50': z_p50,
        'z_p84': z_p84,
        # Original fit values
        'V_original': V_original,
        'sigma_original': sigma_original,
        'chi2_original': chi2_original,
        'degrees': degrees,
        'sps_name': sps_name,
        'n_bootstrap': n_bootstrap,
        'n_failed': n_failed,
        'window': window,
        'seed': seed,
        'z_input': ppxf_inputs['z'],
    }

    if save:
        suffix = f'_{save_suffix}' if save_suffix else ''
        out_path = os.path.join(results_dir, f'ppxf_bootstrap_errors_{sps_name}{suffix}.npz')
        np.savez(out_path, **output)
        print(f"\nSaved: {out_path}")

    return output


# =============================================================================
# Regularization scan
# =============================================================================

def run_regularization_scan(ppxf_inputs, degree, regul_values=None):
    """
    Scan regularization parameter for a given degree.

    Following ppxf prescription: increase regul until chi2 rises by
    sqrt(2 * n_goodpixels) above unregularized chi2.

    Parameters
    ----------
    ppxf_inputs : dict
        Output of setup_ppxf_inputs().
    degree : int
        Polynomial degree.
    regul_values : array or None
        Values to scan. None = np.logspace(-2, 3, 50).

    Returns
    -------
    dict with keys:
        'regul_values' : array
        'chi2'         : array
        'V'            : array
        'sigma'        : array
        'weights'      : list of arrays (template weights per regul)
        'optimal_regul': float (regul where chi2 hits target)
        'chi2_target'  : float
        'chi2_unreg'   : float
    """
    galaxy = ppxf_inputs['galaxy']
    noise = ppxf_inputs['noise']
    velscale = ppxf_inputs['velscale']
    start = ppxf_inputs['start']
    goodpixels = ppxf_inputs['goodpixels']
    lam_gal_rest = ppxf_inputs['lam_gal_rest']
    sps = ppxf_inputs['sps']
    lam_temp = ppxf_inputs['lam_temp']

    if regul_values is None:
        regul_values = np.logspace(-2, 3, 50)

    # Unregularized fit first
    pp0 = ppxf(sps.templates, galaxy, noise, velscale, start,
               goodpixels=goodpixels, plot=False, moments=2, trig=False,
               degree=degree, lam=lam_gal_rest, lam_temp=lam_temp, mdegree=0)
    chi2_unreg = pp0.chi2
    chi2_target = chi2_unreg + np.sqrt(2 * len(goodpixels)) / len(goodpixels)

    print(f"Regularization scan at degree={degree}")
    print(f"  Unregularized chi2/DOF: {chi2_unreg:.2f}")
    print(f"  Target chi2/DOF: {chi2_target:.2f} (delta = sqrt(2*{len(goodpixels)})/{len(goodpixels)} = {np.sqrt(2*len(goodpixels))/len(goodpixels):.4f})")

    chi2_arr = np.zeros(len(regul_values))
    V_arr = np.zeros(len(regul_values))
    sigma_arr = np.zeros(len(regul_values))
    weights_list = []

    for i, regul in enumerate(tqdm(regul_values, desc="Regul scan")):
        pp = ppxf(sps.templates, galaxy, noise, velscale, start,
                  goodpixels=goodpixels, plot=False, moments=2, trig=False,
                  degree=degree, lam=lam_gal_rest, lam_temp=lam_temp, mdegree=0,
                  regul=regul)
        chi2_arr[i] = pp.chi2
        V_arr[i] = pp.sol[0]
        sigma_arr[i] = pp.sol[1]
        weights_list.append(pp.weights.copy())

    # Find optimal regul
    above_target = chi2_arr >= chi2_target
    if np.any(above_target):
        optimal_idx = np.argmax(above_target)
        optimal_regul = regul_values[optimal_idx]
    else:
        optimal_regul = regul_values[-1]
        optimal_idx = -1

    print(f"  Optimal regul: {optimal_regul:.2f}")
    print(f"  Sigma at optimal regul: {sigma_arr[optimal_idx]:.1f} km/s (vs {pp0.sol[1]:.1f} unregularized)")

    return {
        'regul_values': regul_values,
        'chi2': chi2_arr,
        'V': V_arr,
        'sigma': sigma_arr,
        'weights': weights_list,
        'optimal_regul': optimal_regul,
        'chi2_target': chi2_target,
        'chi2_unreg': chi2_unreg,
    }


# =============================================================================
# CLI
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Bootstrap ppxf errors for AGEL0206 velocity dispersion',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/bootstrap_ppxf.py --sps_name fsps --n_bootstrap 500
  python scripts/bootstrap_ppxf.py --sps_name all
  python scripts/bootstrap_ppxf.py --degrees 4,10,16,20 --n_bootstrap 50
        """)
    parser.add_argument('--sps_name', default='fsps',
                        help="Template library: 'fsps', 'emiles', 'xsl', or 'all' (default: fsps)")
    parser.add_argument('--n_bootstrap', type=int, default=500,
                        help='Number of bootstrap iterations (default: 500)')
    parser.add_argument('--window', type=int, default=75,
                        help='Rolling window for local residual scaling (default: 75)')
    parser.add_argument('--seed', type=int, default=42,
                        help='RNG seed (default: 42)')
    parser.add_argument('--ifu_file', default=DEFAULT_IFU,
                        help=f'Path to IFU cube (default: {DEFAULT_IFU})')
    parser.add_argument('--results_dir', default='results',
                        help='Results directory (default: results)')
    parser.add_argument('--degrees', type=str, default=None,
                        help='Comma-separated degrees, e.g. "4,10,16,20" (default: all)')
    args = parser.parse_args()

    degrees = None
    if args.degrees:
        degrees = np.array([int(d) for d in args.degrees.split(',')])

    sps_names = ['fsps', 'emiles', 'xsl'] if args.sps_name == 'all' else [args.sps_name]

    for sps_name in sps_names:
        run_bootstrap(
            ifu_file=args.ifu_file,
            sps_name=sps_name,
            results_dir=args.results_dir,
            degrees=degrees,
            n_bootstrap=args.n_bootstrap,
            window=args.window,
            seed=args.seed,
        )
