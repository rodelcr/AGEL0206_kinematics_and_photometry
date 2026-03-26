"""
Spectroscopic redshift verification from emission and absorption lines.

Fits Gaussian profiles to known spectral lines in the AGEL0206 KCWI IFU
data to independently verify the deflector and source redshifts.

Inspired by Keerthi Vasan's Keck-AGELDR2 redshift measurement pipeline
(Code_ESI/ and Code_Nires/ in 20250910-keerthi-Keck-AGELDR2-main/).

Approach:
  1. Load the integrated deflector spectrum (spatially averaged over deflector region)
  2. Identify expected positions of key absorption/emission lines at trial redshifts
  3. Fit Gaussian profiles (single or multi-line) to measure line centroids
  4. Derive redshift from each line independently and as a weighted average
  5. Report per-line redshifts with uncertainties

Usage:
    from scripts.redshift_verify import verify_redshift, plot_line_fits

Requires: numpy, scipy, astropy, matplotlib
"""

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt


# =============================================================================
# Line dictionaries
# =============================================================================

# Absorption lines (deflector galaxy — early-type stellar features)
ABSORPTION_LINES = {
    'Ca K':    {'lambda': 3933.66, 'label': 'Ca K'},
    'Ca H':    {'lambda': 3968.47, 'label': 'Ca H'},
    'H-delta': {'lambda': 4101.74, 'label': r'H$\delta$'},
    'G-band':  {'lambda': 4304.40, 'label': 'G-band'},
    'H-gamma': {'lambda': 4340.47, 'label': r'H$\gamma$'},
    'Fe4383':  {'lambda': 4383.55, 'label': 'Fe4383'},
    'H-beta':  {'lambda': 4861.33, 'label': r'H$\beta$'},
    'Mg I':    {'lambda': 5175.27, 'label': 'Mg I'},
    'Fe I':    {'lambda': 5270.40, 'label': 'Fe I'},
    'Na D':    {'lambda': 5893.00, 'label': 'Na D'},
}

# Emission lines (source galaxy — star-forming features)
EMISSION_LINES = {
    '[O II] 3726': {'lambda': 3726.03, 'label': '[O II]'},
    '[O II] 3729': {'lambda': 3728.82, 'label': '[O II]'},
    'H-beta':      {'lambda': 4861.33, 'label': r'H$\beta$'},
    '[O III] 4959': {'lambda': 4958.91, 'label': '[O III]'},
    '[O III] 5007': {'lambda': 5006.84, 'label': '[O III]'},
    'H-alpha':     {'lambda': 6562.80, 'label': r'H$\alpha$'},
}


# =============================================================================
# Gaussian fitting functions
# =============================================================================

def gaussian_absorption(x, amp, center, sigma, cont):
    """Single Gaussian absorption line + constant continuum."""
    return cont - amp * np.exp(-0.5 * ((x - center) / sigma)**2)


def gaussian_emission(x, amp, center, sigma, cont):
    """Single Gaussian emission line + constant continuum."""
    return cont + amp * np.exp(-0.5 * ((x - center) / sigma)**2)


def double_gaussian_absorption(x, amp1, amp2, center1, center2, sigma, cont):
    """Two Gaussian absorption lines sharing the same width + continuum.
    Used for Ca H+K doublet."""
    g1 = amp1 * np.exp(-0.5 * ((x - center1) / sigma)**2)
    g2 = amp2 * np.exp(-0.5 * ((x - center2) / sigma)**2)
    return cont - g1 - g2


def oii_doublet(x, amp1, ratio, z, sigma, cont):
    """[O II] 3726/3729 doublet emission with shared redshift and width."""
    z1 = 1.0 + z
    g1 = amp1 * np.exp(-0.5 * ((x - z1 * 3726.03) / sigma)**2)
    g2 = (amp1 / ratio) * np.exp(-0.5 * ((x - z1 * 3728.82) / sigma)**2)
    return cont + g1 + g2


# =============================================================================
# Line fitting
# =============================================================================

def fit_single_line(wavelength, flux, rest_lambda, z_guess, window_angstrom=30.0,
                    emission=False, noise=None):
    """
    Fit a single Gaussian to one spectral line.

    Parameters
    ----------
    wavelength : array
        Observed wavelength array (Angstroms).
    flux : array
        Flux array (same length as wavelength).
    rest_lambda : float
        Rest-frame wavelength of the line (Angstroms).
    z_guess : float
        Initial redshift guess.
    window_angstrom : float
        Half-width of fitting window around expected line position (Angstroms).
    emission : bool
        If True, fit emission line (positive Gaussian). If False, absorption.
    noise : array or None
        Noise spectrum for weighting. If None, uniform weights.

    Returns
    -------
    dict with keys:
        'z_fit'      : fitted redshift
        'z_err'      : redshift uncertainty (from covariance)
        'center'     : fitted line center (Angstroms, observed)
        'sigma'      : fitted Gaussian sigma (Angstroms)
        'amp'        : fitted amplitude
        'cont'       : fitted continuum level
        'popt'       : full parameter vector
        'pcov'       : covariance matrix
        'success'    : bool, whether fit converged
        'wl_window'  : wavelength array in fitting window
        'flux_window': flux array in fitting window
        'model'      : best-fit model in window
    """
    obs_lambda = rest_lambda * (1 + z_guess)

    # Extract fitting window
    mask = (wavelength > obs_lambda - window_angstrom) & \
           (wavelength < obs_lambda + window_angstrom)
    if np.sum(mask) < 10:
        return {'success': False, 'z_fit': np.nan, 'z_err': np.nan,
                'message': 'Too few pixels in window'}

    wl = wavelength[mask]
    fl = flux[mask]
    sigma_weights = noise[mask] if noise is not None else None

    # Initial guesses
    cont_guess = np.median(fl)
    amp_guess = abs(np.max(fl) - np.min(fl)) * 0.5
    sigma_guess = 3.0  # Angstroms (~150 km/s at 6000 Å)

    func = gaussian_emission if emission else gaussian_absorption
    p0 = [amp_guess, obs_lambda, sigma_guess, cont_guess]
    bounds_lo = [0, obs_lambda - 15, 0.5, -np.inf]
    bounds_hi = [np.inf, obs_lambda + 15, 20.0, np.inf]

    try:
        popt, pcov = curve_fit(func, wl, fl, p0=p0,
                                bounds=(bounds_lo, bounds_hi),
                                sigma=sigma_weights, absolute_sigma=True,
                                maxfev=5000)
        perr = np.sqrt(np.diag(pcov))

        center_fit = popt[1]
        center_err = perr[1]
        z_fit = center_fit / rest_lambda - 1
        z_err = center_err / rest_lambda

        model = func(wl, *popt)

        return {
            'z_fit': z_fit,
            'z_err': z_err,
            'center': center_fit,
            'center_err': center_err,
            'sigma': popt[2],
            'amp': popt[0],
            'cont': popt[3],
            'popt': popt,
            'pcov': pcov,
            'success': True,
            'wl_window': wl,
            'flux_window': fl,
            'model': model,
        }
    except (RuntimeError, ValueError) as e:
        return {'success': False, 'z_fit': np.nan, 'z_err': np.nan,
                'message': str(e)}


def fit_cahk_doublet(wavelength, flux, z_guess, window_angstrom=60.0, noise=None):
    """
    Fit the Ca H+K doublet simultaneously with shared velocity width.

    Returns fitted redshift from the Ca K line center.
    """
    obs_cak = 3933.66 * (1 + z_guess)
    obs_cah = 3968.47 * (1 + z_guess)

    center = (obs_cak + obs_cah) / 2
    mask = (wavelength > center - window_angstrom) & \
           (wavelength < center + window_angstrom)
    if np.sum(mask) < 15:
        return {'success': False, 'z_fit': np.nan, 'z_err': np.nan}

    wl = wavelength[mask]
    fl = flux[mask]
    sigma_weights = noise[mask] if noise is not None else None

    cont_guess = np.median(fl)
    amp_guess = abs(np.max(fl) - np.min(fl)) * 0.3
    sigma_guess = 4.0

    p0 = [amp_guess, amp_guess, obs_cak, obs_cah, sigma_guess, cont_guess]
    bounds_lo = [0, 0, obs_cak - 15, obs_cah - 15, 0.5, -np.inf]
    bounds_hi = [np.inf, np.inf, obs_cak + 15, obs_cah + 15, 20, np.inf]

    try:
        popt, pcov = curve_fit(double_gaussian_absorption, wl, fl, p0=p0,
                                bounds=(bounds_lo, bounds_hi),
                                sigma=sigma_weights, absolute_sigma=True,
                                maxfev=5000)
        perr = np.sqrt(np.diag(pcov))

        # Redshift from Ca K center
        z_cak = popt[2] / 3933.66 - 1
        z_cak_err = perr[2] / 3933.66
        # Redshift from Ca H center
        z_cah = popt[3] / 3968.47 - 1
        z_cah_err = perr[3] / 3968.47

        model = double_gaussian_absorption(wl, *popt)

        return {
            'z_cak': z_cak, 'z_cak_err': z_cak_err,
            'z_cah': z_cah, 'z_cah_err': z_cah_err,
            'z_fit': (z_cak + z_cah) / 2,
            'z_err': np.sqrt(z_cak_err**2 + z_cah_err**2) / 2,
            'sigma': popt[4],
            'popt': popt, 'pcov': pcov,
            'success': True,
            'wl_window': wl, 'flux_window': fl, 'model': model,
        }
    except (RuntimeError, ValueError) as e:
        return {'success': False, 'z_fit': np.nan, 'z_err': np.nan,
                'message': str(e)}


# =============================================================================
# Main verification function
# =============================================================================

def verify_redshift(wavelength, flux, z_guess, noise=None,
                    line_dict=None, emission=False,
                    window_angstrom=30.0, wl_range=None):
    """
    Fit all lines in a dictionary and report per-line and weighted-average redshifts.

    Parameters
    ----------
    wavelength : array
        Observed wavelength array (Angstroms).
    flux : array
        Flux array.
    z_guess : float
        Initial redshift guess.
    noise : array or None
        Noise spectrum.
    line_dict : dict or None
        Line dictionary (default: ABSORPTION_LINES if not emission, EMISSION_LINES if emission).
    emission : bool
        Whether to fit emission or absorption lines.
    window_angstrom : float
        Half-width of fitting window per line.
    wl_range : tuple or None
        (min, max) observed wavelength range to restrict which lines are attempted.

    Returns
    -------
    dict with keys:
        'per_line'      : dict of {line_name: fit_result}
        'z_weighted'    : inverse-variance weighted average redshift
        'z_weighted_err': uncertainty on weighted average
        'z_median'      : median of per-line redshifts
        'n_lines_fit'   : number of lines successfully fit
    """
    if line_dict is None:
        line_dict = EMISSION_LINES if emission else ABSORPTION_LINES

    if wl_range is None:
        wl_range = (wavelength.min(), wavelength.max())

    results = {}
    for name, info in line_dict.items():
        obs_wl = info['lambda'] * (1 + z_guess)
        if obs_wl < wl_range[0] or obs_wl > wl_range[1]:
            continue

        result = fit_single_line(wavelength, flux, info['lambda'], z_guess,
                                  window_angstrom=window_angstrom,
                                  emission=emission, noise=noise)
        result['rest_lambda'] = info['lambda']
        result['label'] = info['label']
        results[name] = result

    # Also try Ca H+K doublet if both are in range
    if not emission:
        obs_cak = 3933.66 * (1 + z_guess)
        obs_cah = 3968.47 * (1 + z_guess)
        if obs_cak > wl_range[0] and obs_cah < wl_range[1]:
            cahk = fit_cahk_doublet(wavelength, flux, z_guess,
                                     window_angstrom=60.0, noise=noise)
            cahk['rest_lambda'] = 3933.66
            cahk['label'] = 'Ca H+K doublet'
            results['Ca H+K doublet'] = cahk

    # Compute weighted average
    z_vals = []
    z_errs = []
    for name, r in results.items():
        if r['success'] and np.isfinite(r['z_fit']) and np.isfinite(r['z_err']) and r['z_err'] > 0:
            z_vals.append(r['z_fit'])
            z_errs.append(r['z_err'])

    z_vals = np.array(z_vals)
    z_errs = np.array(z_errs)

    if len(z_vals) > 0:
        weights = 1.0 / z_errs**2
        z_weighted = np.sum(z_vals * weights) / np.sum(weights)
        z_weighted_err = 1.0 / np.sqrt(np.sum(weights))
        z_median = np.median(z_vals)
    else:
        z_weighted = np.nan
        z_weighted_err = np.nan
        z_median = np.nan

    return {
        'per_line': results,
        'z_weighted': z_weighted,
        'z_weighted_err': z_weighted_err,
        'z_median': z_median,
        'n_lines_fit': len(z_vals),
    }


# =============================================================================
# Plotting
# =============================================================================

def plot_line_fits(verify_result, title='Redshift verification'):
    """
    Plot each fitted line with its Gaussian model, and a summary panel.

    Parameters
    ----------
    verify_result : dict
        Output of verify_redshift().
    title : str
        Figure title.
    """
    per_line = verify_result['per_line']
    successful = {k: v for k, v in per_line.items() if v.get('success', False)}
    n = len(successful)
    if n == 0:
        print("No successful line fits to plot.")
        return

    ncols = min(4, n + 1)
    nrows = int(np.ceil((n + 1) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = axes.ravel()

    # Individual line panels
    for i, (name, r) in enumerate(successful.items()):
        ax = axes[i]
        ax.step(r['wl_window'], r['flux_window'], 'k', lw=1, label='Data')
        ax.plot(r['wl_window'], r['model'], 'r', lw=1.5, label='Fit')
        ax.axvline(r.get('center', r['rest_lambda'] * (1 + r['z_fit'])),
                   color='blue', ls='--', alpha=0.5)
        ax.set_title(f"{r.get('label', name)}\nz = {r['z_fit']:.5f} ± {r['z_err']:.5f}",
                     fontsize=10)
        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Flux')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Summary panel
    ax_sum = axes[n]
    z_vals = [r['z_fit'] for r in successful.values()]
    z_errs = [r['z_err'] for r in successful.values()]
    labels = [r.get('label', k) for k, r in successful.items()]

    y_pos = np.arange(len(z_vals))
    ax_sum.errorbar(z_vals, y_pos, xerr=z_errs, fmt='o', color='k', capsize=4)
    ax_sum.set_yticks(y_pos)
    ax_sum.set_yticklabels(labels, fontsize=9)
    ax_sum.axvline(verify_result['z_weighted'], color='r', ls='-', lw=2,
                   label=f"Weighted: {verify_result['z_weighted']:.5f}")
    ax_sum.axvline(verify_result['z_median'], color='b', ls='--', lw=1.5,
                   label=f"Median: {verify_result['z_median']:.5f}")
    ax_sum.set_xlabel('Redshift')
    ax_sum.set_title('Per-line redshifts')
    ax_sum.legend(fontsize=9)
    ax_sum.grid(alpha=0.3)

    # Hide unused axes
    for ax in axes[n + 1:]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    return fig


def plot_spectrum_with_lines(wavelength, flux, z, line_dict=None,
                              emission=False, title='Spectrum with line IDs'):
    """
    Plot full spectrum with vertical lines at expected positions of known lines.

    Parameters
    ----------
    wavelength : array
        Observed wavelength (Angstroms).
    flux : array
        Flux array.
    z : float
        Redshift to use for line positions.
    line_dict : dict or None
        Line dictionary (default: ABSORPTION_LINES or EMISSION_LINES).
    emission : bool
        Whether lines are emission (for label color).
    title : str
        Plot title.
    """
    if line_dict is None:
        line_dict = EMISSION_LINES if emission else ABSORPTION_LINES

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.step(wavelength, flux, 'k', lw=0.8)

    color = 'blue' if emission else 'red'
    for name, info in line_dict.items():
        obs_wl = info['lambda'] * (1 + z)
        if obs_wl < wavelength.min() or obs_wl > wavelength.max():
            continue
        ax.axvline(obs_wl, color=color, ls='--', alpha=0.5, lw=0.8)
        ax.text(obs_wl, ax.get_ylim()[1] * 0.92, info['label'],
                fontsize=8, ha='center', color=color, rotation=45)

    ax.set_xlabel('Observed Wavelength (Å)')
    ax.set_ylabel('Flux')
    ax.set_title(f"{title} (z = {z:.5f})")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    return fig
