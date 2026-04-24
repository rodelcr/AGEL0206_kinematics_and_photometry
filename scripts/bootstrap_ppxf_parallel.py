"""
Parallel version of the wild-bootstrap ppxf inner loop.

The bootstrap iterations inside `run_bootstrap_single_degree` are independent
— each is a fresh ppxf call on a different realization of galaxy_boot. With
N_BOOTSTRAP=500 and ~0.05 s per call, a single degree takes ~25 s sequentially.

This module runs those iterations across a joblib process pool. BLAS threads
are pinned to 1 per worker via threadpoolctl to avoid oversubscription (the
default openblas will otherwise spawn 16 threads per worker on a 16-core Mac).

Measured speedup (n_jobs=8, 16-core Mac): ~5-6x for the bootstrap loop, which
translates to ~3-4x end-to-end because the non-bootstrap ppxf calls (best-fit
per degree) and the per-annulus setup also take some time but are not
parallelized here.

Public API
----------
run_bootstrap_single_degree_parallel(ppxf_inputs, degree, best_fit_spectrum,
                                      n_bootstrap, window, seed, n_jobs)

Drop-in replacement for `run_bootstrap_single_degree` from bootstrap_ppxf.py,
same return dict schema.
"""

import os
import numpy as np
from scipy.ndimage import uniform_filter1d

try:
    import joblib
    _HAVE_JOBLIB = True
except ImportError:
    _HAVE_JOBLIB = False

try:
    from threadpoolctl import threadpool_limits
    _HAVE_THREADPOOLCTL = True
except ImportError:
    _HAVE_THREADPOOLCTL = False

from ppxf.ppxf import ppxf


C_KMS = 299792.458


def _compute_local_residual_scaling(residuals, window=75):
    """Same formula as in bootstrap_ppxf.compute_local_residual_scaling."""
    global_std = np.std(residuals)
    if global_std == 0:
        return np.ones_like(residuals)
    local_mean = uniform_filter1d(residuals, size=window, mode='reflect')
    local_mean_sq = uniform_filter1d(residuals**2, size=window, mode='reflect')
    local_var = np.maximum(local_mean_sq - local_mean**2, 0)
    local_std = np.sqrt(local_var)
    return np.clip(local_std / global_std, 0.2, 5.0)


def _single_bootstrap(iteration, galaxy, noise, best_fit, scaled_residuals,
                      sps_templates, velscale, start, goodpixels,
                      lam_gal_rest, lam_temp, degree, seed, z0):
    """One bootstrap iteration. Runs inside a worker process with BLAS=1.

    Returns (V, sigma, chi2, z) or (nan, nan, nan, nan) on failure.
    """
    # Pin BLAS to 1 thread in this worker (prevents oversubscription)
    if _HAVE_THREADPOOLCTL:
        with threadpool_limits(limits=1):
            return _single_bootstrap_inner(
                iteration, galaxy, noise, best_fit, scaled_residuals,
                sps_templates, velscale, start, goodpixels,
                lam_gal_rest, lam_temp, degree, seed, z0,
            )
    else:
        return _single_bootstrap_inner(
            iteration, galaxy, noise, best_fit, scaled_residuals,
            sps_templates, velscale, start, goodpixels,
            lam_gal_rest, lam_temp, degree, seed, z0,
        )


def _single_bootstrap_inner(iteration, galaxy, noise, best_fit, scaled_residuals,
                            sps_templates, velscale, start, goodpixels,
                            lam_gal_rest, lam_temp, degree, seed, z0):
    rng = np.random.default_rng(seed + iteration)
    signs = rng.choice([-1, 1], size=len(galaxy))
    galaxy_boot = best_fit + scaled_residuals * signs
    try:
        pp = ppxf(sps_templates, galaxy_boot, noise, velscale, start,
                  goodpixels=goodpixels, plot=False, moments=2, trig=False,
                  degree=degree, lam=lam_gal_rest, lam_temp=lam_temp, mdegree=0,
                  quiet=True)
        V, sigma, chi2 = pp.sol[0], pp.sol[1], pp.chi2
        z_i = (1 + z0) * np.exp(V / C_KMS) - 1
        return V, sigma, chi2, z_i
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def run_bootstrap_single_degree_parallel(ppxf_inputs, degree, best_fit_spectrum,
                                          n_bootstrap=500, window=75, seed=None,
                                          n_jobs=8):
    """
    Parallel drop-in for run_bootstrap_single_degree.

    Parameters
    ----------
    ppxf_inputs : dict
        Output of setup_ppxf_inputs_from_spectrum().
    degree : int
        Additive polynomial degree.
    best_fit_spectrum : array shape (n_pixels,)
        Best-fit spectrum from the original (non-bootstrap) ppxf fit.
    n_bootstrap : int
        Number of bootstrap iterations.
    window : int
        Rolling window (pixels) for local residual scaling.
    seed : int
        RNG seed (each bootstrap iteration gets seed + i).
    n_jobs : int
        Number of parallel workers (default 8). Set to 1 to fall back to
        serial execution. Set to -1 to use all cores.

    Returns
    -------
    dict with keys: 'V_samples', 'sigma_samples', 'chi2_samples',
    'z_samples', 'n_failed'. Same schema as run_bootstrap_single_degree.
    """
    if not _HAVE_JOBLIB:
        raise RuntimeError('joblib not installed; cannot parallelize')

    galaxy = ppxf_inputs['galaxy']
    noise = ppxf_inputs['noise']
    velscale = ppxf_inputs['velscale']
    start = ppxf_inputs['start']
    goodpixels = ppxf_inputs['goodpixels']
    lam_gal_rest = ppxf_inputs['lam_gal_rest']
    sps = ppxf_inputs['sps']
    lam_temp = ppxf_inputs['lam_temp']
    z0 = ppxf_inputs.get('z', 0.67511)

    residuals = galaxy - best_fit_spectrum
    scale_factor = _compute_local_residual_scaling(residuals, window=window)
    scaled_residuals = residuals * scale_factor

    if seed is None:
        seed = 0

    sps_templates = sps.templates  # numpy array, picklable

    if n_jobs == 1:
        results = [
            _single_bootstrap_inner(
                i, galaxy, noise, best_fit_spectrum, scaled_residuals,
                sps_templates, velscale, start, goodpixels,
                lam_gal_rest, lam_temp, degree, seed, z0,
            )
            for i in range(n_bootstrap)
        ]
    else:
        results = joblib.Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
            joblib.delayed(_single_bootstrap)(
                i, galaxy, noise, best_fit_spectrum, scaled_residuals,
                sps_templates, velscale, start, goodpixels,
                lam_gal_rest, lam_temp, degree, seed, z0,
            )
            for i in range(n_bootstrap)
        )

    V_samples = np.array([r[0] for r in results])
    sigma_samples = np.array([r[1] for r in results])
    chi2_samples = np.array([r[2] for r in results])
    z_samples = np.array([r[3] for r in results])
    n_failed = int(np.sum(~np.isfinite(sigma_samples)))

    return {
        'V_samples': V_samples,
        'sigma_samples': sigma_samples,
        'chi2_samples': chi2_samples,
        'z_samples': z_samples,
        'n_failed': n_failed,
    }


if __name__ == '__main__':
    # Benchmark: compare serial vs parallel on a synthetic case
    import sys, os as _os
    sys.path.insert(0, _os.path.abspath(_os.path.dirname(__file__)))
    from bootstrap_ppxf import setup_ppxf_inputs, run_bootstrap_single_degree
    from time import perf_counter as clock

    IFU = _os.path.join(_os.path.dirname(__file__), '..',
                        'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits')
    print('Setting up ppxf inputs (FSPS, 190-spaxel) ...')
    inputs = setup_ppxf_inputs(IFU, sps_name='fsps', z=0.67564)

    from ppxf.ppxf import ppxf as _ppxf
    pp = _ppxf(inputs['sps'].templates, inputs['galaxy'], inputs['noise'],
               inputs['velscale'], inputs['start'],
               goodpixels=inputs['goodpixels'], plot=False, moments=2, trig=False,
               degree=20, lam=inputs['lam_gal_rest'],
               lam_temp=inputs['lam_temp'], mdegree=0, quiet=True)
    print(f'original fit: σ={pp.sol[1]:.1f} km/s')

    for label, fn, kw in [
        ('SERIAL (n_jobs=1)', run_bootstrap_single_degree_parallel, {'n_jobs': 1}),
        ('PARALLEL n_jobs=4', run_bootstrap_single_degree_parallel, {'n_jobs': 4}),
        ('PARALLEL n_jobs=8', run_bootstrap_single_degree_parallel, {'n_jobs': 8}),
        ('PARALLEL n_jobs=16', run_bootstrap_single_degree_parallel, {'n_jobs': 16}),
    ]:
        t0 = clock()
        r = fn(inputs, 20, pp.bestfit, n_bootstrap=50, window=75, seed=42, **kw)
        dt = clock() - t0
        sv = r['sigma_samples']; sv = sv[np.isfinite(sv)]
        print(f'{label:<22}  {dt:6.1f}s   σ_median={np.median(sv):.1f}   '
              f'std={np.std(sv):.2f}   n_failed={r["n_failed"]}')
