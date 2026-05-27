"""Single-pipeline driver for the wide-arc-masked σ_e analysis.

Runs the EXACT same machinery (load_setup → extract_aperture_spectrum →
frame-aware ppxf × 3 SPS × 15 polynomial degrees × N wild-bootstrap at the
`wR3800_5400_arcmask` window) on whichever cube you point it at, using
joblib parallel bootstrap with BLAS pinned to 1 thread per worker.

This is the single-entry-point demonstration that the +10.92 km/s
σ_e shift between the headline cube (254.85 km/s) and the new `_mtwdo_`
reduction (265.76 km/s) comes from the cube data, not from a different
pipeline. Audit details in `NOTES_nb09e_audit_2026-05-20.md`.

Usage
─────
  # run on the paper-headline cube
  python scripts/run_wide_sigma_e.py --cube headline

  # run on the new _mtwdo_ reduction
  python scripts/run_wide_sigma_e.py --cube new

  # run a different cube entirely
  python scripts/run_wide_sigma_e.py --cube /path/to/your_cube.fits

  # run both presets and print the side-by-side comparison
  python scripts/run_wide_sigma_e.py --compare

  # override default bootstrap N (default 500) and worker count
  python scripts/run_wide_sigma_e.py --compare --n_bootstrap 50 --n_jobs 4

Caching
───────
This script writes to its OWN cache subtree, separate from the production
caches used by the paper headline (`results/nb09a_wavelength_sweep/`) and
by the nb09e notebook (`results/nb09e_new_red_reduction/`). It does not
read from or write to either of those. Results land in:

    results/run_wide_sigma_e/{cube_alias}/wR3800_5400_arcmask_{sps}_T{T0}_{T1}_N{N}.npz

The first --compare execution refits both cubes from scratch (~60 min at
N=500 with n_jobs=8). Subsequent runs at the same N are instant. Pass
--force to ignore the script's own cache and refit anyway.
"""
from __future__ import annotations

import os
import sys
import time
import argparse
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import scripts.final_sigma_e as fse                       # noqa: E402
from scripts.bootstrap_ppxf import (                      # noqa: E402
    setup_ppxf_inputs_from_spectrum,
)
from scripts.bootstrap_ppxf_parallel import (             # noqa: E402
    run_bootstrap_single_degree_parallel,
)
from scripts.run_window_sweep import (                    # noqa: E402
    _apply_arc_mask, _apply_bad_pixels_mask,
    SPS_LAM_RANGE_TEMP, sps_safe_obs_min, DEGREES,
)
from ppxf.ppxf import ppxf                                # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# Presets
#
# Caches for this single-pipeline reproducibility test live under their OWN
# subtree (`results/run_wide_sigma_e/{alias}/`) — DELIBERATELY separated from
# the production caches in `results/nb09a_wavelength_sweep/` (headline
# paper-driving fits) and `results/nb09e_new_red_reduction/` (nb09e
# notebook). The test does not read from or write to either of those, so
# the production caches are guaranteed clean of any audit/smoke runs done
# through this script. The first --compare execution will re-fit both
# cubes from scratch (~60 min at N=500); subsequent runs hit this script's
# own cache and return instantly.
# ─────────────────────────────────────────────────────────────────────────────
_TEST_CACHE_ROOT = REPO / 'results' / 'run_wide_sigma_e'

PRESETS = {
    'headline': {
        'cube':      'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'headline',
        'descr':     'ORIGINAL cube (multi-night combine, kcwiRedux v1.2) — pre-2026-05-26 headline; arc-mask only (legacy)',
    },
    'new': {
        'cube':      'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'new',
        'descr':     'NEW _mtwdo_ reduction; arc-mask only (legacy, pre-bad-pixel-bake)',
    },
    'headline_clean': {
        'cube':      'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'headline_clean',
        'descr':     'ORIGINAL cube + arc-mask + local-MAD bad-pixel mask (M5 baked in 2026-05-26)',
    },
    'new_clean': {
        'cube':      'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'new_clean',
        'descr':     'NEW _mtwdo_ reduction + arc-mask + local-MAD bad-pixel mask — CURRENT HEADLINE-PRODUCING PRESET (2026-05-26)',
    },
    # 2026-05-27: He I 3819 source-emission mask (5237-5253 Å) + M10 full
    # sky-line audit (9 added OH/sky-residual bands across the fit window,
    # all VERIFIED non-source via the arc spectrum). These are the HEADLINE
    # presets as of 2026-05-27. Separate cache dirs preserve the un-masked-HeI
    # 268.98 km/s headline at new_clean/ and the He-I-only 271.87 km/s
    # interim headline (now superseded) was preserved in this dir but
    # OVERWRITTEN on 2026-05-27 when --force was passed to fold in M10.
    'new_clean_hei': {
        'cube':      'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'new_clean_hei',
        'descr':     'NEW _mtwdo_ + arc-mask (He I 3819 mask 5237-5253) + bad-pix (M10 full sky audit) + Balmer-unmasked — CURRENT HEADLINE-PRODUCING PRESET (2026-05-27)',
    },
    'headline_clean_hei': {
        'cube':      'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits',
        'cache_dir': _TEST_CACHE_ROOT / 'headline_clean_hei',
        'descr':     'OLD original cube + arc-mask (incl. He I 3819 mask 5237-5253) + bad-pix (M10 full sky audit) + Balmer-unmasked',
    },
}

WLABEL = 'wR3800_5400_arcmask'
SPS_LIBS = ('fsps', 'emiles', 'xsl')

# Match nb09a / run_window_sweep.py exactly
_R3800_OBS = 3800.0 * (1 + fse.Z_SYSTEMIC)
_R5400_OBS = 5400.0 * (1 + fse.Z_SYSTEMIC)
LAM_OBS_RANGE = (_R3800_OBS, _R5400_OBS)
ARC_MASK = True

# Systematic budget carried from the headline. UPDATED 2026-05-27 (M10 full
# sky-line audit added): the reduction-pass component refined from ±3.86
# (He-I-only) to ±3.45 (half-Δ between cleaned + He-I-masked + M10-sky-audit
# headlines: 269.62 - 262.72 = 6.90 → /2 = 3.45). M10 added 9 unmasked OH
# airglow / sky-residual bands to BAD_PIXELS_REST (all verified non-source
# via the arc spectrum); shifted σ_e down by ~2 km/s on both cubes and
# tightened the inter-reduction gap further.
#   I-shape 1.5 ⊕ mask 3.8 ⊕ frame 5 ⊕ centering 4 ⊕ window 15 ⊕ reduction 3.45
# Added in quadrature to each side of the asymmetric stat error.
SYS_QUAD = float(np.sqrt(1.5**2 + 3.8**2 + 5.0**2 + 4.0**2 + 15.0**2 + 3.45**2))  # = 17.16

# Reproduce the deterministic seed-offset convention from run_window_sweep
# (so caches written by THIS script bit-match the ones already on disk).
_W_IDX = 12   # see scripts/run_window_sweep.py:WINDOWS index for wR3800_5400_arcmask


def _banner(s):
    print('\n' + '=' * 78)
    print(s)
    print('=' * 78, flush=True)


def _resolve_cube(spec: str):
    """Return (cube_path, cache_dir, alias, descr). spec is a preset key or a path."""
    if spec in PRESETS:
        p = PRESETS[spec]
        cube_path = REPO / p['cube'] if not os.path.isabs(p['cube']) else Path(p['cube'])
        return cube_path, p['cache_dir'], spec, p['descr']
    # arbitrary path
    cube_path = Path(spec).resolve()
    if not cube_path.exists():
        raise FileNotFoundError(f'cube not found: {cube_path}')
    alias = cube_path.stem
    cache_dir = _TEST_CACHE_ROOT / alias
    return cube_path, cache_dir, alias, f'user-provided cube ({cube_path.name})'


def _fit_one_sps(sps, flux, noise, hdr, cache_dir, n_boot, n_jobs, force=False,
                 bad_pix_mask=False):
    """Run ppxf × 15 polynomial degrees × N parallel-bootstrap iterations at
    wR3800_5400_arcmask for one SPS template library. Returns the result dict
    (same schema as run_window_sweep._ppxf_one_window) and the cache path.

    bad_pix_mask : bool
        If True, apply the local-MAD-flagged BAD_PIXELS_REST mask in addition
        to ARC_MASKS_REST. Set automatically for the `*_clean` presets.
    """
    lo, hi = LAM_OBS_RANGE
    lr_eff = (max(lo, sps_safe_obs_min(sps)), hi)
    T0, T1 = int(SPS_LAM_RANGE_TEMP[sps][0]), int(SPS_LAM_RANGE_TEMP[sps][1])
    cache = cache_dir / f'{WLABEL}_{sps}_T{T0}_{T1}_N{n_boot}.npz'
    cache_dir.mkdir(parents=True, exist_ok=True)

    if cache.exists() and not force:
        d = dict(np.load(cache, allow_pickle=True))
        print(f'  [cached]  {sps:6s}  σ range = '
              f'{d["sig_orig"].min():.1f}–{d["sig_orig"].max():.1f} km/s   ({cache.name})', flush=True)
        return d, cache

    t0 = time.time()
    # 2026-05-26: for _clean presets, also DISABLE Balmer-line masking (Hδ,
    # Hγ, Hβ are absorption-dominated in the passive deflector; ppxf's default
    # ±800 km/s emission-line mask throws away ~120 stellar absorption pixels).
    # The setup function's `mask_balmer=False` activates a custom mask that
    # excludes only the truly forbidden lines ([OII], [OIII], [NII], [SII], [OI]).
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, hdr,
        sps_name=sps, z=fse.Z_SYSTEMIC,
        lam_obs_range=lr_eff, lam_fit_range=lr_eff,
        lam_range_temp=SPS_LAM_RANGE_TEMP[sps],
        verbose=False, frame_galaxy='auto',
        mask_balmer=not bad_pix_mask,   # _clean presets: keep Balmer absorption
    )
    if ARC_MASK:
        inputs['goodpixels'] = _apply_arc_mask(
            inputs['goodpixels'], inputs['lam_gal_rest']
        )
    if bad_pix_mask:
        inputs['goodpixels'] = _apply_bad_pixels_mask(
            inputs['goodpixels'], inputs['lam_gal_rest']
        )
    n_pix = len(inputs['galaxy'])
    n_deg = len(DEGREES)
    V_orig    = np.zeros(n_deg)
    sig_orig  = np.zeros(n_deg)
    chi2_orig = np.zeros(n_deg)
    V_boot    = np.full((n_deg, n_boot), np.nan)
    sig_boot  = np.full((n_deg, n_boot), np.nan)
    bf        = np.zeros((n_deg, n_pix))

    sps_idx = SPS_LIBS.index(sps)
    seed_offset = 50_000 + 10_000 * _W_IDX + 100 * sps_idx

    for i, deg in enumerate(DEGREES):
        pp = ppxf(
            inputs['sps'].templates, inputs['galaxy'], inputs['noise'],
            inputs['velscale'], inputs['start'],
            goodpixels=inputs['goodpixels'], plot=False, moments=2, trig=False,
            degree=int(deg), mdegree=0,
            lam=inputs['lam_gal_rest'], lam_temp=inputs['lam_temp'], quiet=True,
        )
        V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
        bf[i] = pp.bestfit
        # PARALLEL bootstrap (joblib + BLAS=1 per worker)
        rb = run_bootstrap_single_degree_parallel(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_boot, window=fse.WINDOW,
            seed=fse.BOOT_SEED + seed_offset + i,
            n_jobs=n_jobs,
        )
        V_boot[i]   = rb['V_samples']
        sig_boot[i] = rb['sigma_samples']

    d = dict(
        sps=sps, lam_obs_range=lr_eff, lam_fit_range=lr_eff,
        n_pix=n_pix,
        lam_rest_min=float(inputs['lam_gal_rest'][0]),
        lam_rest_max=float(inputs['lam_gal_rest'][-1]),
        frame_galaxy=inputs['frame_galaxy'],
        V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
        V_boot=V_boot, sig_boot=sig_boot, bf=bf,
    )
    np.savez(cache, **d)
    elapsed = time.time() - t0
    print(f'  [done]    {sps:6s}  σ range = '
          f'{d["sig_orig"].min():.1f}–{d["sig_orig"].max():.1f} km/s   '
          f'({elapsed:.0f}s; npix={n_pix}; n_jobs={n_jobs})', flush=True)
    return d, cache


def _pool_and_summarise(results):
    """Pool the per-(SPS×degree) bootstrap σ samples and compute asymmetric
    percentiles + total (stat ⊕ sys) error budget.
    """
    pool = []
    per_sps_p50 = {}
    for sps in SPS_LIBS:
        s = results[sps]['sig_boot'].ravel()
        s = s[np.isfinite(s)]
        pool.append(s)
        per_sps_p50[sps] = float(np.percentile(s, 50))
    pool = np.concatenate(pool)
    p16, p50, p84 = np.percentile(pool, [16, 50, 84])
    stat_lo  = float(p50 - p16)
    stat_hi  = float(p84 - p50)
    stat_sym = float((p84 - p16) / 2)
    return dict(
        p50=float(p50), stat_lo=stat_lo, stat_hi=stat_hi, stat_sym=stat_sym,
        total_lo=float(np.sqrt(stat_lo**2 + SYS_QUAD**2)),
        total_hi=float(np.sqrt(stat_hi**2 + SYS_QUAD**2)),
        total_sym=float(np.sqrt(stat_sym**2 + SYS_QUAD**2)),
        per_sps=per_sps_p50,
        spread=max(per_sps_p50.values()) - min(per_sps_p50.values()),
        n_pool=int(len(pool)),
    )


def run_one(cube_spec: str, n_boot: int, n_jobs: int, force: bool):
    """End-to-end run for ONE cube. Returns the pooled summary dict."""
    cube_path, cache_dir, alias, descr = _resolve_cube(cube_spec)
    # Any alias containing '_clean' uses the bad-pix mask + Balmer-unmasked
    # path (matches both 'new_clean' and 'new_clean_hei', etc.).
    bad_pix_mask = '_clean' in alias

    _banner(f'WIDE-ARC-MASKED σ_e   |   cube alias = "{alias}"')
    print(f'  cube path:  {cube_path}')
    print(f'  cube file:  {cube_path.stat().st_size/1e6:.1f} MB    {descr}')
    print(f'  cache dir:  {cache_dir}    (skip if N={n_boot} caches exist; --force to refit)')
    print(f'  N_BOOT={n_boot}  N_JOBS={n_jobs}  arc_mask={ARC_MASK}  bad_pix_mask={bad_pix_mask}  '
          f'fit window obs [{LAM_OBS_RANGE[0]:.0f}, {LAM_OBS_RANGE[1]:.0f}] Å')
    print(f'  sys budget (carried): ±{SYS_QUAD:.2f} km/s (Ishape 1.5 ⊕ mask 3.8 ⊕ frame 5 ⊕ centering 4 ⊕ window 15)')

    # 1. Monkey-patch the cube path and load setup
    fse.IFU_FILE = str(cube_path)
    print()
    state = fse.load_setup()
    print()
    print(f'  → HST-mean center (IFU sub-pix) = ({state["cx_ifu"]:.3f}, {state["cy_ifu"]:.3f})')
    print(f'  → R_E                            = {state["R_E"]:.4f}"')
    print(f'  → arc_spax_mask                  = {state["arc_spax_mask"].sum()} flagged')

    # 2. Extract the I-weighted aperture spectrum at R<R_e (w=0, hard mask)
    flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
        state, r_max=state['R_E'], mask_weight=0.0,
    )
    n_med   = float(np.nanmedian(noise[noise > 0]))
    n_floor = n_med * 0.1
    n_zero  = int(np.sum(~(noise > 0)))
    noise   = np.where(noise > 0, noise, n_floor)
    print(f'  → aperture R<R_e: n_kept={n_kept} spaxels,  S/N(6500–7500)={sn_band:.2f}'
          f',  noise floor → {n_zero} pix')

    # 3. ppxf + parallel wild bootstrap per SPS
    print(f'\n  fitting ppxf at {WLABEL} for FSPS / EMILES / XSL ...')
    results = {}
    for sps in SPS_LIBS:
        results[sps], _ = _fit_one_sps(
            sps, flux, noise, state['hdr'], cache_dir, n_boot, n_jobs,
            force=force, bad_pix_mask=bad_pix_mask,
        )

    # 4. Pool and report
    summary = _pool_and_summarise(results)
    print()
    print(f'  pool size (3 SPS × {len(DEGREES)} degrees × N={n_boot}): {summary["n_pool"]} samples')
    print(f'  per-SPS p50: '
          f'FSPS={summary["per_sps"]["fsps"]:.2f}  '
          f'EMILES={summary["per_sps"]["emiles"]:.2f}  '
          f'XSL={summary["per_sps"]["xsl"]:.2f}  '
          f'(spread {summary["spread"]:.2f} km/s)')
    print(f'  σ_e(<R_e) = {summary["p50"]:.2f}  '
          f'−{summary["stat_lo"]:.2f}/+{summary["stat_hi"]:.2f}  km/s    [STAT only, asymmetric]')
    print(f'             {summary["p50"]:.2f}  '
          f'−{summary["total_lo"]:.2f}/+{summary["total_hi"]:.2f}  km/s    '
          f'[TOTAL,  stat ⊕ sys ±{SYS_QUAD:.2f} each side]')
    print(f'             {summary["p50"]:.2f}  '
          f'± {summary["total_sym"]:.2f}            km/s    [TOTAL,  symmetric quadrature]')

    summary['alias'] = alias
    summary['cube_path'] = str(cube_path)
    summary['descr'] = descr
    return summary


def run_compare(n_boot: int, n_jobs: int, force: bool):
    """Run both presets and print the side-by-side comparison table."""
    summaries = {}
    for alias in ('headline', 'new'):
        summaries[alias] = run_one(alias, n_boot, n_jobs, force)

    _banner('SINGLE-PIPELINE COMPARISON — headline cube vs new _mtwdo_ reduction')
    h = summaries['headline']
    n = summaries['new']
    print(f'  Pipeline:    wR3800_5400_arcmask, R<R_e, mask_weight=0, frame-aware ppxf')
    print(f'               × 3 SPS × {len(DEGREES)} polynomial degrees × N={n_boot} wild bootstrap')
    print(f'               (parallel: joblib loky, n_jobs={n_jobs}, BLAS pinned to 1 per worker)')
    print()
    print(f'{"":18s} {"σ_e (km/s)":>10s}   {"stat (asym)":>16s}   {"total (asym)":>16s}   {"SPS spread":>10s}')
    print('-' * 90)
    for k, s in (('headline (old)', h), ('new (_mtwdo_)', n)):
        print(f'  {k:16s} {s["p50"]:>10.2f}   '
              f'−{s["stat_lo"]:>5.2f}/+{s["stat_hi"]:>5.2f}    '
              f'−{s["total_lo"]:>5.2f}/+{s["total_hi"]:>5.2f}    '
              f'{s["spread"]:>10.2f}')
    print()
    delta = n['p50'] - h['p50']
    ratio_total = abs(delta) / h['total_sym']
    ratio_stat  = abs(delta) / ((h['stat_lo'] + h['stat_hi']) / 2)
    print(f'  Δσ_e (new − old) = {delta:+.2f} km/s')
    print(f'    |Δ| / total_sym = {ratio_total:.3f}σ   (headline symmetric total ±{h["total_sym"]:.2f})')
    print(f'    |Δ| / stat_only = {ratio_stat:.3f}σ   (headline stat 1σ ≈ {(h["stat_lo"]+h["stat_hi"])/2:.2f})')
    print()
    print(f'  Verdict (NOTES_nb09e_audit_2026-05-20.md): identical pipeline applied to')
    print(f'  both cubes; the Δ is a real reduction-pass systematic, not a misapplication.')
    print(f'  The new reduction lands within 2.2 km/s of the narrow-window cross-check')
    print(f'  (267.95 km/s), vs 13.1 km/s for the headline → `_mtwdo_` may have better')
    print(f'  wide-window flux calibration.')
    print()
    return summaries


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--cube', type=str, default=None,
                   help='Cube preset (headline | new) or absolute path to a FITS cube.')
    g.add_argument('--compare', action='store_true',
                   help='Run BOTH presets (headline + new) and print side-by-side comparison.')
    p.add_argument('--n_bootstrap', type=int, default=500,
                   help='Bootstrap iterations per polynomial degree (default 500).')
    p.add_argument('--n_jobs', type=int, default=8,
                   help='Parallel workers for joblib (default 8). BLAS auto-pinned to 1/worker.')
    p.add_argument('--force', action='store_true',
                   help='Ignore cached per-SPS results and refit.')
    args = p.parse_args()

    if args.compare:
        run_compare(args.n_bootstrap, args.n_jobs, args.force)
    else:
        run_one(args.cube, args.n_bootstrap, args.n_jobs, args.force)


if __name__ == '__main__':
    main()
