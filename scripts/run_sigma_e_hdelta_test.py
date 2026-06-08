"""Hδ targeted-masking decision (Task 6, 2026-06-08).

Open TODO (`bootstrap_ppxf.py:_determine_goodpixels_no_balmer`): the headline σ_e
fit keeps Hδ (4101.74 Å), Hγ, Hβ IN the fit as stellar absorption. Hδ is a
higher-order Balmer line whose absorption shape is subtler and more sensitive to
template Balmer-decrement mismatch, and was flagged mid-implementation as a
possible problem region. Decide between:
  (a) keep all Balmer unmasked          [current headline]
  (b) mask Hδ at the ppxf default ±800 km/s
  (c) mask Hδ at a narrow ±300 km/s

Two diagnostics, both on the headline pipeline (NEW `_mtwdo_` cube,
wR3800_5400_arcmask, He I + M10 masks, Balmer-unmasked), holding everything
fixed but the Hδ goodpixels band:
  1. Local-MAD outlier check — is the standardized fit residual at Hδ an
     outlier vs the rest of the fit window? (cheap; no re-fit)
  2. Δσ_e — re-run ppxf × 3 SPS × 15 deg × N bootstrap under (a)/(b)/(c) and
     compare the pooled σ_e.

Usage:  conda activate ISMgas
        python scripts/run_sigma_e_hdelta_test.py --n_bootstrap 500 --n_jobs 8
Summary: results/sigma_e_hdelta_test_N{N}.npz
"""
from __future__ import annotations
import os, sys, time, argparse
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); os.chdir(REPO)

import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum
from scripts.bootstrap_ppxf_parallel import run_bootstrap_single_degree_parallel
from scripts.run_window_sweep import (_apply_arc_mask, _apply_bad_pixels_mask,
                                      SPS_LAM_RANGE_TEMP, sps_safe_obs_min, DEGREES)
from ppxf.ppxf import ppxf

HEADLINE_CUBE = REPO / 'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits'
CACHE_ROOT = REPO / 'results' / 'run_wide_sigma_e'
HDELTA_AIR = 4101.74   # Å, rest (air, NIST)
C_KMS = 299792.458


def _apply_hdelta_mask(goodpixels, lam_rest, width_kms):
    """Drop a ±width_kms band around Hδ (rest frame) from goodpixels."""
    bad = np.ones(len(lam_rest), dtype=bool); bad[goodpixels] = False
    half = HDELTA_AIR * width_kms / C_KMS
    bad |= (lam_rest >= HDELTA_AIR - half) & (lam_rest <= HDELTA_AIR + half)
    return np.where(~bad)[0]


def _fit_variant(sps, flux, noise, hdr, hdelta_width, n_boot, n_jobs):
    """Replicates run_wide_sigma_e._fit_one_sps (headline config: bad_pix_mask,
    Balmer-unmasked) but optionally adds an Hδ mask of ±hdelta_width km/s.
    hdelta_width=None → variant (a), the headline. Returns the result dict."""
    lo, hi = rws.LAM_OBS_RANGE
    lr_eff = (max(lo, sps_safe_obs_min(sps)), hi)
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, hdr, sps_name=sps, z=fse.Z_SYSTEMIC,
        lam_obs_range=lr_eff, lam_fit_range=lr_eff,
        lam_range_temp=SPS_LAM_RANGE_TEMP[sps], verbose=False,
        frame_galaxy='auto', mask_balmer=False)
    inputs['goodpixels'] = _apply_arc_mask(inputs['goodpixels'], inputs['lam_gal_rest'])
    inputs['goodpixels'] = _apply_bad_pixels_mask(inputs['goodpixels'], inputs['lam_gal_rest'])
    if hdelta_width is not None:
        inputs['goodpixels'] = _apply_hdelta_mask(
            inputs['goodpixels'], inputs['lam_gal_rest'], hdelta_width)
    n_deg = len(DEGREES)
    V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
    sig_boot = np.full((n_deg, n_boot), np.nan)
    sps_idx = rws.SPS_LIBS.index(sps)
    seed_offset = 50_000 + 10_000 * rws._W_IDX + 100 * sps_idx
    resid_hd = None
    for i, deg in enumerate(DEGREES):
        pp = ppxf(inputs['sps'].templates, inputs['galaxy'], inputs['noise'],
                  inputs['velscale'], inputs['start'], goodpixels=inputs['goodpixels'],
                  plot=False, moments=2, trig=False, degree=int(deg), mdegree=0,
                  lam=inputs['lam_gal_rest'], lam_temp=inputs['lam_temp'], quiet=True)
        V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
        rb = run_bootstrap_single_degree_parallel(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_boot, window=fse.WINDOW,
            seed=fse.BOOT_SEED + seed_offset + i, n_jobs=n_jobs)
        sig_boot[i] = rb['sigma_samples']
        if i == n_deg // 2 and hdelta_width is None:
            # local-MAD residual diagnostic at the median degree (headline only)
            lamr = inputs['lam_gal_rest']; gp = inputs['goodpixels']
            resid = (inputs['galaxy'] - pp.bestfit)
            std = resid / inputs['noise']
            mad = 1.4826 * np.median(np.abs(std[gp] - np.median(std[gp])))
            half = HDELTA_AIR * 300 / C_KMS
            hd = (lamr >= HDELTA_AIR - half) & (lamr <= HDELTA_AIR + half)
            in_gp = np.zeros(len(lamr), bool); in_gp[gp] = True
            hd_gp = hd & in_gp
            resid_hd = dict(
                z_max_hd=float(np.max(np.abs(std[hd_gp]))) if hd_gp.any() else np.nan,
                z_med_hd=float(np.median(np.abs(std[hd_gp]))) if hd_gp.any() else np.nan,
                mad_global=float(mad), n_hd_pix=int(hd_gp.sum()),
                z_p99_global=float(np.percentile(np.abs(std[gp]), 99)))
    return dict(sig_boot=sig_boot, sig_orig=sig_orig, resid_hd=resid_hd)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n_bootstrap', type=int, default=500)
    p.add_argument('--n_jobs', type=int, default=8)
    args = p.parse_args(); N = args.n_bootstrap

    rws._banner(f'Hδ TARGETED-MASKING DECISION  (N={N})')
    fse.IFU_FILE = str(HEADLINE_CUBE)
    state = fse.load_setup()
    flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
        state, r_max=state['R_E'], mask_weight=0.0)
    n_med = float(np.nanmedian(noise[noise > 0]))
    noise = np.where(noise > 0, noise, n_med * 0.1)
    print(f'\n  aperture R<R_e: n_kept={n_kept}  S/N={sn_band:.2f}\n')

    variants = [('a_unmasked', None), ('b_hd800', 800.0), ('c_hd300', 300.0)]
    summaries = {}
    diag = None
    for key, width in variants:
        print(f'  ── variant {key}  (Hδ mask width = {width})')
        t0 = time.time(); pool = []
        for sps in rws.SPS_LIBS:
            r = _fit_variant(sps, flux, noise, state['hdr'], width, N, args.n_jobs)
            s = r['sig_boot'].ravel(); pool.append(s[np.isfinite(s)])
            if r['resid_hd'] is not None:
                diag = r['resid_hd']; diag['sps'] = sps
        pool = np.concatenate(pool)
        p16, p50, p84 = np.percentile(pool, [16, 50, 84])
        summaries[key] = dict(p50=float(p50), lo=float(p50-p16), hi=float(p84-p50))
        print(f'     σ_e = {p50:.2f}  −{p50-p16:.2f}/+{p84-p50:.2f} km/s  ({time.time()-t0:.0f}s)\n')

    rws._banner('Hδ masking decision — summary')
    if diag:
        print(f'  Local-MAD diagnostic at Hδ (headline fit, {diag["n_hd_pix"]} pix in band, '
              f'SPS={diag["sps"]}):')
        print(f'    max|resid/noise| at Hδ = {diag["z_max_hd"]:.2f}   '
              f'median = {diag["z_med_hd"]:.2f}')
        print(f'    global 99th pct |resid/noise| = {diag["z_p99_global"]:.2f}   '
              f'→ Hδ {"IS" if diag["z_max_hd"] > diag["z_p99_global"] else "is NOT"} an outlier region\n')
    base = summaries['a_unmasked']['p50']
    print(f'{"variant":12s} {"σ_e (km/s)":>11s} {"Δ vs headline":>14s}')
    print('-' * 42)
    for key, _ in variants:
        s = summaries[key]
        print(f'{key:12s} {s["p50"]:>11.2f} {s["p50"]-base:>+14.2f}')
    print('-' * 42)
    dmax = max(abs(summaries[k]['p50']-base) for k, _ in variants)
    print(f'  max |Δσ_e| from masking Hδ = {dmax:.2f} km/s')
    print(f'  decision: {"KEEP Hδ UNMASKED (a)" if dmax < 3.0 else "Hδ mask matters — review"} '
          f'(threshold 3 km/s ≈ stat error)')

    out = REPO / 'results' / f'sigma_e_hdelta_test_N{N}.npz'
    np.savez(out,
        variants=np.array([v[0] for v in variants]),
        p50=np.array([summaries[v[0]]['p50'] for v in variants]),
        lo=np.array([summaries[v[0]]['lo'] for v in variants]),
        hi=np.array([summaries[v[0]]['hi'] for v in variants]),
        max_delta=dmax, headline_p50=base,
        hdelta_z_max=diag['z_max_hd'] if diag else np.nan,
        global_z_p99=diag['z_p99_global'] if diag else np.nan,
        n_bootstrap=N)
    print(f'\n  → wrote {out}')


if __name__ == '__main__':
    main()
