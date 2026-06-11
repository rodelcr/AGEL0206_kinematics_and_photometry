"""σ_e vs R_e GRID at the WIDE arc-masked window — best-mask headline (2026-06-11).

Supersedes the 4-point D7 driver (`run_sigma_e_Re_systematic_wide.py`). After the
photometry/aperture overhaul we adopt the BEST-MASK (single-Sérsic + color/morph gate
+ WCS registration + companion) curve-of-growth R_e = 2.078" as the headline aperture
(was the expert-mask 2.305"). This driver re-measures σ_e on the headline pipeline
(NEW `_mtwdo_` cube, wR3800_5400_arcmask, He I + M10 masks, Balmer-unmasked), varying
ONLY the aperture r_max over a grid that BRACKETS the new headline, holding the arc
mask, centre, cube, frames, and bootstrap seeds fixed.

R_e grid (light-CoG family + Ca H+K/G kinematic estimator):
  best_F140W   1.876"  (global+companion mask CoG)
  best_mean    2.078"  (global+companion mask CoG)        ← NEW HEADLINE
  best_F200LP  2.281"  (global+companion mask CoG)
  exp_F140W    2.168"  (expert mask CoG)        [reuses resys_F140W cache]
  exp_mean     2.305"  (expert mask CoG)        [reuses new_clean_hei cache]
  exp_F200LP   2.441"  (expert mask CoG)        [reuses resys_F200LP cache]
  CaHK_G       2.902"  (Ca H+K + G-band depth I-map CoG)  [reuses resys_CaHK_G cache]

Headline σ_e = σ_e(best_mean = 2.078").
R_e systematic = peak-to-peak/2 of σ_e across the full grid (the rising σ(R) gradient
to CaHK+G is folded in, per the 2026-06-08 user decision).

Best-mask R_e is computed inline (NOT hardcoded) via the same CoG used by
`scripts/re_mask_sensitivity.py` (global color/morph mask ∪ 2-R_e companion mask).

Usage:  conda activate ISMgas
        python scripts/run_sigma_e_Re_grid.py --n_bootstrap 500 --n_jobs 8
Caches: results/run_wide_sigma_e/resys_best_{F140W,mean,F200LP}/  (+ reused dirs above)
Summary: results/sigma_e_Re_grid_N{N}.npz
"""
from __future__ import annotations
import os, sys, time, argparse
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); os.chdir(REPO)

import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws
from scripts.final_sigma_e import curve_of_growth, find_center, RA_DEFL, DEC_DEFL

HEADLINE_CUBE = REPO / 'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits'
CACHE_ROOT = REPO / 'results' / 'run_wide_sigma_e'
VDI = REPO / '..' / 'velocity_dispersion_from_IFU'
HST_CUT = {"F140W": VDI / "AGEL020613-011417A_F140W_WFC3_cutout_L3.fits",
           "F200LP": VDI / "AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits"}


def best_mask_Re():
    """Best-mask (global color/morph ∪ companion) CoG R_e — same recipe as
    re_mask_sensitivity.py. Returns (Re_F140W, Re_F200LP, mean)."""
    sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
    apm = np.load("results/aperture_2re_masks.npz", allow_pickle=True)
    res = {}
    for n in ("F140W", "F200LP"):
        with fits.open(HST_CUT[n]) as h:
            img = h[0].data.astype(float)
            ps = abs(WCS(h[0].header).proj_plane_pixel_scales()[0].value) * 3600
        mask = (sysz[f"{n}_global_mask"].astype(bool)
                | apm[f"{n}_2Re_companion"].astype(bool))
        w = WCS(fits.getheader(HST_CUT[n]))
        xc, yc = find_center(img, mask, w, RA_DEFL, DEC_DEFL, 3.0, ps)
        res[n] = float(curve_of_growth(img, (xc, yc), ps, mask=mask))
    return res["F140W"], res["F200LP"], 0.5 * (res["F140W"] + res["F200LP"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n_bootstrap', type=int, default=500)
    p.add_argument('--n_jobs', type=int, default=8)
    p.add_argument('--force', action='store_true')
    args = p.parse_args(); N = args.n_bootstrap

    rws._banner(f'σ_e vs R_e GRID @ WIDE arc-masked window  (best-mask headline, N={N})')
    fse.IFU_FILE = str(HEADLINE_CUBE)
    state = fse.load_setup()

    bF140, bF200, bMEAN = best_mask_Re()
    print(f'\n  Best-mask CoG R_e:  F140W={bF140:.3f}"  F200LP={bF200:.3f}"  '
          f'mean={bMEAN:.3f}"  ← NEW HEADLINE')
    print(f'  Expert-mask CoG R_e: F140W={state["Re_140"]:.3f}"  '
          f'F200LP={state["Re_200"]:.3f}"  mean={state["R_E"]:.3f}"')
    print(f'  CaHK+G R_e: {state["Re_cahk"]:.3f}"\n')

    # (label, R_e, cache_dir, is_headline)
    estimators = [
        ('best_F140W',  bF140,            CACHE_ROOT / 'resys_best_F140W', False),
        ('best_mean',   bMEAN,            CACHE_ROOT / 'resys_best_mean',  True),
        ('best_F200LP', bF200,            CACHE_ROOT / 'resys_best_F200LP',False),
        ('exp_F140W',   state['Re_140'],  CACHE_ROOT / 'resys_F140W',      False),
        ('exp_mean',    state['R_E'],     CACHE_ROOT / 'new_clean_hei',    False),
        ('exp_F200LP',  state['Re_200'],  CACHE_ROOT / 'resys_F200LP',     False),
        ('CaHK_G',      state['Re_cahk'], CACHE_ROOT / 'resys_CaHK_G',     False),
    ]

    summaries = {}
    for label, r_e, cache_dir, is_hl in estimators:
        flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
            state, r_max=r_e, mask_weight=0.0)
        n_med = float(np.nanmedian(noise[noise > 0]))
        noise = np.where(noise > 0, noise, n_med * 0.1)
        tag = '  ← HEADLINE' if is_hl else ''
        print(f'  ── R_e[{label}] = {r_e:.3f}"  (n_kept={n_kept}, S/N={sn_band:.2f})'
              f'  cache={cache_dir.name}{tag}')
        t0 = time.time(); results = {}
        for sps in rws.SPS_LIBS:
            results[sps], _ = rws._fit_one_sps(
                sps, flux, noise, state['hdr'], cache_dir, N, args.n_jobs,
                force=args.force, bad_pix_mask=True)
        s = rws._pool_and_summarise(results)
        s['R_e'] = float(r_e); s['n_kept'] = int(n_kept); s['is_headline'] = is_hl
        summaries[label] = s
        print(f'     σ_e = {s["p50"]:.2f}  −{s["stat_lo"]:.2f}/+{s["stat_hi"]:.2f} km/s'
              f'  (SPS spread {s["spread"]:.2f}; {time.time()-t0:.0f}s)\n')

    hl = summaries['best_mean']['p50']
    p50s = {k: summaries[k]['p50'] for k in summaries}
    vals = np.array(list(p50s.values()))
    pk2pk = float(vals.max() - vals.min()); sys_re = pk2pk / 2.0
    # light-CoG-only family (exclude CaHK+G) as a cross-check
    light = [k for k in p50s if k != 'CaHK_G']
    lvals = np.array([p50s[k] for k in light])
    sys_re_light = float(lvals.max() - lvals.min()) / 2.0
    # ADOPTED (user 2026-06-11): best-mask light family only — the three best-mask
    # CoG R_e {best_F140W, best_mean, best_F200LP}. Consistent with "best mask
    # throughout": excludes the superseded expert-mask R_e AND the CaHK+G tracer
    # (a different I-map definition at 2.90", kept as a noted cross-check, NOT an
    # error on the photometric R_e).
    bvals = np.array([p50s[k] for k in ('best_F140W', 'best_mean', 'best_F200LP')])
    sys_re_bestlight = float(bvals.max() - bvals.min()) / 2.0
    cahk_dev = float(p50s['CaHK_G'] - hl)  # CaHK+G cross-check deviation
    # asymmetric about the headline
    dlo = hl - vals.min(); dhi = vals.max() - hl

    rws._banner('σ_e vs R_e grid — summary  (best-mask headline)')
    print(f'{"R_e source":12s} {"R_e [\"]":>8s} {"n_kept":>7s} {"σ_e (km/s)":>11s} {"Δ vs HL":>9s}')
    print('-' * 52)
    for label, _, _, is_hl in estimators:
        s = summaries[label]
        star = ' *' if is_hl else '  '
        print(f'{label:12s} {s["R_e"]:>8.3f} {s["n_kept"]:>7d} {s["p50"]:>11.2f} '
              f'{s["p50"]-hl:>+9.2f}{star}')
    print('-' * 52)
    print(f'  HEADLINE σ_e = σ_e(R_e={summaries["best_mean"]["R_e"]:.3f}") = {hl:.2f} '
          f'−{summaries["best_mean"]["stat_lo"]:.2f}/+{summaries["best_mean"]["stat_hi"]:.2f} km/s')
    print(f'  full-grid peak-to-peak = {pk2pk:.2f}  →  full R_e-source sys = ±{sys_re:.2f} km/s')
    print(f'  (about headline: −{dlo:.2f}/+{dhi:.2f})')
    print(f'  ADOPTED: best-mask light family ±{sys_re_bestlight:.2f} km/s  '
          f'(best F140W/mean/F200LP only)')
    print(f'  light-no-CaHK (incl. expert pts) = ±{sys_re_light:.2f} km/s   [cross-check]')
    print(f'  CaHK+G cross-check deviation = {cahk_dev:+.2f} km/s at 2.902"')
    print(f'  (prior D7 4-pt expert-headline: ±6.13 about 269.62)')

    out = REPO / 'results' / f'sigma_e_Re_grid_N{N}.npz'
    labels = [e[0] for e in estimators]
    np.savez(out,
        labels=np.array(labels),
        R_e=np.array([summaries[k]['R_e'] for k in labels]),
        p50=np.array([summaries[k]['p50'] for k in labels]),
        stat_lo=np.array([summaries[k]['stat_lo'] for k in labels]),
        stat_hi=np.array([summaries[k]['stat_hi'] for k in labels]),
        n_kept=np.array([summaries[k]['n_kept'] for k in labels]),
        is_headline=np.array([summaries[k]['is_headline'] for k in labels]),
        headline_sigma_e=hl,
        headline_R_e=summaries['best_mean']['R_e'],
        headline_stat_lo=summaries['best_mean']['stat_lo'],
        headline_stat_hi=summaries['best_mean']['stat_hi'],
        peak_to_peak=pk2pk, sys_re=sys_re, sys_re_light=sys_re_light,
        sys_re_bestlight=sys_re_bestlight, cahk_dev=cahk_dev,
        d_lo=dlo, d_hi=dhi, n_bootstrap=N)
    print(f'\n  → wrote {out}')


if __name__ == '__main__':
    main()
