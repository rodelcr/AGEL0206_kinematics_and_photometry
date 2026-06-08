"""D7 R_e-source systematic at the WIDE arc-masked window (Task 5, 2026-06-08).

D7 (TESTS_AND_DIAGNOSTICS) quantifies how σ_e moves when the aperture radius R_e
is taken from different estimators. It was measured at the NARROW window
(spread = 16.9 km/s: mean=268, F140W=265, F200LP=276, CaHK+G=282) and **never
re-measured at the wide arc-masked window** that drives the headline — so it is
NOT folded into the M11 wide budget. This driver re-measures it on the headline
pipeline (NEW `_mtwdo_` cube, wR3800_5400_arcmask, He I + M10 masks,
Balmer-unmasked), varying ONLY the aperture r_max while holding the arc mask,
centre, cube, frames, and bootstrap seeds fixed.

R_e estimators (from final_sigma_e.load_setup on the headline cube):
  mean    2.305"  — headline (F140W+F200LP CoG mean)   [reuses new_clean_hei cache]
  F140W   2.168"
  F200LP  2.441"
  CaHK_G  2.902"  — Ca H+K + G-band absorption-depth I-map CoG

The wide-window D7 systematic = peak-to-peak/2 across these σ_e(R_e) values.

Usage:  conda activate ISMgas
        python scripts/run_sigma_e_Re_systematic_wide.py --n_bootstrap 500 --n_jobs 8
Caches: results/run_wide_sigma_e/resys_{label}/...   (mean reuses new_clean_hei/)
Summary: results/sigma_e_Re_systematic_wide_N{N}.npz
"""
from __future__ import annotations
import os, sys, time, argparse
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); os.chdir(REPO)

import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws

HEADLINE_CUBE = REPO / 'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits'
CACHE_ROOT = REPO / 'results' / 'run_wide_sigma_e'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n_bootstrap', type=int, default=500)
    p.add_argument('--n_jobs', type=int, default=8)
    p.add_argument('--force', action='store_true')
    args = p.parse_args(); N = args.n_bootstrap

    rws._banner(f'D7 R_e-SOURCE SYSTEMATIC @ WIDE arc-masked window  (N={N})')
    fse.IFU_FILE = str(HEADLINE_CUBE)
    state = fse.load_setup()

    # R_e estimators straight off the headline load_setup.
    estimators = [
        ('mean',   state['R_E'],     CACHE_ROOT / 'new_clean_hei'),  # reuse headline cache
        ('F140W',  state['Re_140'],  CACHE_ROOT / 'resys_F140W'),
        ('F200LP', state['Re_200'],  CACHE_ROOT / 'resys_F200LP'),
        ('CaHK_G', state['Re_cahk'], CACHE_ROOT / 'resys_CaHK_G'),
    ]
    print(f'\n  R_e estimators:  mean={state["R_E"]:.3f}"  F140W={state["Re_140"]:.3f}"  '
          f'F200LP={state["Re_200"]:.3f}"  CaHK+G={state["Re_cahk"]:.3f}"\n')

    summaries = {}
    for label, r_e, cache_dir in estimators:
        flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
            state, r_max=r_e, mask_weight=0.0)
        n_med = float(np.nanmedian(noise[noise > 0]))
        noise = np.where(noise > 0, noise, n_med * 0.1)
        print(f'  ── R_e[{label}] = {r_e:.3f}"  (n_kept={n_kept}, S/N={sn_band:.2f})  cache={cache_dir.name}')
        t0 = time.time(); results = {}
        for sps in rws.SPS_LIBS:
            results[sps], _ = rws._fit_one_sps(
                sps, flux, noise, state['hdr'], cache_dir, N, args.n_jobs,
                force=args.force, bad_pix_mask=True)
        s = rws._pool_and_summarise(results)
        s['R_e'] = float(r_e); s['n_kept'] = int(n_kept)
        summaries[label] = s
        print(f'     σ_e = {s["p50"]:.2f}  −{s["stat_lo"]:.2f}/+{s["stat_hi"]:.2f} km/s'
              f'  (SPS spread {s["spread"]:.2f}; {time.time()-t0:.0f}s)\n')

    p50s = {k: summaries[k]['p50'] for k in summaries}
    vals = np.array(list(p50s.values()))
    pk2pk = float(vals.max() - vals.min()); sys_re = pk2pk / 2.0
    mean_p50 = p50s['mean']

    rws._banner('D7 wide-window R_e-source systematic — summary')
    print(f'{"R_e source":10s} {"R_e [\"]":>8s} {"n_kept":>7s} {"σ_e (km/s)":>11s} {"Δ vs mean":>10s}')
    print('-' * 56)
    for label, _, _ in [(e[0], e[1], e[2]) for e in estimators]:
        s = summaries[label]
        print(f'{label:10s} {s["R_e"]:>8.3f} {s["n_kept"]:>7d} {s["p50"]:>11.2f} {s["p50"]-mean_p50:>+10.2f}')
    print('-' * 56)
    print(f'  peak-to-peak = {pk2pk:.2f} km/s   →   D7(wide) = ±{sys_re:.2f} km/s')
    print(f'  (compare D7(narrow) = ±8.45 km/s [spread 16.9/2]; headline = mean = {mean_p50:.2f})')

    out = REPO / 'results' / f'sigma_e_Re_systematic_wide_N{N}.npz'
    np.savez(out,
        labels=np.array([e[0] for e in estimators]),
        R_e=np.array([summaries[e[0]]['R_e'] for e in estimators]),
        p50=np.array([summaries[e[0]]['p50'] for e in estimators]),
        stat_lo=np.array([summaries[e[0]]['stat_lo'] for e in estimators]),
        stat_hi=np.array([summaries[e[0]]['stat_hi'] for e in estimators]),
        n_kept=np.array([summaries[e[0]]['n_kept'] for e in estimators]),
        peak_to_peak=pk2pk, sys_re=sys_re, mean_p50=mean_p50, n_bootstrap=N)
    print(f'\n  → wrote {out}')


if __name__ == '__main__':
    main()
