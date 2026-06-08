"""σ_e arc-masking systematic — spectroscopic analogue of the ±0.16 dex M★ one.

The photometry-side M★ budget (2026-05-29, `scripts/photometry_systematics.py`)
quoted an explicit ±0.16 dex *masking-approach* systematic: the headline M★
moves when you swap the hand-painted "expert" arc mask for objective /
IR-extended alternatives. This driver is the matching test on the σ_e side:
**reproject each masking approach onto the IFU grid, override the per-spaxel
arc mask used by the I-weighted aperture extraction, and re-run the EXACT
headline σ_e pipeline** (`run_wide_sigma_e._fit_one_sps`, frame-aware ppxf ×
3 SPS × 15 polynomial degrees × N wild-bootstrap at wR3800_5400_arcmask on the
NEW `_mtwdo_` cube with M10 sky + He I masks, Balmer-unmasked).

Only `state['arc_spax_mask']` is varied. R_E (2.305"), the HST-mean centre, the
wavelength-direction goodpixels (arc/sky/CR masks), the cube, the SPS frames,
and the bootstrap seeds are all held FIXED — so the σ_e spread across approaches
is a clean isolate of the spatial-arc-mask choice, just as the M★ systematic
isolated it on the photometry side.

Masking approaches (all live on the 500×500 F200LP grid → reprojected with the
same WCS + order-0 map_coordinates path as the headline arc_spax_mask):
  expert   — F200LP_expert      (hand-painted; == current headline mask)
  sersic   — F200LP_m_sersic    (objective Sérsic-residual auto reproduction)
  perband  — F200LP_perband_mask(F200LP-located, per-band recipe)
  global   — F200LP_global_mask (IR-extended union → largest footprint)

Sources: results/arc_mask_verification.npz (expert, sersic) and
results/photometry_systematics.npz (perband, global).

Usage
─────
  conda activate ISMgas
  # smoke (validate wiring + reprojected spaxel counts)
  python scripts/run_sigma_e_mask_systematic.py --n_bootstrap 50
  # production
  python scripts/run_sigma_e_mask_systematic.py --n_bootstrap 500 --n_jobs 8

Caches: results/run_wide_sigma_e/masksys_{approach}/wR3800_5400_arcmask_{sps}_T*_N{N}.npz
Summary npz: results/sigma_e_mask_systematic_N{N}.npz
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
os.chdir(REPO)

from astropy.io import fits                                  # noqa: E402
from astropy.wcs import WCS                                  # noqa: E402
from scipy.ndimage import map_coordinates                    # noqa: E402

import scripts.final_sigma_e as fse                          # noqa: E402
import scripts.run_wide_sigma_e as rws                       # noqa: E402

# NEW _mtwdo_ headline cube + its clean (bad-pix + Balmer-unmasked) path.
HEADLINE_CUBE = REPO / 'raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits'
CACHE_ROOT = REPO / 'results' / 'run_wide_sigma_e'

# (approach key, npz file, key in npz, human description)
APPROACHES = [
    ('expert',  'results/arc_mask_verification.npz',  'F200LP_expert',
     'hand-painted F200LP mask (== current headline arc_spax_mask)'),
    ('sersic',  'results/arc_mask_verification.npz',  'F200LP_m_sersic',
     'objective Sérsic-residual auto reproduction (k=3)'),
    ('perband', 'results/photometry_systematics.npz', 'F200LP_perband_mask',
     'F200LP-located per-band recipe'),
    ('global',  'results/photometry_systematics.npz', 'F200LP_global_mask',
     'IR-extended union (largest footprint)'),
]


def _reproject_f200_to_ifu(mask_f200, wcs_f200, wcs_ifu, ny, nx):
    """Same path as final_sigma_e.load_setup §4: IFU spaxel centres → world →
    F200LP pixels → nearest-neighbour sample of the (boolean) mask."""
    yy_, xx_ = np.mgrid[:ny, :nx]
    ra_s, dec_s = wcs_ifu.pixel_to_world_values(xx_.ravel(), yy_.ravel())
    xh, yh = wcs_f200.world_to_pixel_values(ra_s, dec_s)
    spax = map_coordinates(
        mask_f200.astype(float), [yh, xh], order=0, mode='constant', cval=0.0
    ).reshape(ny, nx).astype(bool)
    return spax


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--n_bootstrap', type=int, default=500)
    p.add_argument('--n_jobs', type=int, default=8)
    p.add_argument('--force', action='store_true')
    args = p.parse_args()
    N = args.n_bootstrap

    rws._banner(f'σ_e ARC-MASKING SYSTEMATIC  (N={N})  — spectroscopic analogue of the ±0.16 dex M★ one')

    # 1. Load the headline state ONCE (expert mask drives centre + R_E; both held fixed).
    fse.IFU_FILE = str(HEADLINE_CUBE)
    state = fse.load_setup()
    ny, nx = state['ny'], state['nx']
    R_E = state['R_E']
    wcs_ifu = WCS(state['hdr'], naxis=2)
    wcs_f200 = WCS(fits.getheader(fse.HST_F200LP))
    expert_loaded = state['arc_spax_mask'].copy()
    print(f'\n  R_E held fixed at {R_E:.4f}"   |   headline arc_spax_mask = {expert_loaded.sum()} spaxels')
    print(f'  cube: {HEADLINE_CUBE.name}   |   bad_pix_mask=True (Balmer-unmasked, M10+HeI goodpixels)\n')

    summaries = {}
    spaxel_counts = {}
    for key, npz, mkey, descr in APPROACHES:
        d = np.load(REPO / npz, allow_pickle=True)
        mask_f200 = np.asarray(d[mkey]).astype(bool)
        arc_spax = _reproject_f200_to_ifu(mask_f200, wcs_f200, wcs_ifu, ny, nx)
        spaxel_counts[key] = int(arc_spax.sum())

        # Sanity: 'expert' reprojection must reproduce the headline mask bit-for-bit.
        if key == 'expert':
            match = bool(np.array_equal(arc_spax, expert_loaded))
            print(f'  [expert] reprojected mask matches headline arc_spax_mask: {match} '
                  f'({arc_spax.sum()} vs {expert_loaded.sum()} spaxels)')

        # Override only the spatial arc mask, re-extract the aperture spectrum.
        st = dict(state)
        st['arc_spax_mask'] = arc_spax
        flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
            st, r_max=R_E, mask_weight=0.0,
        )
        n_med = float(np.nanmedian(noise[noise > 0]))
        noise = np.where(noise > 0, noise, n_med * 0.1)

        cache_dir = CACHE_ROOT / f'masksys_{key}'
        print(f'\n  ── approach "{key}": {descr}')
        print(f'     F200LP mask px = {mask_f200.sum():6d} → IFU spaxels flagged = {arc_spax.sum():3d}'
              f'   |   aperture n_kept={n_kept}  S/N={sn_band:.2f}')
        t0 = time.time()
        results = {}
        for sps in rws.SPS_LIBS:
            results[sps], _ = rws._fit_one_sps(
                sps, flux, noise, state['hdr'], cache_dir, N, args.n_jobs,
                force=args.force, bad_pix_mask=True,
            )
        summ = rws._pool_and_summarise(results)
        summ['n_spaxels_arc'] = int(arc_spax.sum())
        summ['n_kept'] = int(n_kept)
        summaries[key] = summ
        print(f'     σ_e(<R_e) = {summ["p50"]:.2f}  −{summ["stat_lo"]:.2f}/+{summ["stat_hi"]:.2f} km/s'
              f'   (SPS spread {summ["spread"]:.2f}; {time.time()-t0:.0f}s)')

    # 2. Systematic = peak-to-peak / 2 across approaches (matches M★ convention).
    p50s = {k: summaries[k]['p50'] for k in summaries}
    vals = np.array(list(p50s.values()))
    pk2pk = float(vals.max() - vals.min())
    sys_mask = pk2pk / 2.0
    # also the directed shifts vs expert (the headline)
    expert_p50 = p50s['expert']

    rws._banner('σ_e MASKING-APPROACH SYSTEMATIC — summary')
    print(f'{"approach":10s} {"arc spaxels":>11s} {"n_kept":>7s} {"σ_e (km/s)":>11s} '
          f'{"stat asym":>14s} {"Δ vs expert":>12s}')
    print('-' * 74)
    for key, _, _, _ in APPROACHES:
        s = summaries[key]
        print(f'{key:10s} {s["n_spaxels_arc"]:>11d} {s["n_kept"]:>7d} {s["p50"]:>11.2f} '
              f'  −{s["stat_lo"]:>5.2f}/+{s["stat_hi"]:>5.2f} {s["p50"]-expert_p50:>+12.2f}')
    print('-' * 74)
    print(f'  peak-to-peak               = {pk2pk:.2f} km/s')
    print(f'  σ_e masking systematic     = peak-to-peak/2 = ±{sys_mask:.2f} km/s')
    print(f'  (headline = expert approach = {expert_p50:.2f} km/s)')

    out = REPO / 'results' / f'sigma_e_mask_systematic_N{N}.npz'
    np.savez(
        out,
        approaches=np.array([a[0] for a in APPROACHES]),
        p50=np.array([summaries[a[0]]['p50'] for a in APPROACHES]),
        stat_lo=np.array([summaries[a[0]]['stat_lo'] for a in APPROACHES]),
        stat_hi=np.array([summaries[a[0]]['stat_hi'] for a in APPROACHES]),
        n_spaxels_arc=np.array([summaries[a[0]]['n_spaxels_arc'] for a in APPROACHES]),
        n_kept=np.array([summaries[a[0]]['n_kept'] for a in APPROACHES]),
        per_sps=np.array([summaries[a[0]]['per_sps'] for a in APPROACHES], dtype=object),
        peak_to_peak=pk2pk, sys_mask=sys_mask, expert_p50=expert_p50,
        n_bootstrap=N, window=rws.WLABEL,
    )
    print(f'\n  → wrote {out}')
    return summaries


if __name__ == '__main__':
    main()
