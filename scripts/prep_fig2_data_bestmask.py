"""Rebuild figure2_data.npz for the best-mask R_e=2.097" aperture (2026-06-11).

Figure 2 (left inset) of the ApJL shows the I-weighted aperture spectrum + the 45
ppxf best fits (3 SPS × 15 deg) + the pooled σ_e posterior. After "best mask
throughout" the headline aperture is r<R_e=2.097" (was 2.305"). This reconstructs
the self-contained npz the figure loads, pulling the galaxy/goodpixels from the SAME
ppxf setup `run_wide_sigma_e._fit_one_sps` uses, and the best-fits + σ-pool from the
already-computed `resys_best_mean` N=500 caches.

Output: results/run_wide_sigma_e/resys_best_mean/figure2_data.npz
Usage:  conda activate ISMgas; python scripts/prep_fig2_data_bestmask.py
"""
import os, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); os.chdir(REPO)

import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws
from scripts.run_sigma_e_Re_grid import best_mask_Re, HEADLINE_CUBE
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum
from scripts.run_window_sweep import ARC_MASKS_REST, BAD_PIXELS_REST

CACHE = REPO / 'results' / 'run_wide_sigma_e' / 'resys_best_mean'
SPS = ('fsps', 'emiles', 'xsl')


def main():
    fse.IFU_FILE = str(HEADLINE_CUBE)
    state = fse.load_setup()
    _, _, R_E = best_mask_Re()
    print(f'  best-mask R_e = {R_E:.3f}"')
    flux, noise, n_kept, sn_band = fse.extract_aperture_spectrum(
        state, r_max=R_E, mask_weight=0.0)
    n_med = float(np.nanmedian(noise[noise > 0]))
    noise = np.where(noise > 0, noise, n_med * 0.1)

    # galaxy/lam/goodpixels from the SAME setup _fit_one_sps uses (FSPS representative,
    # matching the original single-spectrum figure2_data; bad_pix_mask=True → Balmer kept)
    sps0 = 'fsps'
    lo, hi = rws.LAM_OBS_RANGE
    lr_eff = (max(lo, rws.sps_safe_obs_min(sps0)), hi)
    inp = setup_ppxf_inputs_from_spectrum(
        flux, noise, state['hdr'], sps_name=sps0, z=fse.Z_SYSTEMIC,
        lam_obs_range=lr_eff, lam_fit_range=lr_eff,
        lam_range_temp=rws.SPS_LAM_RANGE_TEMP[sps0],
        verbose=False, frame_galaxy='auto', mask_balmer=False)
    gp = inp['goodpixels']
    if rws.ARC_MASK:
        gp = rws._apply_arc_mask(gp, inp['lam_gal_rest'])
    gp = rws._apply_bad_pixels_mask(gp, inp['lam_gal_rest'])

    galaxy = np.asarray(inp['galaxy'])
    lam_gal_rest = np.asarray(inp['lam_gal_rest'])
    goodpixels = np.asarray(gp)
    n_pix = len(galaxy)
    print(f'  galaxy n_pix={n_pix}  goodpixels={len(goodpixels)}  (orig fig2: 2574/2170)')

    # best-fits + σ pool from the resys_best_mean caches
    all_bf, pool, per_sps = [], [], {}
    for sps in SPS:
        c = sorted(CACHE.glob(f'wR3800_5400_arcmask_{sps}_T*_N500.npz'))
        assert len(c) == 1, f'{sps}: {c}'
        d = np.load(c[0], allow_pickle=True)
        bf = d['bf']
        if bf.shape[1] != n_pix:                       # frame/grid guard
            print(f'  ! {sps} bf n_pix {bf.shape[1]} != {n_pix}; trimming/padding to match')
            m = min(bf.shape[1], n_pix)
            bf2 = np.zeros((bf.shape[0], n_pix)); bf2[:, :m] = bf[:, :m]; bf = bf2
        all_bf.append(bf)
        s = d['sig_boot'].ravel(); s = s[np.isfinite(s)]
        pool.append(s); per_sps[sps] = float(np.percentile(s, 50))
    all_best_fits = np.vstack(all_bf)                  # (45, n_pix)
    sigma_pool = np.concatenate(pool)
    p16, p50, p84 = (float(x) for x in np.percentile(sigma_pool, [16, 50, 84]))
    print(f'  σ_e pool: {p50:.2f} −{p50-p16:.2f}/+{p84-p50:.2f}  ({len(sigma_pool)} samples)')
    print(f'  per-SPS p50: ' + '  '.join(f'{k}={v:.1f}' for k, v in per_sps.items()))

    out = CACHE / 'figure2_data.npz'
    np.savez(out,
        galaxy=galaxy, lam_gal_rest=lam_gal_rest, goodpixels=goodpixels,
        all_best_fits=all_best_fits, sigma_pool=sigma_pool,
        sigma_p16=p16, sigma_p50=p50, sigma_p84=p84,
        per_sps_p50=per_sps, n_pix=n_pix, n_kept=int(n_kept), sn_band=float(sn_band),
        cube_path=str(HEADLINE_CUBE), lam_obs_range=np.array(lr_eff),
        arc_masks_rest=np.array(ARC_MASKS_REST),
        bad_pixels_rest=np.array(BAD_PIXELS_REST),
        mask_balmer=False, R_e=R_E)
    print(f'\n  → {out}')


if __name__ == '__main__':
    main()
