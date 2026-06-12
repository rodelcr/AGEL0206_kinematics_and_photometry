"""Fold the arc-free PSF-matched weight maps into the I-shape systematic — POOLED footing.

The I-shape term (±2.27) is peak-to-peak/2 of the 10-shape σ_e spread 266.83–271.37, where
each shape's σ_e is the 3-SPS × 15-deg POOLED central. ishape_arcfree_psf_test.py measured the
4 arc-free maps at single-FSPS only (which runs ~2 km/s high), so this re-measures them on the
pooled footing and recomputes the I-shape spread over all 14 shapes.

Output: results/ishape_arcfree_pooled.npz (+ printed new I-shape value)
Usage: conda activate ISMgas; python scripts/ishape_arcfree_pooled.py
"""
import os, sys
import numpy as np
from astropy.wcs import WCS
from ppxf.ppxf import ppxf

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, str(REPO)); os.chdir(REPO)
import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum
from scripts.ishape_arcfree_psf_test import main as _unused  # ensures module imports cleanly

# existing 10-shape pooled spread (M11)
EXIST_MIN, EXIST_MAX = 266.83, 271.37
R_E = 2.097


def pooled_sigma_e(state, Imap):
    """3-SPS × 15-deg pooled median σ_e for a given I-weight map (σ is frame-invariant)."""
    st = dict(state); st["ifu_band"] = np.clip(Imap, 0, None)
    flux, noise, n_kept, sn = fse.extract_aperture_spectrum(st, r_max=R_E, mask_weight=0.0)
    n_med = float(np.nanmedian(noise[noise > 0]))
    noise = np.where(noise > 0, noise, n_med * 0.1)
    sig_all = []
    for sps in rws.SPS_LIBS:
        lo, hi = rws.LAM_OBS_RANGE
        lr = (max(lo, rws.sps_safe_obs_min(sps)), hi)
        inp = setup_ppxf_inputs_from_spectrum(flux, noise, state["hdr"], sps_name=sps,
            z=fse.Z_SYSTEMIC, lam_obs_range=lr, lam_fit_range=lr,
            lam_range_temp=rws.SPS_LAM_RANGE_TEMP[sps], verbose=False,
            frame_galaxy='auto', mask_balmer=False)
        gp = rws._apply_arc_mask(inp['goodpixels'], inp['lam_gal_rest']) if rws.ARC_MASK else inp['goodpixels']
        gp = rws._apply_bad_pixels_mask(gp, inp['lam_gal_rest'])
        for deg in rws.DEGREES:
            pp = ppxf(inp['sps'].templates, inp['galaxy'], inp['noise'], inp['velscale'],
                      inp['start'], goodpixels=gp, plot=False, moments=2, trig=False,
                      degree=int(deg), mdegree=0, lam=inp['lam_gal_rest'],
                      lam_temp=inp['lam_temp'], quiet=True)
            sig_all.append(pp.sol[1])
    return float(np.median(sig_all))


def main():
    from scripts.ishape_arcfree_psf_test import SEEING_FWHM, SYSZ, reproject
    from scripts.mask_method_comparison import load_band
    from scripts.photometry_systematics import fit_sersic
    from scipy.ndimage import gaussian_filter

    fse.IFU_FILE = "raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits"
    st = fse.load_setup()
    wcs_ifu = WCS(st["hdr"], naxis=2); ny, nx = st["ny"], st["nx"]

    maps = {"IFU white-light (headline, calib)": st["ifu_band"]}
    for band in ("F140W", "F200LP"):
        b = load_band(band); ps = b["pix_scale"]; w = b["wcs"]
        mask = SYSZ[f"{band}_global_mask"].astype(bool)
        sig_pix = (SEEING_FWHM / 2.3548) / ps
        model, sky = fit_sersic(b, mask)
        filled = np.where(mask, model, np.nan_to_num(b["img"]) - sky)
        maps[f"{band} Sérsic→PSF-conv"] = reproject(gaussian_filter(model, sig_pix, mode="nearest"),
                                                    w, wcs_ifu, ny, nx)
        maps[f"{band} filled→PSF-conv"] = reproject(gaussian_filter(np.clip(filled, 0, None), sig_pix,
                                                    mode="nearest"), w, wcs_ifu, ny, nx)

    print(f"\n  {'I-weight map':36s} {'σ_e POOLED':>11s}")
    print("  " + "-" * 50)
    pooled = {}
    for name, Imap in maps.items():
        s = pooled_sigma_e(st, Imap); pooled[name] = s
        print(f"  {name:36s} {s:>11.2f}")

    calib = pooled["IFU white-light (headline, calib)"]
    arcfree = [pooled[k] for k in pooled if "PSF-conv" in k]
    new_min = min(EXIST_MIN, *arcfree); new_max = max(EXIST_MAX, *arcfree)
    new_ishape = (new_max - new_min) / 2
    print(f"\n  calibration: pooled IFU white-light = {calib:.2f}  (pooled headline = 267.31; "
          f"offset {calib-267.31:+.2f})")
    print(f"  arc-free pooled range: {min(arcfree):.2f}–{max(arcfree):.2f}")
    print(f"  existing 10-shape spread: {EXIST_MIN}–{EXIST_MAX}  (I-shape ±2.27)")
    print(f"  NEW 14-shape spread: {new_min:.2f}–{new_max:.2f}  →  I-shape ±{new_ishape:.2f}")
    if new_max <= EXIST_MAX + 0.01 and new_min >= EXIST_MIN - 0.01:
        print(f"  → arc-free maps fall WITHIN the existing spread; I-shape UNCHANGED at ±2.27.")
    np.savez("results/ishape_arcfree_pooled.npz",
             names=np.array(list(pooled.keys())), sigma_e=np.array(list(pooled.values())),
             new_min=new_min, new_max=new_max, new_ishape=new_ishape, calib=calib)
    print("\n  → results/ishape_arcfree_pooled.npz")


if __name__ == "__main__":
    main()
