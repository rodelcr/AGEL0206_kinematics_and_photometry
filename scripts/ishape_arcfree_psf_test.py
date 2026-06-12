"""ARC-FREE, PSF-matched I-weight maps → σ_e (2026-06-12, user proposal).

`ishape_psf_convolved_test.py` convolved the RAW HST images, which still contain the arc, so
the +7-8 km/s it found was arc LEAKAGE, not a clean resolution effect. The correct
resolution-matched deflector weight removes the arc FIRST, then convolves to the KCWI PSF:
  (a) Sérsic-conv:  the single-Sérsic deflector MODEL (arc-free by construction) convolved.
  (b) filled-conv:  the HST image with the masked arc region FILLED by the Sérsic model,
                    then convolved.
Both give a smooth, arc-free, hole-free deflector light map at the KCWI ~1.27" resolution —
candidate HEADLINE weightings. Compared against the IFU white-light headline and the raw /
raw-convolved HST maps.

Output: results/ishape_arcfree_psf_test.npz (+ printed table)
Usage: conda activate ISMgas; python scripts/ishape_arcfree_psf_test.py
"""
import os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import gaussian_filter

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, str(REPO)); os.chdir(REPO)
import scripts.final_sigma_e as fse
from scripts.mask_method_comparison import load_band
from scripts.photometry_systematics import fit_sersic
from scripts.ishape_psf_convolved_test import reproject, sigma_e_for_weight, SEEING_FWHM

SYSZ = np.load("results/photometry_systematics.npz", allow_pickle=True)


def main():
    fse.IFU_FILE = "raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits"
    st = fse.load_setup()
    wcs_ifu = WCS(st["hdr"], naxis=2); ny, nx = st["ny"], st["nx"]

    maps = {"IFU white-light (headline)": st["ifu_band"]}
    for band in ("F140W", "F200LP"):
        b = load_band(band); ps = b["pix_scale"]; w = b["wcs"]
        mask = SYSZ[f"{band}_global_mask"].astype(bool)
        sig_pix = (SEEING_FWHM / 2.3548) / ps
        # single-Sérsic deflector model over the cutout (arc-free: fit on masked data)
        model, sky = fit_sersic(b, mask)
        data = np.nan_to_num(b["img"]) - sky
        filled = np.where(mask, model, data)          # arc region replaced by the model
        # (a) Sérsic-conv, (b) filled-conv — convolve to KCWI PSF on the HST grid, then reproject
        maps[f"{band} Sérsic→PSF-conv"] = reproject(gaussian_filter(model, sig_pix, mode="nearest"),
                                                    w, wcs_ifu, ny, nx)
        maps[f"{band} filled→PSF-conv"] = reproject(gaussian_filter(np.clip(filled, 0, None), sig_pix,
                                                    mode="nearest"), w, wcs_ifu, ny, nx)
        # raw + raw-conv (the confounded one) for reference
        maps[f"{band} raw (arc-contaminated conv)"] = reproject(np.nan_to_num(b["img"]), w, wcs_ifu, ny, nx)
        print(f"  {band}: σ_psf {sig_pix:.1f} pix; Sérsic n={'(fit)':>5s}")

    print(f"\n  {'I-weight map':36s} {'σ_e (FSPS)':>10s} {'Δ vs headline':>14s}")
    print("  " + "-" * 64)
    res, hl = {}, None
    for name, Imap in maps.items():
        s, nk, sn = sigma_e_for_weight(st, Imap)
        if hl is None:
            hl = s
        res[name] = s
        print(f"  {name:36s} {s:>10.2f} {s-hl:>+13.2f}")

    print("\n  ── ARC-FREE PSF-matched shift vs headline (the clean test) ──")
    for band in ("F140W", "F200LP"):
        for kind in ("Sérsic→PSF-conv", "filled→PSF-conv"):
            print(f"    {band} {kind:18s}: {res[f'{band} {kind}']:.2f}  (Δ headline {res[f'{band} {kind}']-hl:+.2f})")
    print(f"  (headline IFU white-light σ_e[FSPS] = {hl:.2f}; pooled headline = 267.31 ± 12.79; I-shape ±2.27)")

    np.savez("results/ishape_arcfree_psf_test.npz",
             names=np.array(list(res.keys())), sigma_e=np.array(list(res.values())),
             headline_fsps=hl, seeing_fwhm=SEEING_FWHM)
    print("\n  → results/ishape_arcfree_psf_test.npz")


if __name__ == "__main__":
    main()
