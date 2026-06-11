"""Masked-spaxel fraction within R_e: direct HST-reprojected mask vs KCWI-PSF-convolved
(2026-06-12).

The headline σ_e arc mask is the F200LP HST mask reprojected (nearest-neighbour) onto the
IFU grid → it flags only the spaxels that DIRECTLY overlap HST arc pixels. But KCWI has
~1.27″ seeing (FWHM, GUIDFWHM), so the arc light bleeds into neighbouring spaxels. This
quantifies how many more spaxels (and how much more I-weight) inside R_e would be flagged
if the mask were convolved with the KCWI PSF before thresholding — i.e. the seeing-smeared
arc footprint vs the geometric one.

Output: results/psf_mask_fraction.npz (+ printed table)
Usage: conda activate ISMgas; python scripts/psf_mask_fraction.py
"""
import os, sys
import numpy as np
from scipy.ndimage import gaussian_filter, binary_dilation

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
import scripts.final_sigma_e as fse

SEEING_FWHM = 1.27          # arcsec (KCWI GUIDFWHM)


def main():
    fse.IFU_FILE = "raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits"
    st = fse.load_setup()
    ps = st["pix_scale_ifu"]                       # arcsec/spaxel (~0.30)
    R_E = 2.097                                     # pinned headline best-mask R_e (load_setup's
                                                    # st['R_E']=2.305 is the expert-mask value)
    r = st["r_spax"]                               # radius map (arcsec)
    mask = st["arc_spax_mask"].astype(bool)        # direct HST-reprojected F200LP mask
    Iw = st["ifu_band"]                            # white-light I-weight map

    in_re = r < R_E
    N = int(in_re.sum())
    Itot = float(Iw[in_re].sum())

    # direct (geometric) mask
    d_cnt = int((mask & in_re).sum())
    d_Iw = float(Iw[mask & in_re].sum())

    # KCWI-PSF-convolved mask: spread the binary mask by the seeing Gaussian, then a spaxel
    # is "PSF-contaminated" if a fraction > thr of its PSF-weighted light comes from the
    # masked (arc) region. sigma in spaxels:
    sig = (SEEING_FWHM / ps) / 2.3548
    conv = gaussian_filter(mask.astype(float), sig, mode="constant")

    print(f"  IFU pix scale = {ps:.3f}\"/spax;  seeing FWHM = {SEEING_FWHM}\" = {SEEING_FWHM/ps:.2f} spax "
          f"(σ_psf = {sig:.2f} spax)")
    print(f"  R_e = {R_E:.3f}\";  spaxels within R_e: N = {N}\n")
    print(f"  {'mask':30s} {'spaxels':>9s} {'% by count':>11s} {'% by I-weight':>14s}")
    print(f"  {'direct (HST reproject)':30s} {d_cnt:>9d} {100*d_cnt/N:>10.1f}% {100*d_Iw/Itot:>13.1f}%")
    out = {"N_in_Re": N, "R_e": R_E, "ps": ps, "seeing_fwhm": SEEING_FWHM,
           "direct_count": d_cnt, "direct_frac_count": d_cnt / N,
           "direct_frac_Iw": d_Iw / Itot}

    # --- cleaner "grow the mask by the seeing": binary dilation by the PSF radius ---
    yy, xx = np.mgrid[-6:7, -6:7]
    for rad_as, lab in ((SEEING_FWHM / 2, "FWHM/2 = %.2f\"" % (SEEING_FWHM / 2)),
                        (SEEING_FWHM, "FWHM = %.2f\"" % SEEING_FWHM)):
        rad = rad_as / ps
        se = (xx**2 + yy**2) <= rad**2
        dil = binary_dilation(mask, structure=se) & in_re
        c = int(dil.sum()); iw = float(Iw[dil].sum())
        print(f"  {'dilate by '+lab:30s} {c:>9d} {100*c/N:>10.1f}% {100*iw/Itot:>13.1f}%")
        out[f"dilate_count_{rad_as:.2f}"] = c
        out[f"dilate_frac_count_{rad_as:.2f}"] = c / N
        out[f"dilate_frac_Iw_{rad_as:.2f}"] = iw / Itot
    for thr in (0.5, 0.25, 0.1):
        pm = (conv > thr) & in_re
        c = int(pm.sum()); iw = float(Iw[pm].sum())
        print(f"  {'PSF-conv (>%.2f arc light)'%thr:28s} {c:>9d} {100*c/N:>10.1f}% {100*iw/Itot:>13.1f}%")
        out[f"psf_count_thr{thr}"] = c
        out[f"psf_frac_count_thr{thr}"] = c / N
        out[f"psf_frac_Iw_thr{thr}"] = iw / Itot
    print(f"\n  Takeaway: the direct geometric mask flags {100*out['direct_frac_count']:.0f}% by count "
          f"({100*out['direct_frac_Iw']:.0f}% by I-weight). Accounting for the {SEEING_FWHM}\" seeing,"
          f" the arc contaminates much more: dilating by FWHM/2 → {100*out['dilate_frac_count_%.2f'%(SEEING_FWHM/2)]:.0f}%/"
          f"{100*out['dilate_frac_Iw_%.2f'%(SEEING_FWHM/2)]:.0f}%, by full FWHM → "
          f"{100*out['dilate_frac_count_%.2f'%SEEING_FWHM]:.0f}%/{100*out['dilate_frac_Iw_%.2f'%SEEING_FWHM]:.0f}%."
          f" (A hard PSF-grown mask removes nearly the whole aperture → infeasible; hence I-weighting"
          f" + the mask-weight systematic instead.) The convolve>0.5 cut (1%) is a thin-arc dilution"
          f" artifact, not 'less masking'.")
    np.savez("results/psf_mask_fraction.npz", **out)
    print("  → results/psf_mask_fraction.npz")


if __name__ == "__main__":
    main()
