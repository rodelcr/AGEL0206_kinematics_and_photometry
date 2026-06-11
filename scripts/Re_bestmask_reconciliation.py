"""Best-mask R_e reconciliation + method systematic (2026-06-11).

Script form of notebook 14, re-run on the BEST mask (single-Sérsic + color/morph
gate + WCS registration; the HST companion mask is empty at 2 R_e). Adopting the
best mask throughout pulls the HST F140W+F200LP curve-of-growth R_e from the
expert-mask 2.305" to 2.097".

Method family (mean of F140W + F200LP), all on the best mask:
  raw CoG (r_max=6")          — headline; photutils-validated to ±0.002"
  sky-sub CoG (r_max=6")      — subtract median sky in a 10-14" arc-free annulus
  single-Sérsic r_eff         — astropy Sersic2D fit (amplitude=I(r_eff) gotcha handled)
R_e method systematic = peak-to-peak/2 across the three.

Output: results/Re_cog_reconciliation_bestmask.npz  (+ printed table)
Usage:  conda activate ISMgas; python scripts/Re_bestmask_reconciliation.py
"""
import os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
import scripts.final_sigma_e as fse
import scripts.measure_Re as mre
from scripts.final_sigma_e import (curve_of_growth, find_center, RA_DEFL, DEC_DEFL,
                                   HST_F140W, HST_F200LP)
from scripts.sersic_total_photometry import fit_sersic2d, R_E_INIT

BANDS = ["F200LP", "F140W"]
IMG = {"F200LP": HST_F200LP, "F140W": HST_F140W}
SYSZ = np.load("results/photometry_systematics.npz", allow_pickle=True)
APM = np.load("results/aperture_2re_masks.npz", allow_pickle=True)


def cog_sky(img, cen, ps, mask, r_max, sky, step=0.08):
    xc, yc = cen; ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]; r = np.hypot(xx - xc, yy - yc) * ps
    edges = np.arange(0, r_max + step, step); rm, Im = [], []
    im = img - sky
    for j in range(len(edges) - 1):
        ann = (r >= edges[j]) & (r < edges[j + 1])
        if mask is not None:
            ann = ann & ~mask
        if ann.sum() == 0:
            continue
        Im.append(float(np.mean(im[ann]))); rm.append(0.5 * (edges[j] + edges[j + 1]))
    return mre.measure_Re_from_profile(np.array(rm), np.array(Im))[0]


def main():
    fam = {}
    for b in BANDS:
        with fits.open(IMG[b]) as h:
            img = h[0].data.astype(float); w = WCS(h[0].header)
        ps = float(np.abs(proj_plane_pixel_scales(w)[0]) * 3600)
        mask = (SYSZ[f"{b}_global_mask"].astype(bool)
                | APM[f"{b}_2Re_companion"].astype(bool))
        xc, yc = find_center(img, mask, w, RA_DEFL, DEC_DEFL, 3.0, ps)
        # raw + sky-sub CoG
        raw6 = float(curve_of_growth(img, (xc, yc), ps, mask=mask, r_max=6.0))
        yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
        r = np.hypot(xx - xc, yy - yc) * ps
        sky = float(np.median(img[(r > 10) & (r < 14) & ~mask]))
        sky6 = float(cog_sky(img, (xc, yc), ps, mask, 6.0, sky))
        # single-Sérsic r_eff on best mask
        fit, _, _ = fit_sersic2d(img, mask, (xc, yc), R_E_INIT, 6.0, ps)
        sreff = float(fit.r_eff.value) * ps
        snn = float(fit.n.value)
        fam[b] = dict(raw=raw6, sky=sky6, sersic=sreff, n=snn, sky_level=sky)
        print(f"  {b}: raw CoG={raw6:.3f}\"  sky-sub CoG={sky6:.3f}\"  "
              f"Sérsic r_eff={sreff:.3f}\" (n={snn:.2f})  sky={sky:.2e}")

    hl_raw = 0.5 * (fam["F200LP"]["raw"] + fam["F140W"]["raw"])
    hl_sky = 0.5 * (fam["F200LP"]["sky"] + fam["F140W"]["sky"])
    hl_ser = 0.5 * (fam["F200LP"]["sersic"] + fam["F140W"]["sersic"])
    family = np.array([hl_raw, hl_sky, hl_ser])
    re_sys = float((family.max() - family.min()) / 2)

    print(f"\n  HST-mean R_e (best mask):  raw CoG = {hl_raw:.3f}\"  [HEADLINE]   "
          f"sky-sub = {hl_sky:.3f}\"   Sérsic = {hl_ser:.3f}\"")
    print(f"  R_e method systematic = ±{re_sys:.3f}\"  (~{100*re_sys/hl_raw:.0f}% of headline)")
    print(f"  (expert-mask headline was 2.305\"; method sys was ±0.082\")")

    np.savez("results/Re_cog_reconciliation_bestmask.npz",
             bands=np.array(BANDS),
             raw=np.array([fam[b]["raw"] for b in BANDS]),
             skysub=np.array([fam[b]["sky"] for b in BANDS]),
             sersic_reff=np.array([fam[b]["sersic"] for b in BANDS]),
             sersic_n=np.array([fam[b]["n"] for b in BANDS]),
             headline_mean_raw=hl_raw, skysub_mean=hl_sky, sersic_mean=hl_ser,
             re_systematic_arcsec=re_sys)
    print("\n  → results/Re_cog_reconciliation_bestmask.npz")


if __name__ == "__main__":
    main()
