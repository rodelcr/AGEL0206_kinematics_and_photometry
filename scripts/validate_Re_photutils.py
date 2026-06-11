"""Validate our R_e / curve-of-growth machinery against photutils (2026-06-11).

Our pipeline (scripts/final_sigma_e.curve_of_growth + measure_Re.measure_Re_from_profile):
  centre = photutils centroid_2dg  →  azimuthally-averaged mean intensity I(r) in
  0.08" annuli (masked pixels DROPPED from each annulus' mean)  →  integrate
  I(r)·2πr dr  →  half-light radius by interpolation.

The azimuthal-mean step effectively FILLS masked azimuthal wedges with the annulus
mean (assumes the galaxy is azimuthally symmetric) — the standard treatment for a
half-light radius when a mask removes wedges of an otherwise-symmetric profile.

Two INDEPENDENT photutils paths cross-check this:
  P1  CurveOfGrowth  — direct 2D circular-aperture SUM F(<r). With no mask this is
      the ground-truth enclosed flux; it must agree with our integration. With a
      mask photutils DROPS masked pixels from the sum (does NOT fill) → expected to
      diverge from ours exactly by the masked-wedge flux; quantifies the fill choice.
  P2  RadialProfile  — photutils' own azimuthal mean-per-annulus. Integrating it the
      same way we integrate I(r) replicates OUR methodology with an independent
      annulus engine → isolates whether our hand-rolled annuli are correct.

Centre cross-check: centroid_2dg vs centroid_com vs centroid_quadratic.

Usage:  conda activate ISMgas; python scripts/validate_Re_photutils.py
Output: results/validate_Re_photutils.npz  (+ printed table)
"""
import os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.profiles import CurveOfGrowth, RadialProfile
from photutils.centroids import centroid_2dg, centroid_com, centroid_quadratic
from scipy.integrate import cumulative_trapezoid

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.final_sigma_e import (curve_of_growth, find_center,
                                   HST_F140W, HST_F200LP, HST_F140W_MASK,
                                   HST_F200LP_MASK, RA_DEFL, DEC_DEFL)

BANDS = {"F140W": (HST_F140W, HST_F140W_MASK),
         "F200LP": (HST_F200LP, HST_F200LP_MASK)}
SYSZ = np.load("results/photometry_systematics.npz", allow_pickle=True)
APM = np.load("results/aperture_2re_masks.npz", allow_pickle=True)
R_MAX_AS = 6.0          # match curve_of_growth default
R_STEP_AS = 0.08


def half_light_from_cumulative(r, cum):
    total = cum[-1]
    if total <= 0:
        return np.nan
    return float(np.interp(total / 2.0, cum, r))


def photutils_cog(img, cen, ps, mask):
    """P1: direct circular-aperture enclosed flux (photutils CurveOfGrowth)."""
    radii_pix = np.arange(R_STEP_AS, R_MAX_AS + R_STEP_AS, R_STEP_AS) / ps
    cog = CurveOfGrowth(np.nan_to_num(img), cen, radii_pix,
                        mask=(None if mask is None else mask))
    return half_light_from_cumulative(cog.radii * ps, cog.profile)


def photutils_radprof_integrated(img, cen, ps, mask):
    """P2: photutils azimuthal mean-per-annulus, integrated I·2πr dr (our method)."""
    edges_pix = np.arange(0, (R_MAX_AS + R_STEP_AS) / ps, R_STEP_AS / ps)
    rp = RadialProfile(np.nan_to_num(img), cen, edges_pix,
                       mask=(None if mask is None else mask))
    r = rp.radius * ps                      # annulus centres, arcsec
    I = rp.profile
    good = np.isfinite(I)
    r, I = r[good], I[good]
    if r[0] > 0:
        r = np.concatenate([[0], r]); I = np.concatenate([[I[0]], I])
    integrand = I * 2 * np.pi * r
    cum = np.concatenate([[0], cumulative_trapezoid(integrand, r)])
    return half_light_from_cumulative(r, cum)


def main():
    rows = []
    out = {}
    for bname, (ipath, mpath) in BANDS.items():
        with fits.open(ipath) as h:
            img = h[0].data.astype(float)
            w = WCS(h[0].header)
            ps = float(np.abs(proj_plane_pixel_scales(w)[0]) * 3600)
        masks = {
            "unmasked": None,
            "expert":   fits.getdata(mpath).astype(bool),
            "best":     (SYSZ[f"{bname}_global_mask"].astype(bool)
                         | APM[f"{bname}_2Re_companion"].astype(bool)),
        }
        # shared centre (centroid_2dg, same as headline pipeline) per mask
        for mlabel, mask in masks.items():
            cx, cy = find_center(img, (mask if mask is not None
                                       else np.zeros(img.shape, bool)),
                                 w, RA_DEFL, DEC_DEFL, 3.0, ps)
            ours = float(curve_of_growth(img, (cx, cy), ps, mask=mask))
            p1 = photutils_cog(img, (cx, cy), ps, mask)
            p2 = photutils_radprof_integrated(img, (cx, cy), ps, mask)
            rows.append((bname, mlabel, ours, p1, p2, p1 - ours, p2 - ours))
            out[f"{bname}_{mlabel}"] = np.array([ours, p1, p2])

        # centroid cross-check (best mask, the headline mask)
        mb = masks["best"]
        cx0, cy0 = find_center(img, mb, w, RA_DEFL, DEC_DEFL, 3.0, ps)
        half = int(np.ceil(3.0 / ps))
        y1, y2 = max(0, int(cy0) - half), min(img.shape[0], int(cy0) + half + 1)
        x1, x2 = max(0, int(cx0) - half), min(img.shape[1], int(cx0) + half + 1)
        sub = np.nan_to_num(img[y1:y2, x1:x2]); subm = mb[y1:y2, x1:x2]
        c2dg = centroid_2dg(sub, mask=subm)
        ccom = centroid_com(sub, mask=subm)
        cquad = centroid_quadratic(sub, mask=subm)
        d_com = np.hypot(*(np.array(c2dg) - np.array(ccom))) * ps
        d_quad = np.hypot(*(np.array(c2dg) - np.array(cquad))) * ps
        out[f"{bname}_centroid_delta_com_quad_arcsec"] = np.array([d_com, d_quad])
        print(f"  [{bname}] centroid_2dg vs com = {d_com*1000:.1f} mas; "
              f"vs quadratic = {d_quad*1000:.1f} mas")

    print("\n  R_e cross-validation (arcsec):")
    print(f"  {'band':7s} {'mask':9s} {'ours':>7s} {'P1 CoG':>8s} {'P2 radprof':>10s} "
          f"{'ΔP1':>7s} {'ΔP2':>7s}")
    print("  " + "-" * 56)
    for bn, ml, ours, p1, p2, d1, d2 in rows:
        print(f"  {bn:7s} {ml:9s} {ours:>7.3f} {p1:>8.3f} {p2:>10.3f} "
              f"{d1:>+7.3f} {d2:>+7.3f}")

    # headline-relevant summary: best-mask mean R_e, all 3 methods
    bm = {bn: (ours, p1, p2) for bn, ml, ours, p1, p2, *_ in rows if ml == "best"}
    mean_ours = 0.5 * (bm["F140W"][0] + bm["F200LP"][0])
    mean_p1 = 0.5 * (bm["F140W"][1] + bm["F200LP"][1])
    mean_p2 = 0.5 * (bm["F140W"][2] + bm["F200LP"][2])
    print("\n  Best-mask mean R_e:  ours={:.3f}\"  P1(CoG-drop)={:.3f}\"  "
          "P2(radprof-int)={:.3f}\"".format(mean_ours, mean_p1, mean_p2))
    print("  (headline = ours = {:.3f}\"; P2 should match ours, P1 differs by the "
          "masked-wedge flux)".format(mean_ours))
    out["best_mean_R_e"] = np.array([mean_ours, mean_p1, mean_p2])

    np.savez("results/validate_Re_photutils.npz",
             rows=np.array(rows, dtype=object), **out)
    print("\n  → results/validate_Re_photutils.npz")


if __name__ == "__main__":
    main()
