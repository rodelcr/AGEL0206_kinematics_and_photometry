"""Matched 2 R_e aperture + companion masking for the total-light photometry.

Rationale (2026-06-10, user):
  - The per-band hand-drawn apertures are inconsistent (F200LP 1.19" ~0.5 R_e vs IR
    2.66" ~1 R_e) and capture only ~17-48% of the light -> biased, mismatched colors.
  - Fix: ONE elliptical aperture at 2 R_e (a=4.61"), same physical region in every band
    (deflector axis-ratio + PA, HST-mean centre), capturing ~80% of the light; a Sersic
    correction (+~0.25 mag) extends to total. S/N stays high to 2 R_e (>100 all bands).
  - A 2 R_e aperture ENCLOSES field companions the old tight aperture excluded -> they
    must be masked. We default to the GLOBAL arc mask (user preference) and add a
    segmentation-based companion mask: detect field sources, exclude the deflector core
    and arc, mask the rest. (2.5 R_e was rejected: triples companions for only 0.036 dex.)

This module builds the masks + aperture and renders a verification figure. Photometry
+ M* re-run are wired separately once the masks are confirmed.

Usage: conda activate ISMgas; python scripts/aperture_2re_companions.py
Output: results/figures/aperture_2re_companions.png + results/aperture_2re_masks.npz
"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.visualization import AsinhStretch, ManualInterval
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.ndimage import binary_dilation

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
# Headline R_e from the BEST-MASK (single-Sérsic + color/morph gate + WCS reg)
# F140W+F200LP curve-of-growth = 2.097" (was the expert-mask 2.305"). The 2-R_e
# companion mask is EMPTY in HST (companions live only in deep JWST), so it does
# not enter the HST R_e — global color/morph mask alone gives 2.097" (self-consistent
# fixed point). Adopted 2026-06-11 ("best mask throughout"); photutils-validated to
# ±0.002" (scripts/validate_Re_photutils.py). See scripts/re_mask_sensitivity.py.
RE = 2.097
AXIS_RATIO = 0.75              # deflector b/a (from light); applied consistently
CORE_EXCL_NRE = 1.25           # don't call anything inside this a deflector/arc "companion"
COMP_RGROW = 0.7               # arcsec: fixed companion-mask circle radius


def aperture_ellipse(b, nre):
    """Matched elliptical aperture (pixels) at nre*R_e, HST-mean deflector centre,
    deflector PA. Returns (xc, yc, a_pix, b_pix, theta_rad)."""
    p = json.load(open(b["cfg"]["params"]))
    xc, yc = b["wcs"].world_to_pixel_values(p["target_ra"], p["target_dec"])
    theta = p["pixel_geometry"]["theta_rad"]
    a_pix = nre * RE / b["pix_scale"]
    return float(xc), float(yc), a_pix, a_pix * AXIS_RATIO, theta


def in_aperture(shape, xc, yc, a, bb, theta, scale=1.0):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    ct, st = np.cos(-theta), np.sin(-theta)
    xr = (xx - xc) * ct - (yy - yc) * st
    yr = (xx - xc) * st + (yy - yc) * ct
    return (xr / (a * scale))**2 + (yr / (bb * scale))**2 <= 1.0


def companion_mask(b, arc_mask, xc, yc, a, bb, theta):
    """Field-companion mask: DAOStarFinder peaks masked with FIXED circles.

    (deblend_sources needs skimage, absent in ISMgas; flux-grow leaks into the arc.)
    A peak is masked as a field companion iff it is (a) beyond CORE_EXCL_NRE R_e from
    the deflector (neither deflector nor inner arc), (b) within 1.25x the aperture, (c)
    not already an arc pixel, and (d) reasonably ROUND (|roundness|<0.7) so elongated
    arc substructure is rejected. Each kept peak masks a FIXED COMP_RGROW-radius circle
    (no flux-grow → cannot leak into the extended arc/deflector). Beyond 1.25 R_e the
    deflector is faint, so 6sigma peaks there are real compact sources."""
    img = np.nan_to_num(b["img"]); ps = b["pix_scale"]
    mean, med, std = sigma_clipped_stats(img, sigma=3.0)
    t = DAOStarFinder(fwhm=max(3.0, 0.25 / ps), threshold=6 * std, roundlo=-0.7, roundhi=0.7)(img - med)
    out = np.zeros(img.shape, bool)
    if t is None:
        return out
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    reach = in_aperture(img.shape, xc, yc, a, bb, theta, scale=1.25)
    kept = []
    for r in t:
        px, py = float(r["xcentroid"]), float(r["ycentroid"])
        ipy, ipx = int(round(py)), int(round(px))
        if not (0 <= ipy < img.shape[0] and 0 <= ipx < img.shape[1]):
            continue
        if np.hypot(px - xc, py - yc) * ps < CORE_EXCL_NRE * RE:     # deflector / inner arc
            continue
        if not reach[ipy, ipx] or arc_mask[ipy, ipx]:               # outside aper / arc pixel
            continue
        out |= np.hypot(xx - px, yy - py) * ps < COMP_RGROW          # fixed circle
        kept.append((px, py))
    return binary_dilation(out, iterations=1)


def main():
    bands = {n: load_band(n) for n in ORDER}
    sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
    # rows: 2 R_e (headline) and 2.5 R_e (aperture cross-check)
    fig, axes = plt.subplots(2, 4, figsize=(20, 11))
    save = {}
    for row, nre in enumerate((2.0, 2.5)):
        for ax, n in zip(axes[row], ORDER):
            b = bands[n]; ps = b["pix_scale"]
            xc, yc, a, bb, theta = aperture_ellipse(b, nre)
            arc = sysz[f"{n}_global_mask"].astype(bool)     # global mask (user default)
            comp = companion_mask(b, arc, xc, yc, a, bb, theta)
            tag = f"{n}_{nre:g}Re"
            save[f"{tag}_companion"] = comp
            save[f"{tag}_full"] = arc | comp
            save[f"{tag}_aperture"] = np.array([xc, yc, a, bb, theta])
            half = int(np.ceil((nre + 0.6) * RE / ps)); iy, ix = int(yc), int(xc)
            sl = (slice(max(0, iy - half), iy + half), slice(max(0, ix - half), ix + half))
            y0, x0 = sl[0].start, sl[1].start
            img = np.nan_to_num(b["img"])[sl]
            lo, hi = np.nanpercentile(img, 35), np.nanpercentile(img, 99.3)
            ax.imshow(AsinhStretch(0.05)(ManualInterval(lo, hi)(img)), origin="lower", cmap="gray")
            ax.contourf(arc[sl].astype(float), levels=[.5, 1.5], colors="lime", alpha=0.28)
            ax.contour(arc[sl].astype(float), levels=[.5], colors="lime", linewidths=0.9)
            if comp.any():
                ax.contourf(comp[sl].astype(float), levels=[.5, 1.5], colors="red", alpha=0.45)
                ax.contour(comp[sl].astype(float), levels=[.5], colors="red", linewidths=1.0)
            ax.add_patch(Ellipse((xc - x0, yc - y0), 2 * a, 2 * bb, angle=np.degrees(theta),
                                 fill=False, ec="yellow", lw=1.6, ls="--"))
            ax.set_title(f"{n} @ {nre:g} R_e: lime=arc  red=companions({int(comp.sum())}px)",
                         fontsize=9)
            ax.set_xlim(0, img.shape[1]); ax.set_ylim(0, img.shape[0])
            ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Matched aperture (yellow) + global arc mask (lime) + companion mask (red) — "
                 "top: 2 R_e (headline), bottom: 2.5 R_e (cross-check)", fontsize=13, y=1.005)
    plt.tight_layout()
    plt.savefig("results/figures/aperture_2re_companions.png", dpi=150, bbox_inches="tight")
    np.savez("results/aperture_2re_masks.npz", **save)
    print("Saved → results/figures/aperture_2re_companions.png + results/aperture_2re_masks.npz")
    for nre in (2.0, 2.5):
        for n in ORDER:
            print(f"  {n} @ {nre:g}Re: companions={int(save[f'{n}_{nre:g}Re_companion'].sum())}px")


if __name__ == "__main__":
    main()
