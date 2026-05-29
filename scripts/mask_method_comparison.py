"""Mask-method comparison + Sersic fill-in systematic, all four photometry bands.

For each band (F200LP, F140W, F150W2, F322W2) and each of three arc masks:
    - expert      : the band's own hand-painted mask
    - hst_reproj  : union of F200LP & F140W expert masks reprojected onto this grid
    - sersic      : independent per-band Sersic-residual mask (k=3)
re-fit a NEW per-band Sersic2D with those masked pixels EXCLUDED ("photometry
re-informs the Sersic"), then report three flux estimates per (band, mask):
    raw    : photutils aperture sum with masked pixels excluded  (UNDER-counts the
             deflector light hidden under the arc)
    filled : replace ONLY the masked pixels inside the aperture with the re-fit
             Sersic model, then sum  (nb08 'model-filled' — recovers under-arc light)
    total  : analytic integrate-to-infinity Sersic total (aperture-independent)

Motivation: the independent per-band Sersic over-masks on F140W/JWST because a single
Sersic under-fits the bright central galaxy (see scripts/principled_mask_photometry.py).
Reprojecting the HST-defined arc footprint sidesteps that; the fill-in quantifies how
much deflector light the mask removes from under the arc.

Usage:  conda activate ISMgas; python scripts/mask_method_comparison.py
Output: results/mask_method_comparison.npz + results/figures/maskcompare_{band}.png
"""
import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from photutils.aperture import EllipticalAperture, EllipticalAnnulus, ApertureStats

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.sersic_total_photometry import (load_hst, load_jwst, find_center_2dg,
                                             sersic_total_flux_analytic)
from scripts.arc_mask_verification import (_fit_sersic_local, sky_sigma,
                                           sersic_residual_mask, core_exclude)
from scripts.principled_mask_photometry import BANDS, RA_DEFL, DEC_DEFL, R_E_INIT

FIG_DIR = "results/figures"
K_SIGMA, SNR_FLOOR, CORE_R = 3.0, 3.0, 0.4
ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]


def load_band(name):
    cfg = BANDS[name]
    if cfg["kind"] == "hst":
        b = load_hst(cfg["image"], cfg["mask"], cfg["photflam"], cfg["photplam"])
        b["ab_zp"] = (-2.5 * np.log10(cfg["photflam"]) - 21.10
                      - 5 * np.log10(cfg["photplam"]) + 18.69)
    else:
        b = load_jwst(cfg["image"], cfg["mask"])
        b["ab_zp"] = b["zp_ab"]
    b["name"], b["cfg"] = name, cfg
    cx, cy, _ = find_center_2dg(b["img"], b["mask"], b["wcs"], RA_DEFL, DEC_DEFL,
                                box_arcsec=3.0, pix_scale=b["pix_scale"])
    b["cx"], b["cy"] = cx, cy
    return b


def reproject_mask(src_mask, src_wcs, dst_wcs, dst_shape):
    """Nearest-neighbour reprojection of a boolean mask onto dst grid via WCS."""
    yy, xx = np.mgrid[:dst_shape[0], :dst_shape[1]]
    ra, dec = dst_wcs.pixel_to_world_values(xx.ravel(), yy.ravel())
    xs, ys = src_wcs.world_to_pixel_values(ra, dec)
    out = map_coordinates(src_mask.astype(float), [ys, xs], order=0,
                          mode="constant", cval=0.0)
    return out.reshape(dst_shape) > 0.5


def _aperture(b, p):
    pg = p["pixel_geometry"]
    xc, yc = b["wcs"].world_to_pixel_values(p["target_ra"], p["target_dec"])
    aper = EllipticalAperture((float(xc), float(yc)), a=pg["a"], b=pg["b"], theta=pg["theta_rad"])
    ann = EllipticalAnnulus((float(xc), float(yc)), a_in=pg["annulus_a_in"],
                            a_out=pg["annulus_a_out"],
                            b_out=pg["annulus_a_out"] * (pg["b"] / pg["a"]), theta=pg["theta_rad"])
    return aper, ann


def flux_mag(b, data, aper, ann, mask):
    aps = ApertureStats(np.nan_to_num(data), aper, mask=mask)
    ans = ApertureStats(np.nan_to_num(data), ann, mask=mask)
    net = aps.sum - ans.median * aps.sum_aper_area.value
    return float(-2.5 * np.log10(net) + b["ab_zp"]) if net > 0 else float("nan")


def refit_render(b, mask, ellip_init, theta_init):
    """NEW per-band Sersic excluding `mask`; return (model_full sky-sub, fit)."""
    fit, (x1, y1) = _fit_sersic_local(b["img"], mask, (b["cx"], b["cy"]), R_E_INIT, 5.0,
                                      b["pix_scale"], ellip_init, theta_init)
    yy, xx = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    model = np.clip(np.asarray(fit(xx - x1, yy - y1), dtype=float), 0, None)
    return model, fit


def analyze(name, hst_master):
    b = load_band(name)
    ps, cx, cy = b["pix_scale"], b["cx"], b["cy"]
    sky_med, sky_sig = sky_sigma(b["img"], cx, cy, ps)
    p = json.load(open(b["cfg"]["params"]))
    ellip_init = float(p["geometry"]["ellipticity"])
    theta_init = float(p["pixel_geometry"]["theta_rad"])
    aper, ann = _aperture(b, p)
    expert = b["mask"].astype(bool)

    # HST-reproj mask = union of F200LP & F140W expert masks reprojected here
    hst_reproj = np.zeros(b["img"].shape, bool)
    for (m_src, w_src) in hst_master:
        hst_reproj |= reproject_mask(m_src, w_src, b["wcs"], b["img"].shape)

    # independent per-band Sersic-residual mask (first-pass fit, arc down-weighted by expert)
    model0, _ = refit_render(b, expert, ellip_init, theta_init)
    resid0 = np.nan_to_num(b["img"]) - sky_med - model0
    flux = np.nan_to_num(b["img"]) - sky_med
    excl = core_exclude(b["img"].shape, cx, cy, ps, CORE_R)
    sersic_mask = sersic_residual_mask(resid0, sky_sig, K_SIGMA, flux, SNR_FLOOR, excl)

    masks = {"expert": expert, "hst_reproj": hst_reproj, "sersic": sersic_mask}
    rows = {}
    models = {}
    for mlabel, M in masks.items():
        model, fit = refit_render(b, M, ellip_init, theta_init)   # photometry re-informs Sersic
        models[mlabel] = (model, fit)
        raw = flux_mag(b, b["img"], aper, ann, M)
        data_filled = np.where(M, model + sky_med, np.nan_to_num(b["img"]))
        filled = flux_mag(b, data_filled, aper, ann, None)
        F_tot = sersic_total_flux_analytic(fit.amplitude.value, fit.r_eff.value,
                                           fit.n.value, fit.ellip.value)
        total = float(-2.5 * np.log10(F_tot) + b["ab_zp"]) if F_tot > 0 else float("nan")
        rows[mlabel] = dict(n_px=int(np.asarray(M).sum()), raw=raw, filled=filled,
                            total=total, fill_corr=filled - raw,
                            sersic_n=float(fit.n.value), sersic_reff=float(fit.r_eff.value * ps),
                            sersic_ellip=float(fit.ellip.value))
        print(f"  {name:<7} {mlabel:<11} n_px={rows[mlabel]['n_px']:>7d} "
              f"raw={raw:.4f} filled={filled:.4f} (fill {filled-raw:+.4f}) total={total:.4f} "
              f"[n={fit.n.value:.2f} re={fit.r_eff.value*ps:.2f}\" e={fit.ellip.value:.2f}]")

    # figure: data + the three masks as contours, and the expert-refit filled image
    half = int(np.ceil(6.0 / ps))
    sl = (slice(max(0, int(cy - half)), int(cy + half)),
          slice(max(0, int(cx - half)), int(cx + half)))
    data_ss = flux
    vmax = float(np.nanpercentile(np.abs(data_ss[sl]), 99))
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    im = ax[0].imshow(data_ss[sl], origin="lower", cmap="magma", vmin=0, vmax=vmax)
    for M, c, lbl in ((expert, "cyan", "expert"), (hst_reproj, "lime", "HST-reproj"),
                      (sersic_mask, "red", "indep-Sersic")):
        ax[0].contour(M[sl].astype(float), levels=[0.5], colors=c, linewidths=0.8)
    ax[0].set_title(f"{name} data + masks (cyan=expert, lime=HSTreproj, red=Sersic)", fontsize=9)
    em, _ = models["expert"]
    filled_img = np.where(expert, em + sky_med, np.nan_to_num(b["img"]))
    ax[1].imshow((filled_img - sky_med)[sl], origin="lower", cmap="magma", vmin=0, vmax=vmax)
    ax[1].set_title(f"{name} expert-mask Sersic fill-in", fontsize=9)
    ax[2].imshow((em)[sl], origin="lower", cmap="magma", vmin=0, vmax=vmax)
    ax[2].set_title(f"{name} re-fit Sersic model", fontsize=9)
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    plt.colorbar(im, ax=ax[0], fraction=0.046)
    plt.tight_layout()
    out = f"{FIG_DIR}/maskcompare_{name}.png"
    plt.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print(f"  saved {out}")
    return dict(name=name, pivot=b["cfg"]["pivot_AA"], ap_headline=b["cfg"]["ap_mag"],
                rows=rows)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    # HST master masks for reprojection (F200LP detailed + F140W)
    f200 = load_band("F200LP"); f140 = load_band("F140W")
    hst_master = [(f200["mask"].astype(bool), f200["wcs"]),
                  (f140["mask"].astype(bool), f140["wcs"])]
    print("Mask-method comparison (raw photutils / Sersic fill-in / analytic total)\n")
    res = [analyze(n, hst_master) for n in ORDER]

    print("\n" + "=" * 96)
    print("SUMMARY  — AB mag per (band, mask method)")
    print("=" * 96)
    print(f"{'band':<8}{'headline':>9}  {'method':<11}{'raw':>9}{'filled':>9}{'total':>9}"
          f"{'fill_corr':>10}{'n_px':>9}")
    for r in res:
        for mlabel, d in r["rows"].items():
            print(f"{r['name']:<8}{r['ap_headline']:>9.4f}  {mlabel:<11}{d['raw']:>9.4f}"
                  f"{d['filled']:>9.4f}{d['total']:>9.4f}{d['fill_corr']:>+10.4f}{d['n_px']:>9d}")

    np.savez("results/mask_method_comparison.npz",
             results=np.array(res, dtype=object), order=np.array(ORDER))
    print("\nSaved → results/mask_method_comparison.npz")


if __name__ == "__main__":
    main()
