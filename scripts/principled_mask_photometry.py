"""Principled (Sersic-residual k=3) arc masking applied to ALL FOUR photometry
bands — F200LP, F140W (HST) + F150W2, F322W2 (JWST) — as an audit + the basis for
adopting it as the headline photometry standard.

Two purposes:
  1. AUDIT: verify our aperture photometry reproduces the official headline AB mags
     under the EXPERT mask (sanity that the Δmag deltas below are trustworthy).
  2. EXTEND: apply the Sersic-residual k=3 mask (the standard adopted from
     scripts/arc_mask_verification.py) to every band, including JWST, so we
     understand the change in each before touching the headline M★/R_e.

Reuses validated building blocks:
  - scripts.sersic_total_photometry.load_hst / load_jwst / find_center_2dg
  - scripts.arc_mask_verification._fit_sersic_local / sky_sigma /
    sersic_residual_mask / core_exclude

Usage
-----
    conda activate ISMgas
    python scripts/principled_mask_photometry.py

Output: results/principled_mask_photometry.npz + results/figures/principled_{band}.png
"""
import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from photutils.aperture import EllipticalAperture, EllipticalAnnulus, ApertureStats

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.sersic_total_photometry import load_hst, load_jwst, find_center_2dg
from scripts.arc_mask_verification import (_fit_sersic_local, sky_sigma,
                                           sersic_residual_mask, core_exclude)

VDI = "../velocity_dispersion_from_IFU"
RA_DEFL, DEC_DEFL = 31.55611, -1.23817
R_E_INIT = 2.305
FIG_DIR = "results/figures"
K_SIGMA = 3.0          # adopted Sersic-residual threshold (see arc_mask_verification)
SNR_FLOOR = 3.0
CORE_R = 0.4

# Headline aperture AB mags (nb02) — the audit target
BANDS = {
    "F200LP": dict(kind="hst", image=f"{VDI}/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits",
                   mask=f"{VDI}/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits",
                   params="example_outputs/AGEL020613-011417A_F200LP_WFC3_cutout_L3_params.json",
                   photflam=5.1851e-20, photplam=4923.48, pivot_AA=4972.0, ap_mag=22.6126),
    "F140W": dict(kind="hst", image=f"{VDI}/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits",
                  mask=f"{VDI}/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits",
                  params="example_outputs/AGEL020613-011417A_F140W_WFC3_cutout_L3_params.json",
                  photflam=1.4737e-20, photplam=13922.9, pivot_AA=13923.0, ap_mag=19.1335),
    "F150W2": dict(kind="jwst", image=f"{VDI}/jw05594-o101_t103_nircam_clear-f150w2_i2d.fits",
                   mask=f"{VDI}/jw05594-o101_t103_nircam_clear-f150w2_i2d_mask.fits",
                   params="example_outputs/jw05594-o101_t103_nircam_clear-f150w2_i2d_params.json",
                   pivot_AA=16720.0, ap_mag=18.9425),
    "F322W2": dict(kind="jwst", image=f"{VDI}/jw05594-o101_t103_nircam_clear-f322w2_i2d.fits",
                   mask=f"{VDI}/jw05594-o101_t103_nircam_clear-f322w2_i2d_mask.fits",
                   params="example_outputs/jw05594-o101_t103_nircam_clear-f322w2_i2d_params.json",
                   pivot_AA=32470.0, ap_mag=18.6042),
}


def load_band(name):
    cfg = BANDS[name]
    if cfg["kind"] == "hst":
        b = load_hst(cfg["image"], cfg["mask"], cfg["photflam"], cfg["photplam"])
        # ZP consistent with photometry_masking_HST (reproduces the headline mags)
        b["ab_zp"] = (-2.5 * np.log10(cfg["photflam"]) - 21.10
                      - 5 * np.log10(cfg["photplam"]) + 18.69)
    else:
        b = load_jwst(cfg["image"], cfg["mask"])
        b["ab_zp"] = b["zp_ab"]                      # -6.10 - 2.5 log10(PIXAR_SR)
    b["name"], b["cfg"] = name, cfg
    cx, cy, _ = find_center_2dg(b["img"], b["mask"], b["wcs"], RA_DEFL, DEC_DEFL,
                                box_arcsec=3.0, pix_scale=b["pix_scale"])
    b["cx"], b["cy"] = cx, cy
    return b


def aperture_mag(b, mask):
    """Elliptical-aperture AB mag (params.json geometry, band ZP). mask True=excluded."""
    p = json.load(open(b["cfg"]["params"]))
    pg = p["pixel_geometry"]
    xc, yc = b["wcs"].world_to_pixel_values(p["target_ra"], p["target_dec"])
    aper = EllipticalAperture((float(xc), float(yc)), a=pg["a"], b=pg["b"], theta=pg["theta_rad"])
    ann = EllipticalAnnulus((float(xc), float(yc)), a_in=pg["annulus_a_in"],
                            a_out=pg["annulus_a_out"],
                            b_out=pg["annulus_a_out"] * (pg["b"] / pg["a"]), theta=pg["theta_rad"])
    data = np.nan_to_num(b["img"])
    tot = None if mask is None else np.asarray(mask).astype(bool)
    aps = ApertureStats(data, aper, mask=tot)
    ans = ApertureStats(data, ann, mask=tot)
    net = aps.sum - ans.median * aps.sum_aper_area.value
    return float(-2.5 * np.log10(net) + b["ab_zp"]) if net > 0 else float("nan")


def analyze(name):
    b = load_band(name)
    img, ps, cx, cy = b["img"], b["pix_scale"], b["cx"], b["cy"]
    sky_med, sky_sig = sky_sigma(img, cx, cy, ps)
    pg = json.load(open(b["cfg"]["params"]))
    ellip_init = float(pg["geometry"]["ellipticity"])
    theta_init = float(pg["pixel_geometry"]["theta_rad"])
    box = 5.0
    fit, (x1, y1) = _fit_sersic_local(img, b["mask"], (cx, cy), R_E_INIT, box, ps,
                                      ellip_init, theta_init)
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    model = np.clip(np.asarray(fit(xx - x1, yy - y1), dtype=float), 0, None)
    resid = np.nan_to_num(img) - sky_med - model
    flux = np.nan_to_num(img) - sky_med
    excl = core_exclude(img.shape, cx, cy, ps, CORE_R)
    m_sersic = sersic_residual_mask(resid, sky_sig, K_SIGMA, flux, SNR_FLOOR, excl)
    expert = b["mask"].astype(bool)

    mags = {lbl: aperture_mag(b, m) for lbl, m in
            (("none", None), ("expert", expert), ("sersic", m_sersic))}
    # over-masking: fraction of masked flux from the smooth model
    mm = m_sersic
    md = np.clip(flux[mm], 0, None).sum()
    mf = float(np.clip(model[mm], 0, None).sum() / md) if md > 0 else float("nan")
    audit = mags["expert"] - b["cfg"]["ap_mag"]   # should be ~0 if photometry reproduces headline

    print(f"\n{name} ({b['cfg']['kind']}, {ps:.4f}\"/pix): "
          f"Sersic n={fit.n.value:.2f} reff={fit.r_eff.value*ps:.2f}\" ellip={fit.ellip.value:.2f}")
    print(f"  AUDIT expert vs headline ap_mag: {mags['expert']:.4f} vs {b['cfg']['ap_mag']:.4f} "
          f"(Δ {audit:+.4f})")
    print(f"  none={mags['none']:.4f}  expert={mags['expert']:.4f}  "
          f"sersic-k3={mags['sersic']:.4f}  (Δ sersic−expert {mags['sersic']-mags['expert']:+.4f})")
    print(f"  sersic mask: n_px={int(mm.sum())}  model_flux_frac={mf:.2f}")

    # audit figure: data / model / residual(S/N) + mask contour, zoomed
    half = int(np.ceil((box + 1) / ps))
    sl = (slice(max(0, int(cy - half)), int(cy + half)),
          slice(max(0, int(cx - half)), int(cx + half)))
    data_ss = flux
    vmax = float(np.nanpercentile(np.abs(data_ss[sl]), 99))
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.6))
    for a, arr, t, cm, vr in (
        (ax[0], data_ss, "data (sky-sub)", "magma", (0, vmax)),
        (ax[1], model, "Sersic model", "magma", (0, vmax)),
        (ax[2], resid, "residual", "RdBu_r", (-6 * sky_sig, 6 * sky_sig)),
        (ax[3], resid / sky_sig, "residual/σ_sky +mask", "RdBu_r", (-6, 6))):
        im = a.imshow(arr[sl], origin="lower", cmap=cm, vmin=vr[0], vmax=vr[1])
        a.set_title(f"{name} — {t}", fontsize=10); a.set_xticks([]); a.set_yticks([])
        plt.colorbar(im, ax=a, fraction=0.046)
    ax[3].contour(m_sersic[sl].astype(float), levels=[0.5], colors="lime", linewidths=0.7)
    plt.tight_layout()
    out = f"{FIG_DIR}/principled_{name}.png"
    plt.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print(f"  saved {out}")

    return dict(name=name, kind=b["cfg"]["kind"], pivot=b["cfg"]["pivot_AA"],
                ap_headline=b["cfg"]["ap_mag"], audit_delta=audit,
                mag_none=mags["none"], mag_expert=mags["expert"], mag_sersic=mags["sersic"],
                d_sersic_expert=mags["sersic"] - mags["expert"],
                n_sersic_px=int(mm.sum()), model_flux_frac=mf,
                sersic_n=float(fit.n.value), sersic_ellip=float(fit.ellip.value))


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    order = ["F200LP", "F140W", "F150W2", "F322W2"]
    res = [analyze(n) for n in order]

    print("\n" + "=" * 78)
    print("4-BAND SUMMARY  (principled = Sersic-residual k=3)")
    print("=" * 78)
    print(f"{'band':<8}{'pivot':>7}{'headline':>9}{'audit Δ':>9}{'none':>9}"
          f"{'expert':>9}{'sersic':>9}{'ser−exp':>9}{'mdlfrac':>8}")
    for r in res:
        print(f"{r['name']:<8}{r['pivot']:>7.0f}{r['ap_headline']:>9.4f}{r['audit_delta']:>+9.4f}"
              f"{r['mag_none']:>9.4f}{r['mag_expert']:>9.4f}{r['mag_sersic']:>9.4f}"
              f"{r['d_sersic_expert']:>+9.4f}{r['model_flux_frac']:>8.2f}")

    save = {r["name"]: np.array([r], dtype=object) for r in res}
    save["mag_headline"] = np.array([r["ap_headline"] for r in res])
    save["mag_expert"] = np.array([r["mag_expert"] for r in res])
    save["mag_sersic"] = np.array([r["mag_sersic"] for r in res])
    save["pivot"] = np.array([r["pivot"] for r in res])
    save["filter_names"] = np.array(order)
    np.savez("results/principled_mask_photometry.npz", **save)
    print("\nSaved → results/principled_mask_photometry.npz")


if __name__ == "__main__":
    main()
