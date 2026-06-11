"""Matched-aperture total-light photometry + M* at 1.0 / 2.0 / 2.5 R_e.

For each aperture radius nre*R_e (consistent elliptical aperture, deflector PA/axis-
ratio, HST-mean centre, ALL bands the same physical region) and each band:
  - global arc mask (user default) UNION the companion mask (scripts.aperture_2re_companions)
  - F_raw    = sum(data - sky) over [aperture & ~mask]          (empirical aperture flux)
  - F_filled = same but masked pixels replaced by the single-Sersic model (recovers the
               deflector light hidden under the arc/companions WITHIN the aperture)
  - aperture correction to total: +(-2.5 log10 frac), frac = Sersic-enclosed fraction
    within nre R_e from the fitted n  ->  mag_total = mag_filled + apcorr
Three AB-mag vectors per radius: aperture-raw, filled (within-aperture), total.
Bagpipes (nb02 priors, 10% errors) on aperture-raw and total -> M*.

Compares 1 R_e vs 2 R_e: the aperture-raw M* differs (captures less light), the TOTAL
should converge -> validates the aperture correction. 2.5 R_e is the upper cross-check.

Usage: conda activate ISMgas; python scripts/aperture_matched_photometry.py
Output: results/aperture_matched_photometry.npz
"""
import os, sys, json
import numpy as np
from scipy.special import gammainc
from scipy.optimize import brentq

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band
from scripts.arc_mask_verification import sky_sigma
from scripts.photometry_systematics import fit_sersic, reg_shift_for
from scripts.mask_method_comparison import reproject_mask
from scripts.aperture_2re_companions import aperture_ellipse, companion_mask, in_aperture, RE
from scripts.bagpipes_sersic_refit import run_bagpipes_for_mags

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
NRES = [1.0, 2.0, 2.5]


def bn(n):
    return brentq(lambda b: gammainc(2 * n, b) - 0.5, 0.1, 30)


def enclosed_frac(n, nre):
    """Fraction of TOTAL Sersic flux within nre*R_e (R_e = half-light)."""
    return float(gammainc(2 * n, bn(n) * nre ** (1.0 / n)))


def ab(net, zp):
    return float(-2.5 * np.log10(net) + zp) if net > 0 else float("nan")


def main():
    bands = {n: load_band(n) for n in ORDER}
    f200 = bands["F200LP"]; src200 = (f200["mask"].astype(bool), f200["wcs"])
    sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
    ser = np.load("results/sersic_total_photometry.npz", allow_pickle=True)
    sfn = list(ser["filter_names"])
    pivot = np.array([bands[n]["cfg"]["pivot_AA"] for n in ORDER])

    # per-band single-Sersic model (sky-sub) + sky, fit once (registered F200 seed mask)
    models, skys = {}, {}
    for n in ORDER:
        b = bands[n]
        M = reproject_mask(*src200, b["wcs"], b["img"].shape)
        sh = reg_shift_for(f200, b)
        if sh != (0, 0):
            from scipy.ndimage import shift as ndshift
            M = ndshift(M.astype(float), sh, order=0, mode="constant") > 0.5
        model, sky_med = fit_sersic(b, M)
        models[n] = model; skys[n] = sky_med

    out = {"filter_names": np.array(ORDER), "pivot": pivot, "nres": np.array(NRES)}
    for nre in NRES:
        mags_raw, mags_filled, mags_total, mags_raw_apcorr = [], [], [], []
        for n in ORDER:
            b = bands[n]; ps = b["pix_scale"]; zp = b["ab_zp"]
            xc, yc, a, bb, theta = aperture_ellipse(b, nre)
            arc = sysz[f"{n}_global_mask"].astype(bool)
            comp = companion_mask(b, arc, xc, yc, a, bb, theta)
            mask = arc | comp
            aper = in_aperture(b["img"].shape, xc, yc, a, bb, theta)
            data = np.nan_to_num(b["img"]) - skys[n]
            model = models[n]
            sky = skys[n]
            # local sky re-estimate in an arc-free ring (robust)
            _, _ = sky_sigma(b["img"], xc, yc, ps)
            ap_un = aper & ~mask
            F_raw = float(np.sum(data[ap_un]))
            filled_img = np.where(mask, model, data)
            F_filled = float(np.sum(filled_img[aper]))
            # aperture correction to total: ADD the SAME single-Sersic model's light from
            # OUTSIDE the aperture (additive, model-consistent; handles the elliptical
            # aperture + profile shape exactly, no R_e/r_eff convention issue). model is
            # the sky-subtracted Sersic rendered over the full cutout (~ total to >10 R_e).
            F_model_within = float(np.sum(model[aper]))
            F_model_total = float(np.sum(model))
            F_beyond = max(0.0, F_model_total - F_model_within)
            mags_raw.append(ab(F_raw, zp))
            mags_filled.append(ab(F_filled, zp))
            mags_total.append(ab(F_filled + F_beyond, zp))
            # MOST-EMPIRICAL total: raw aperture (masked pixels NOT filled) + only the
            # model's beyond-aperture wings. Models only the unavoidable outer light;
            # the masked interior is accepted as a (small) loss, not model-filled.
            mags_raw_apcorr.append(ab(F_raw + F_beyond, zp))
        out[f"mag_raw_{nre:g}"] = np.array(mags_raw)
        out[f"mag_filled_{nre:g}"] = np.array(mags_filled)
        out[f"mag_total_{nre:g}"] = np.array(mags_total)
        out[f"mag_raw_apcorr_{nre:g}"] = np.array(mags_raw_apcorr)
        print(f"\n=== {nre:g} R_e ===")
        for i, n in enumerate(ORDER):
            print(f"  {n}: raw {mags_raw[i]:.3f}  filled {mags_filled[i]:.3f}  total {mags_total[i]:.3f}")

    # Bagpipes M* for all four estimators, at each radius:
    #   raw    = empirical aperture (contaminants masked, not filled)
    #   filled = + masked pixels filled with the single-Sersic model (within aperture)
    #   total  = + the Sersic model's light beyond the aperture (aperture-corrected total)
    #   (Sersic-total, the pure-model estimator, is in results/bagpipes_sersic_refit.npz)
    print("\nBagpipes M* ...")
    for nre in NRES:
        for kind in ("raw", "raw_apcorr", "filled", "total"):
            mags = out[f"mag_{kind}_{nre:g}"]
            if not np.all(np.isfinite(mags)):
                print(f"  {nre:g}Re {kind}: non-finite mag, skip")
                continue
            run = f"AGEL0206_aper_{kind}_{nre:g}Re"
            samp = run_bagpipes_for_mags(mags, pivot, run, err_frac=0.1)
            p = np.percentile(samp, [16, 50, 84])
            out[f"logM_{kind}_{nre:g}"] = p
            out[f"samples_logM_{kind}_{nre:g}"] = samp
            print(f"  {nre:g}Re {kind}: logM* = {p[1]:.3f} (-{p[1]-p[0]:.3f}/+{p[2]-p[1]:.3f})")

    np.savez("results/aperture_matched_photometry.npz", **out)
    print("\nSaved → results/aperture_matched_photometry.npz")


if __name__ == "__main__":
    main()
