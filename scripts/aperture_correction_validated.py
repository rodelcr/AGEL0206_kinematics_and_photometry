"""Validated-fit aperture correction for the deflector total-light M* (2026-06-15).

ROOT-CAUSE FIX (audit 2026-06-15): the aperture-corrected estimators in
`aperture_matched_photometry.py` build the under-arc FILL and the BEYOND-aperture
correction from a Sersic model fit with `fit_sersic` (photometry_systematics) on the
AUTO color/morph mask, with loose bounds (r_eff up to 2.5x init=2.305", n<=8, box 6").
That fit is biased: r_eff -> 2.2-2.7" (vs the photutils-validated CoG R_e=2.097" and the
published per-band Sersic table 1.6-2.0") and ellipticity rails to circular. It puts only
70-79% of the model light inside 2 R_e (vs 84-90% for the validated fit), inflating the
beyond-aperture correction by ~0.15 mag/band and pushing log M*(total) to 11.50, ABOVE the
pure-model sersic_total (11.43).

FIX: hold the Sersic SHAPE fixed to the photutils-validated per-band table values
(n, r_eff, ellip, PA from results/sersic_parameter_table.npz) and fit ONLY amplitude+sky
to the data (auto mask). Use this validated-shape model for the fill + beyond, exactly
mirroring aperture_matched_photometry.py; keep the auto mask for the empirical F_raw.

Outputs results/aperture_correction_validated.npz with the validated total/filled/raw_apcorr
mags + Bagpipes M* posteriors + the SED model spectrum + mass-weighted age (headline-
consistent), and the apcorr-model systematic = logM(auto total) - logM(validated total).

Run from repo root:  python scripts/aperture_correction_validated.py
"""
from __future__ import annotations
import numpy as np
from astropy.modeling.models import Sersic2D

from scripts.mask_method_comparison import load_band
from scripts.aperture_2re_companions import aperture_ellipse, companion_mask, in_aperture
from scripts.bagpipes_sersic_refit import fit_and_extract

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
NRE = 2                       # headline aperture (2 R_e)
FITBOX_AS = 6.0              # amplitude+sky fit box radius (arcsec)

_TAB = np.load("results/sersic_parameter_table.npz", allow_pickle=True)
def tab(band, key):
    return float({k: v for k, v in _TAB[f"{band}_cen"]}[key])


def ab(net, zp):
    return float(-2.5 * np.log10(net) + zp) if net > 0 else float("nan")


def main():
    sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
    bands = {n: load_band(n) for n in ORDER}
    pivot = np.array([bands[n]["cfg"]["pivot_AA"] for n in ORDER])

    mags_raw, mags_filled, mags_total, mags_rawapc = [], [], [], []
    diag = {}
    for n in ORDER:
        b = bands[n]; ps = b["pix_scale"]; zp = b["ab_zp"]
        cx, cy = b["cx"], b["cy"]
        xc, yc, a, bb, th_ap = aperture_ellipse(b, NRE)
        arc = sysz[f"{n}_global_mask"].astype(bool)
        comp = companion_mask(b, arc, xc, yc, a, bb, th_ap)
        mask = arc | comp
        aper = in_aperture(b["img"].shape, xc, yc, a, bb, th_ap)

        # validated SHAPE (amplitude=1), centred at the deflector, table n/r_eff/ellip/PA
        reff_px = tab(n, "reff_as") / ps
        nser = tab(n, "n"); ell = tab(n, "_ell"); pa = tab(n, "pa")
        yy, xx = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
        shape = np.clip(np.asarray(
            Sersic2D(amplitude=1.0, r_eff=reff_px, n=nser, x_0=cx, y_0=cy,
                     ellip=ell, theta=np.radians(pa))(xx, yy), dtype=float), 0, None)

        # fit amplitude + flat sky linearly on unmasked pixels in a circular fit box
        fitbox = in_aperture(b["img"].shape, cx, cy, FITBOX_AS / ps, FITBOX_AS / ps, 0.0)
        sel = fitbox & ~mask & np.isfinite(b["img"])
        A = np.vstack([shape[sel], np.ones(int(sel.sum()))]).T
        (amp, sky), *_ = np.linalg.lstsq(A, np.nan_to_num(b["img"])[sel], rcond=None)
        model = amp * shape                          # sky-subtracted validated model
        data = np.nan_to_num(b["img"]) - sky

        ap_un = aper & ~mask
        F_raw = float(np.sum(data[ap_un]))
        filled_img = np.where(mask, model, data)
        F_filled = float(np.sum(filled_img[aper]))
        F_within = float(np.sum(model[aper])); F_tot = float(np.sum(model))
        F_beyond = max(0.0, F_tot - F_within)

        mags_raw.append(ab(F_raw, zp))
        mags_filled.append(ab(F_filled, zp))
        mags_total.append(ab(F_filled + F_beyond, zp))
        mags_rawapc.append(ab(F_raw + F_beyond, zp))
        enc = F_within / F_tot
        diag[n] = dict(reff_as=tab(n, "reff_as"), n=nser, ell=ell, amp=float(amp),
                       sky=float(sky), enclosed=enc, apcorr_mag=-2.5 * np.log10(enc),
                       m_tot_table=tab(n, "m_tot"))
        print(f"{n:7s} reff={tab(n,'reff_as'):.2f}\" n={nser:.2f} ell={ell:.2f} | "
              f"enclosed={enc:.3f} apcorr={-2.5*np.log10(enc):.3f} | "
              f"total={mags_total[-1]:.3f} (sersic-tab m_tot={tab(n,'m_tot'):.3f})")

    mags_raw = np.array(mags_raw); mags_filled = np.array(mags_filled)
    mags_total = np.array(mags_total); mags_rawapc = np.array(mags_rawapc)

    out = dict(filter_names=np.array(ORDER), pivot=pivot,
               mag_raw_2=mags_raw, mag_filled_2=mags_filled,
               mag_total_2=mags_total, mag_raw_apcorr_2=mags_rawapc)
    for n in ORDER:
        for k, vv in diag[n].items():
            out[f"diag_{n}_{k}"] = vv

    # Bagpipes M* on the validated total (HEADLINE) + filled + raw_apcorr; also pull the
    # SED model spectrum + mass-weighted age from the validated TOTAL fit (headline-consistent).
    print("\nBagpipes M* (validated)...")
    res_total = fit_and_extract(mags_total, pivot, "AGEL0206_aperVAL_total_2Re", err_frac=0.1)
    res_fill = fit_and_extract(mags_filled, pivot, "AGEL0206_aperVAL_filled_2Re", err_frac=0.1)
    res_rapc = fit_and_extract(mags_rawapc, pivot, "AGEL0206_aperVAL_rawapc_2Re", err_frac=0.1)
    for tag, res in [("total", res_total), ("filled", res_fill), ("raw_apcorr", res_rapc)]:
        s = res["stellar_mass"]; p = np.percentile(s, [16, 50, 84])
        out[f"logM_{tag}_2"] = p
        out[f"samples_logM_{tag}_2"] = s
        print(f"  {tag:11s}: logM* = {p[1]:.3f} (-{p[1]-p[0]:.3f}/+{p[2]-p[1]:.3f})")

    # headline-consistent SED model spectrum + photometry + age (from the TOTAL fit)
    out["sed_wav_obs"] = res_total["wav_obs"]
    out["sed_spec_p16"] = res_total["spec_p16"]
    out["sed_spec_p50"] = res_total["spec_p50"]
    out["sed_spec_p84"] = res_total["spec_p84"]
    out["sed_eff_wavs"] = res_total["eff_wavs"]
    out["sed_flam"] = res_total["flam"]
    out["sed_flam_err"] = res_total["flam_err"]
    out["sed_model_phot_p16"] = res_total["model_phot_p16"]
    out["sed_model_phot_p50"] = res_total["model_phot_p50"]
    out["sed_model_phot_p84"] = res_total["model_phot_p84"]
    out["sed_z_p50"] = res_total["z_p50"]
    age = res_total["mass_weighted_age"]
    if age is not None:
        a16, a50, a84 = np.percentile(np.asarray(age, float), [16, 50, 84])
        out["mass_weighted_age"] = np.asarray(age, float)
        print(f"  mass-weighted age (validated total): {a50:.2f} +{a84-a50:.2f}/-{a50-a16:.2f} Gyr")

    # apcorr-model systematic: auto-fit total (old central) vs validated total
    am = np.load("results/aperture_matched_photometry.npz", allow_pickle=True)
    logM_auto_total = float(np.median(am["samples_logM_total_2"]))
    logM_val_total = float(np.median(out["samples_logM_total_2"]))
    out["sys_apcorr_model_dex"] = abs(logM_auto_total - logM_val_total)
    print(f"\nauto total logM*={logM_auto_total:.3f}  validated total logM*={logM_val_total:.3f}")
    print(f"apcorr-model systematic = {out['sys_apcorr_model_dex']:.3f} dex")

    np.savez("results/aperture_correction_validated.npz", **out)
    print("\nSaved -> results/aperture_correction_validated.npz")


if __name__ == "__main__":
    main()
