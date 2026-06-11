"""Validate the aperture-corrected photometry machinery against established codes
(2026-06-11). petrofit/statmorph are absent from ISMgas, so the references are:

  * photutils.aperture.EllipticalAperture  — exact sub-pixel aperture photometry
    (industry standard) vs our hard-edge pixel-centre `in_aperture` boolean sum.
  * scipy.special incomplete-gamma — the EXACT analytic Sérsic enclosed-light law
    L(<R)/L_tot = P(2n, b_n (R/R_e)^{1/n})  and  b_n = gammaincinv(2n, 1/2)
    (Graham & Driver 2005) — vs our Ciotti&Bertin b_n approximation, our analytic
    total-flux formula, and the rendered-model aperture correction.

Four checks per band (HST F200LP/F140W, JWST F150W2/F322W2), best (global) mask:

  1. b_n approximation:  Ciotti&Bertin 1999  vs  exact gammaincinv(2n, 0.5).
  2. Sérsic total flux:  our `sersic_total_flux_analytic`  vs  (a) independent
     closed-form with EXACT b_n + scipy gamma,  (b) numeric render-to-20 R_e.
     Also: production `sum(model_full)` (rendered over the finite cutout) vs the
     to-infinity total → the cutout-truncation loss in the apcorr denominator.
  3. Enclosed-light fraction (circular):  rendered sum(model[r<2R_e])/total  vs
     analytic gammainc(2n, b_n 2^{1/n}) → validates the aperture correction.
  4. Empirical aperture flux:  our `in_aperture` boolean sum  vs  photutils
     EllipticalAperture(method='exact') and ('center'), on the real band image
     (validates F_raw) and on the rendered model (validates F_within).

Output: results/validate_apcorr_established.npz (+ printed table)
Usage:  conda activate ISMgas; python scripts/validate_apcorr_established.py
"""
import os, sys
import numpy as np
from scipy.special import gamma as Gamma, gammainc, gammaincinv
from astropy.modeling.models import Sersic2D
from photutils.aperture import EllipticalAperture, aperture_photometry

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band
from scripts.sersic_total_photometry import (fit_sersic2d, sersic_total_flux_analytic,
                                             b_n as b_n_approx, R_E_INIT)
from scripts.aperture_2re_companions import aperture_ellipse, in_aperture, RE, AXIS_RATIO

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
SYSZ = np.load("results/photometry_systematics.npz", allow_pickle=True)
NRE = 2.0


def render_full(fit, shape, x1, y1):
    """Render the subimage-frame fit back onto the full image (production framing:
    fit_sersic shifts the full grid by the subimage origin x1,y1 before evaluating)."""
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    return np.clip(np.asarray(fit(xx - x1, yy - y1), dtype=float), 0, None)


def subimage_origin(cx, cy, shape, box_arcsec, ps):
    """The (x1,y1) fit_sersic2d uses for its subimage box (must match exactly)."""
    half = int(np.ceil(box_arcsec / ps))
    y1 = max(0, int(cy) - half); x1 = max(0, int(cx) - half)
    return x1, y1


def main():
    rows = {}
    print(f"{'band':7s} {'n':>5s} {'bn_appx':>8s} {'bn_exact':>8s} {'Δbn%':>6s} "
          f"{'Ftot_an/num':>11s} {'cutoutloss%':>11s} {'encl2Re_mod':>11s} "
          f"{'encl2Re_an':>10s} {'ap_exact/ctr':>12s} {'ap/ours':>8s}")
    print("-" * 110)
    for bn_ in ORDER:
        b = load_band(bn_)
        mask = SYSZ[f"{bn_}_global_mask"].astype(bool)
        ps = b["pix_scale"]
        fit, _, _ = fit_sersic2d(b["img"], mask, (b["cx"], b["cy"]), R_E_INIT, 6.0, ps)
        amp, reff, n, ell = (float(fit.amplitude.value), float(fit.r_eff.value),
                             float(fit.n.value), float(fit.ellip.value))

        # --- 1. b_n ---
        bn_a = float(b_n_approx(n))
        bn_e = float(gammaincinv(2 * n, 0.5))
        dbn = 100 * (bn_a - bn_e) / bn_e

        # --- 2. total flux: ours vs exact-bn closed-form vs numeric-20Re ---
        F_analytic = float(sersic_total_flux_analytic(amp, reff, n, ell))
        F_exact = float(2 * np.pi * n * Gamma(2 * n) * np.exp(bn_e) / bn_e ** (2 * n)
                        * amp * reff ** 2 * (1 - ell))
        half = int(max(400, 20 * reff))
        yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
        F_numeric = float(np.sum(fit(xx + fit.x_0.value, yy + fit.y_0.value)))
        # production model_full (rendered over the finite cutout, correct frame) — apcorr denom.
        x1, y1 = subimage_origin(b["cx"], b["cy"], b["img"].shape, 6.0, ps)
        model_full = render_full(fit, b["img"].shape, x1, y1)
        F_cutout = float(np.sum(model_full))
        cutout_loss = 100 * (F_numeric - F_cutout) / F_numeric

        # --- 3. enclosed-light within the MODEL's OWN ellipse at semi-major 2 R_e ---
        # (matches the analytic Sérsic law exactly for any ellipticity)
        mx, my = fit.x_0.value + x1, fit.y_0.value + y1   # model centre in full frame
        Re_pix = float(fit.r_eff.value)   # the MODEL's own r_eff (px) — match the analytic law
        thf = float(fit.theta.value)
        yyf, xxf = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
        ct, st = np.cos(-thf), np.sin(-thf)
        xr = (xxf - mx) * ct - (yyf - my) * st
        yr = (xxf - mx) * st + (yyf - my) * ct
        rell = np.hypot(xr, yr / max(1e-3, (1 - ell)))     # elliptical radius (semi-major units)
        encl_mod = float(np.sum(model_full[rell < NRE * Re_pix]) / F_numeric)
        encl_an = float(gammainc(2 * n, bn_e * NRE ** (1.0 / n)))   # Sérsic enclosed-light law

        # --- 4. aperture flux: ours (hard-edge) vs photutils EllipticalAperture ---
        xc, yc, a_pix, b_pix, theta = aperture_ellipse(b, NRE)
        ap = EllipticalAperture((xc, yc), a_pix, b_pix, theta=theta)
        # on the rendered MODEL (validates F_within summation method)
        ap_exact = float(aperture_photometry(model_full, ap, method='exact')['aperture_sum'][0])
        ap_ctr = float(aperture_photometry(model_full, ap, method='center')['aperture_sum'][0])
        ours_within = float(np.sum(model_full[in_aperture(b["img"].shape, xc, yc, a_pix, b_pix, theta)]))
        ap_exact_ctr = ap_exact / ap_ctr
        ap_over_ours = ap_exact / ours_within
        # on the real DATA image (validates F_raw) — sky-subtracted, masked pixels excluded
        rsky = np.hypot(xxf - xc, yyf - yc) * ps
        ring = (rsky > 2 * R_E_INIT) & ~mask
        sky = float(np.median(b["img"][ring])) if ring.any() else 0.0
        data = np.nan_to_num(b["img"]) - sky
        amask = in_aperture(b["img"].shape, xc, yc, a_pix, b_pix, theta) & ~mask
        Fraw_ours = float(np.sum(data[amask]))
        Fraw_exact = float(aperture_photometry(data, ap, mask=mask,
                                               method='exact')['aperture_sum'][0])

        rows[bn_] = dict(n=n, bn_a=bn_a, bn_e=bn_e, dbn=dbn,
                         F_analytic=F_analytic, F_exact=F_exact, F_numeric=F_numeric,
                         F_cutout=F_cutout, cutout_loss=cutout_loss,
                         encl_mod=encl_mod, encl_an=encl_an,
                         ap_exact_ctr=ap_exact_ctr, ap_over_ours=ap_over_ours,
                         Fraw_ours=Fraw_ours, Fraw_exact=Fraw_exact)
        print(f"{bn_:7s} {n:5.2f} {bn_a:8.4f} {bn_e:8.4f} {dbn:5.2f}% "
              f"{F_analytic/F_numeric:11.4f} {cutout_loss:10.2f}% {encl_mod:11.4f} "
              f"{encl_an:10.4f} {ap_exact_ctr:12.4f} {ap_over_ours:8.4f}")

    print("-" * 110)
    # worst-case summaries
    dbn_max = max(abs(rows[b]["dbn"]) for b in ORDER)
    an_num_max = max(abs(rows[b]["F_analytic"] / rows[b]["F_numeric"] - 1) for b in ORDER)
    encl_max = max(abs(rows[b]["encl_mod"] - rows[b]["encl_an"]) for b in ORDER)
    apraw_max = max(abs(rows[b]["Fraw_exact"] / rows[b]["Fraw_ours"] - 1) for b in ORDER)
    cutloss_max = max(rows[b]["cutout_loss"] for b in ORDER)
    print(f"  b_n approx error      : ≤ {dbn_max:.2f}%   (Ciotti&Bertin vs exact gammaincinv)")
    print(f"  analytic vs numeric   : ≤ {100*an_num_max:.2f}%   (our total-flux formula vs render-20R_e)")
    print(f"  enclosed-frac (model vs analytic gammainc) : ≤ {encl_max:.4f}  (Δ in enclosed fraction at 2 R_e)")
    print(f"  F_raw ours vs photutils exact-aperture     : ≤ {100*apraw_max:.2f}%  (hard-edge vs sub-pixel)")
    print(f"  cutout-truncation loss in sum(model_full)  : ≤ {cutloss_max:.2f}%  (apcorr denominator vs ∞)")
    print(f"    → mag bias from cutout truncation        : ≤ {2.5*np.log10(1/(1-cutloss_max/100)):.4f} mag")

    np.savez("results/validate_apcorr_established.npz",
             rows=np.array([(b, rows[b]) for b in ORDER], dtype=object),
             dbn_max=dbn_max, an_num_max=an_num_max, encl_max=encl_max,
             apraw_max=apraw_max, cutloss_max=cutloss_max)
    print("\n  → results/validate_apcorr_established.npz")


if __name__ == "__main__":
    main()
