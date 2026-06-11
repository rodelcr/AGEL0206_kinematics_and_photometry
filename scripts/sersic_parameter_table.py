"""Single-component Sérsic parameters per band — appendix table (2026-06-11).

Fits the production single-Sérsic deflector model (`sersic_total_photometry.fit_sersic2d`,
the same fitter used for the aperture-corrected photometry) on the BEST (global
color/morph + WCS-reg) mask, for all four bands (F200LP, F140W, F150W2, F322W2), and
tabulates the publication parameters:

  r_eff  [″ and kpc]   effective (half-light) radius
  n                    Sérsic index
  b/a = 1 − ellip      axis ratio
  PA   [deg E of N]    major-axis position angle, via the band WCS
  μ_e  [mag/arcsec²]   surface brightness at r_eff (AB) = -2.5 log10(I_e/pixarea) + ZP
  m_tot[AB]            Sérsic total magnitude (integrated to ∞; Graham & Driver 2005)

Uncertainties: parametric bootstrap — add residual-scaled Gaussian noise to the image
and re-fit N times; 1σ = half the 16–84 percentile span. (Avoids the railed-bound
covariance pathologies of LevMar; the ellipticity in particular can sit near a bound.)

Outputs: results/sersic_parameter_table.npz, results/sersic_parameter_table.md (markdown
appendix table), results/sersic_parameter_table.tex (deluxetable for the manuscript).
Usage: conda activate ISMgas; python scripts/sersic_parameter_table.py [--n-boot 80]
"""
import os, sys, argparse
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band
from scripts.sersic_total_photometry import (fit_sersic2d, sersic_total_flux_analytic,
                                             R_E_INIT)

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
SYSZ = np.load("results/photometry_systematics.npz", allow_pickle=True)
KPC_PER_ARCSEC = 7.04   # at z=0.67564
PIVOT = {n: float(load_band(n)["cfg"]["pivot_AA"]) for n in ORDER}


def sky_pa_deg(wcs, xc, yc, theta_img, ps):
    """Position angle (deg E of N) of the major axis, from the image-frame theta."""
    step = 5.0 / ps  # 5" along the major axis
    x2, y2 = xc + step * np.cos(theta_img), yc + step * np.sin(theta_img)
    (ra0, dec0), (ra1, dec1) = (wcs.pixel_to_world_values(xc, yc),
                                wcs.pixel_to_world_values(x2, y2))
    dra = (ra1 - ra0) * np.cos(np.radians(dec0))
    ddec = (dec1 - dec0)
    pa = np.degrees(np.arctan2(dra, ddec)) % 180.0   # E of N, mod 180 (axis, not vector)
    return float(pa)


def fit_params(b, mask, ps, zp, wcs):
    fit, _, info = fit_sersic2d(b["img"], mask, (b["cx"], b["cy"]), R_E_INIT, 6.0, ps)
    amp, reff, n, ell, th = (float(fit.amplitude.value), float(fit.r_eff.value),
                             float(fit.n.value), float(fit.ellip.value),
                             float(fit.theta.value))
    F = sersic_total_flux_analytic(amp, reff, n, ell)
    m_tot = -2.5 * np.log10(F) + zp
    mu_e = -2.5 * np.log10(amp / ps**2) + zp        # I_e per pixel → per arcsec²
    pa = sky_pa_deg(wcs, b["cx"], b["cy"], th, ps)
    return dict(reff_as=reff * ps, reff_kpc=reff * ps * KPC_PER_ARCSEC, n=n,
                ba=1 - ell, pa=pa, mu_e=mu_e, m_tot=m_tot,
                res_med=info["res_med_pct"], _reff_px=reff, _amp=amp, _ell=ell)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=80)
    args = ap.parse_args()
    rng = np.random.default_rng(20260611)

    rows = {}
    for bn in ORDER:
        b = load_band(bn); ps = b["pix_scale"]; zp = b["ab_zp"]
        wcs = b["wcs"]
        mask = SYSZ[f"{bn}_global_mask"].astype(bool)
        cen = fit_params(b, mask, ps, zp, wcs)

        # residual scatter for the parametric bootstrap (sky-sub, unmasked, in the box)
        fit0, (x1, y1), info0 = fit_sersic2d(b["img"], mask, (b["cx"], b["cy"]),
                                             R_E_INIT, 6.0, ps)
        resid = info0["residual"][~info0["sub_mask"]]
        sig = float(np.std(resid[np.isfinite(resid)]))
        half = int(np.ceil(6.0 / ps))
        boot = {k: [] for k in ("reff_as", "reff_kpc", "n", "ba", "pa", "mu_e", "m_tot")}
        for _ in range(args.n_boot):
            img2 = b["img"].copy().astype(float)
            y2 = min(img2.shape[0], int(b["cy"]) + half + 1)
            x2 = min(img2.shape[1], int(b["cx"]) + half + 1)
            yy0, xx0 = max(0, int(b["cy"]) - half), max(0, int(b["cx"]) - half)
            noise = rng.normal(0, sig, size=img2[yy0:y2, xx0:x2].shape)
            img2[yy0:y2, xx0:x2] += noise
            bb = dict(b); bb["img"] = img2
            try:
                p = fit_params(bb, mask, ps, zp, wcs)
                for k in boot:
                    boot[k].append(p[k])
            except Exception:
                continue
        err = {k: float((np.percentile(v, 84) - np.percentile(v, 16)) / 2)
               if len(v) > 5 else float("nan") for k, v in boot.items()}
        # The parametric bootstrap captures only pixel noise → unrealistically tight on
        # shape (b/a±0.00, PA±1°). Floor the shape errors with the injection-recovery
        # SCATTER from scripts/validate_sersic_fitter_synthetic.py (the realistic
        # systematic for these depths/S/N): b/a ±0.06, PA ±20°, n ±0.10.
        try:
            syn = np.load("results/validate_sersic_fitter_synthetic.npz", allow_pickle=True)
            ba_floor = float(syn["ba_scat"]); n_floor = float(syn["n_scat"])
        except Exception:
            ba_floor, n_floor = 0.06, 0.10
        err["ba"] = max(err["ba"], ba_floor)
        err["pa"] = max(err["pa"], 20.0)
        err["n"] = max(err["n"], min(n_floor, 0.15))
        rows[bn] = dict(cen=cen, err=err, n_boot_ok=len(boot["n"]),
                        circular=bool(cen["ba"] > 0.97))
        print(f"  {bn}: r_eff={cen['reff_as']:.3f}±{err['reff_as']:.3f}\"  "
              f"n={cen['n']:.2f}±{err['n']:.2f}  b/a={cen['ba']:.2f}±{err['ba']:.2f}  "
              f"PA={cen['pa']:.0f}±{err['pa']:.0f}°  m={cen['m_tot']:.2f}±{err['m_tot']:.2f}  "
              f"(boot {len(boot['n'])}/{args.n_boot})")

    def ba_str(bn, tex=False):
        c, e = rows[bn]["cen"], rows[bn]["err"]
        pm = r"\pm" if tex else "±"
        if rows[bn]["circular"]:
            return r"$\gtrsim$0.95$^{\dagger}$" if tex else "≳0.95†"
        return f"${c['ba']:.2f}{pm}{e['ba']:.2f}$" if tex else f"{c['ba']:.2f}{pm}{e['ba']:.2f}"

    def pa_str(bn, tex=False):
        c, e = rows[bn]["cen"], rows[bn]["err"]
        if rows[bn]["circular"]:
            return r"\nodata" if tex else "—†"
        pm = r"\pm" if tex else "±"
        return f"${c['pa']:.0f}{pm}{e['pa']:.0f}$" if tex else f"{c['pa']:.0f}{pm}{e['pa']:.0f}"

    # ---- markdown table ----
    md = ["| Band | $\\lambda_{\\rm piv}$ (Å) | $r_e$ (″) | $r_e$ (kpc) | $n$ | $b/a$ | "
          "PA (°) | $\\mu_e$ (mag/″²) | $m_{\\rm AB}^{\\rm tot}$ |",
          "|---|---|---|---|---|---|---|---|---|"]
    for bn in ORDER:
        c, e = rows[bn]["cen"], rows[bn]["err"]
        md.append(f"| {bn} | {PIVOT[bn]:.0f} | {c['reff_as']:.3f}±{e['reff_as']:.3f} | "
                  f"{c['reff_kpc']:.2f}±{e['reff_kpc']:.2f} | {c['n']:.2f}±{e['n']:.2f} | "
                  f"{ba_str(bn)} | {pa_str(bn)} | "
                  f"{c['mu_e']:.2f} | {c['m_tot']:.2f}±{e['m_tot']:.2f} |")
    md.append("")
    md.append("† F200LP: the deflector is faint in this band and the single-Sérsic ellipticity "
              "is unconstrained (consistent with circular); the IR bands set the shape.")
    md_txt = "\n".join(md)
    open("results/sersic_parameter_table.md", "w").write(md_txt + "\n")

    # ---- LaTeX deluxetable ----
    tex = [r"\begin{deluxetable}{lcccccccc}",
           r"\tablecaption{Single-component S\'ersic parameters of the AGEL\,J020613$-$011417 "
           r"deflector, fit per band on the best (color/morph-gated, WCS-registered) arc mask. "
           r"Uncertainties are 1$\sigma$ from a parametric bootstrap, with the $b/a$, PA, and $n$ "
           r"errors floored by the injection-recovery scatter (\S validation); the mask-choice "
           r"systematic (Sec.~2.1.1b) dominates and is not included. $^{\dagger}$F200LP: the deflector "
           r"is too faint here to constrain the ellipticity (consistent with circular).\label{tab:sersic}}",
           r"\tablehead{\colhead{Band} & \colhead{$\lambda_{\rm piv}$} & \colhead{$r_e$} & "
           r"\colhead{$r_e$} & \colhead{$n$} & \colhead{$b/a$} & \colhead{PA} & "
           r"\colhead{$\mu_e$} & \colhead{$m_{\rm AB}^{\rm tot}$} \\ "
           r"\colhead{} & \colhead{(\AA)} & \colhead{($''$)} & \colhead{(kpc)} & \colhead{} & "
           r"\colhead{} & \colhead{(deg)} & \colhead{(mag\,arcsec$^{-2}$)} & \colhead{(AB)}}",
           r"\startdata"]
    for bn in ORDER:
        c, e = rows[bn]["cen"], rows[bn]["err"]
        tex.append(f"{bn} & {PIVOT[bn]:.0f} & ${c['reff_as']:.3f}\\pm{e['reff_as']:.3f}$ & "
                   f"${c['reff_kpc']:.2f}\\pm{e['reff_kpc']:.2f}$ & ${c['n']:.2f}\\pm{e['n']:.2f}$ & "
                   f"{ba_str(bn, tex=True)} & {pa_str(bn, tex=True)} & "
                   f"${c['mu_e']:.2f}$ & ${c['m_tot']:.2f}\\pm{e['m_tot']:.2f}$ \\\\")
    tex += [r"\enddata", r"\end{deluxetable}"]
    open("results/sersic_parameter_table.tex", "w").write("\n".join(tex) + "\n")

    np.savez("results/sersic_parameter_table.npz",
             bands=np.array(ORDER),
             **{f"{bn}_{kind}": np.array(list(rows[bn][kind].items()), dtype=object)
                for bn in ORDER for kind in ("cen", "err")})
    print("\n" + md_txt)
    print("\n  → results/sersic_parameter_table.{md,tex,npz}")


if __name__ == "__main__":
    main()
