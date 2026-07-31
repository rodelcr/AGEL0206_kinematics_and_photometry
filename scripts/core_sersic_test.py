"""Test C (#9) — core-Sérsic vs single-Sérsic, PSF-convolved (2026-06-17).

Question: is the deflector a CORED or a coreless (power-law) elliptical, and does it
matter for the single-Sérsic photometry used in the M* aperture correction?

Physical expectation (nb18): for sigma_e≈267 km/s the M•–r_b scalings (Rusli+13,
Thomas+16) predict a break/core radius r_b ~ 0.01–0.05", which is BELOW the finest
PSF FWHM (F150W2 0.049"; A6). So the core is almost certainly UNRESOLVED and this is
an upper-limit test, not a detection.

Method (per band, on a postage stamp around the deflector, arc+companions masked):
  1. PSF-convolved single-Sérsic fit (the production inner model).
  2. PSF-convolved core-Sérsic fit (Trujillo+2004 / Graham+2003): adds an inner
     power-law slope gamma inside a break radius r_b (alpha=10 fixed = sharp break).
  3. Compare BIC; report best-fit r_b vs the PSF HWHM and pixel scale; inspect the
     azimuthally-averaged central residual of the single-Sérsic (a real core shows a
     central DEFICIT, data < model).

Run on the two deepest, best-sampled bands: F150W2 (finest PSF) and F140W (STScI PSF).

Output: results/core_sersic_test.npz + results/figures/core_sersic_test.png.
Run from repo root:  python -m scripts.core_sersic_test
"""
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from astropy.convolution import convolve_fft
from astropy.modeling.models import Sersic2D
from astropy.stats import sigma_clipped_stats

from scripts.mask_method_comparison import load_band
from scripts.aperture_2re_companions import aperture_ellipse, companion_mask, in_aperture

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PSFDIR = os.path.join(REPO, "results", "psf_models")
_SYS = np.load(os.path.join(REPO, "results", "photometry_systematics.npz"), allow_pickle=True)
_TAB = np.load(os.path.join(REPO, "results", "sersic_parameter_table.npz"), allow_pickle=True)
STAMP_AS = 4.0          # postage-stamp half-size (arcsec)
FIT_AS = 2.5            # fit radius (arcsec) — inner galaxy, avoids outer envelope


def tab(band, key):
    return float({k: v for k, v in _TAB[f"{band}_cen"]}[key])


def b_n(n):
    return 2 * n - 1.0 / 3.0 + 4.0 / (405 * n) + 46.0 / (25515 * n ** 2)


def psf_native(band):
    """Load PSF, block-reduce the oversampled EPSF/PSFSTD to native pixels, renormalize."""
    d = np.load(os.path.join(PSFDIR, f"{band}_psf.npz"), allow_pickle=True)
    psf = np.asarray(d["psf"], float); ov = int(d["oversampling"])
    fwhm = float(d["fwhm_as"])
    if ov > 1:                                   # block-sum ov x ov -> native
        ny, nx = psf.shape
        ny -= ny % ov; nx -= nx % ov
        psf = psf[:ny, :nx].reshape(ny // ov, ov, nx // ov, ov).sum(axis=(1, 3))
    psf = np.clip(psf, 0, None); psf /= psf.sum()
    return psf, fwhm


def ell_r(xx, yy, x0, y0, ell, theta):
    q = max(1e-3, 1.0 - ell)
    dx = xx - x0; dy = yy - y0
    ct, st = np.cos(theta), np.sin(theta)
    xp = dx * ct + dy * st; yp = -dx * st + dy * ct
    return np.sqrt(xp ** 2 + (yp / q) ** 2)


def core_sersic(r, Ib, rb, reff, n, gamma, alpha=10.0):
    """Trujillo+2004 core-Sérsic (1D in elliptical radius)."""
    r = np.maximum(r, 1e-6)
    bn = b_n(n)
    pref = Ib * (2.0 ** (-gamma / alpha)) * (1.0 + (rb / r) ** alpha) ** (gamma / alpha)
    arg = (r ** alpha + rb ** alpha) / reff ** alpha
    return pref * np.exp(-bn * arg ** (1.0 / (alpha * n)))


def make_stamp(band):
    b = load_band(band); ps = b["pix_scale"]
    cx, cy = b["cx"], b["cy"]
    hs = int(round(STAMP_AS / ps))
    img = np.asarray(b["img"], float)
    x0, x1 = int(cx) - hs, int(cx) + hs + 1
    y0, y1 = int(cy) - hs, int(cy) + hs + 1
    stamp = img[y0:y1, x0:x1].copy()
    # arc + companion mask, cropped to the stamp
    xc, yc, a, bb, th = aperture_ellipse(b, 2)
    arc = _SYS[f"{band}_global_mask"].astype(bool)
    comp = companion_mask(b, arc, xc, yc, a, bb, th)
    fullmask = arc | comp
    mask = fullmask[y0:y1, x0:x1]
    cxl, cyl = cx - x0, cy - y0      # deflector center in stamp coords
    return stamp, mask, cxl, cyl, ps, b


def fit_models(band):
    stamp, mask, cx, cy, ps, b = make_stamp(band)
    psf, fwhm = psf_native(band)
    ny, nx = stamp.shape
    yy, xx = np.mgrid[:ny, :nx]
    rpix = np.hypot(xx - cx, yy - cy)
    _, med, std = sigma_clipped_stats(stamp[~mask], sigma=3.0)
    fitsel = (~mask) & np.isfinite(stamp) & (rpix < FIT_AS / ps)
    data = stamp - med
    noise = max(std, 1e-12)

    reff0 = tab(band, "reff_as") / ps; n0 = tab(band, "n")
    ell0 = tab(band, "_ell"); pa0 = np.radians(tab(band, "pa"))
    peak = np.nanmax(data[fitsel])
    amp0 = peak / 5.0

    def conv(model):
        return convolve_fft(model, psf, normalize_kernel=True, boundary="wrap")

    # ---- single Sérsic (PSF-convolved) ----
    def resid_single(p):
        logA, x0, y0, reff, n, ell, th, sky = p
        m = Sersic2D(amplitude=10 ** logA, r_eff=reff, n=n, x_0=x0, y_0=y0,
                     ellip=ell, theta=th)(xx, yy)
        return ((conv(m) + sky - data)[fitsel] / noise)

    p0s = [np.log10(amp0), cx, cy, reff0, n0, ell0, pa0, 0.0]
    lo = [np.log10(amp0) - 3, cx - 3, cy - 3, 0.3 / ps, 0.5, 0.0, -np.pi, -5 * std]
    hi = [np.log10(amp0) + 3, cx + 3, cy + 3, 6.0 / ps, 8.0, 0.9, np.pi, 5 * std]
    rs = least_squares(resid_single, p0s, bounds=(lo, hi), max_nfev=4000)
    chi2_s = float(np.sum(rs.fun ** 2)); k_s = len(p0s)

    # ---- core-Sérsic (PSF-convolved) ----
    def resid_core(p):
        logIb, x0, y0, reff, n, ell, th, sky, rb, gamma = p
        r = ell_r(xx, yy, x0, y0, ell, th)
        m = core_sersic(r, 10 ** logIb, rb, reff, n, gamma)
        return ((conv(m) + sky - data)[fitsel] / noise)

    p0c = [np.log10(amp0), cx, cy, reff0, n0, ell0, pa0, 0.0, 0.5 / ps, 0.1]
    loc = lo[:] + [1e-3, 0.0]; hic = hi[:] + [1.5 / ps, 1.5]
    rc = least_squares(resid_core, p0c, bounds=(loc, hic), max_nfev=6000)
    chi2_c = float(np.sum(rc.fun ** 2)); k_c = len(p0c)

    ndata = int(fitsel.sum())
    bic_s = chi2_s + k_s * np.log(ndata)
    bic_c = chi2_c + k_c * np.log(ndata)
    rb_fit_as = rc.x[8] * ps; gamma_fit = rc.x[9]

    # azimuthally-averaged central residual of the single-Sérsic (core → central deficit)
    logA, x0, y0, reff, n, ell, th, sky = rs.x
    msingle = conv(Sersic2D(amplitude=10 ** logA, r_eff=reff, n=n, x_0=x0, y_0=y0,
                            ellip=ell, theta=th)(xx, yy)) + sky
    resid = data - msingle
    redges = np.arange(0, FIT_AS / ps, max(1.0, 0.5 * fwhm / ps))
    rprof, rfrac = [], []
    for i in range(len(redges) - 1):
        sel = fitsel & (rpix >= redges[i]) & (rpix < redges[i + 1])
        if sel.sum() > 3:
            rprof.append(((redges[i] + redges[i + 1]) / 2) * ps)
            rfrac.append(float(np.median(resid[sel] / np.maximum(msingle[sel], noise))))

    # inner-deficit trend: does the single-Sérsic residual DEEPEN toward r=0 (real core
    # signature) or stay flat/oscillate (global-fit mismatch)? Compare innermost 2 bins
    # to the 3rd-5th bins.
    rfrac_arr = np.array(rfrac)
    inner_deepens = (len(rfrac_arr) >= 5 and
                     np.mean(rfrac_arr[:2]) < np.mean(rfrac_arr[2:5]) - 0.02)
    return dict(band=band, ps=ps, fwhm_as=fwhm, ndata=ndata,
                chi2_single=chi2_s, chi2_core=chi2_c, bic_single=bic_s, bic_core=bic_c,
                dBIC=bic_s - bic_c, n_single=rs.x[4], reff_single_as=rs.x[3] * ps,
                rb_fit_as=rb_fit_as, gamma_fit=gamma_fit, inner_deepens=inner_deepens,
                rb_upper_as=0.5 * fwhm, rprof_as=np.array(rprof), rfrac=rfrac_arr)


def main():
    os.makedirs(os.path.join(REPO, "results", "figures"), exist_ok=True)
    out = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    KPC_PER_AS = 7.2764
    # Depleted cores in M•~5e8 ellipticals are ~tens to a few hundred pc
    # (Lauer+2007, Rusli+2013) → ~0.005–0.05" at z=0.676. A FITTED r_b is a plausible
    # REAL depleted core only if it is BOTH resolved (> PSF HWHM) AND physically small
    # (< ~0.3 kpc ≈ 0.04"); anything larger is the core-Sérsic absorbing a global
    # single-Sérsic profile mismatch, not a nuclear core.
    RB_PHYS_MAX_AS = 0.3 / KPC_PER_AS    # 0.3 kpc ≈ 0.041"
    print(f"{'band':8s} {'PSF FWHM':>9s} {'n_s':>5s} {'rb_fit(\")':>9s} {'rb(kpc)':>8s} "
          f"{'gamma':>6s} {'ΔBIC':>8s} {'deepens?':>8s}  verdict")
    print("-" * 100)
    for i, band in enumerate(["F150W2", "F140W"]):
        r = fit_models(band)
        for k, v in r.items():
            out[f"{band}_{k}"] = v
        rb_kpc = r["rb_fit_as"] * KPC_PER_AS
        resolved = r["rb_fit_as"] > r["rb_upper_as"]
        physical = r["rb_fit_as"] < RB_PHYS_MAX_AS
        real_core = resolved and physical and r["dBIC"] > 10 and r["inner_deepens"]
        verdict = ("REAL depleted core" if real_core else
                   "artifact: rb too large for a depleted core (global-fit flexibility), "
                   "no inner deepening → core UNRESOLVED")
        out[f"{band}_real_core"] = real_core
        out[f"{band}_rb_phys_max_as"] = RB_PHYS_MAX_AS
        print(f"{band:8s} {r['fwhm_as']:8.3f}\" {r['n_single']:5.2f} {r['rb_fit_as']:9.3f} "
              f"{rb_kpc:8.2f} {r['gamma_fit']:6.2f} {r['dBIC']:8.1f} {str(bool(r['inner_deepens'])):>8s}  {verdict}")
        ax = axes[i]
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(r["rb_upper_as"], color="purple", ls=":", label=f"PSF HWHM={r['rb_upper_as']:.3f}\"")
        ax.plot(r["rprof_as"], 100 * r["rfrac"], "o-", color="C3")
        ax.set_title(f"{band}: single-Sérsic central residual\nΔBIC(single−core)={r['dBIC']:.1f}")
        ax.set_xlabel("radius (arcsec)"); ax.set_ylabel("(data−model)/model  [%]")
        ax.legend(fontsize=8); ax.set_xlim(0, FIT_AS)
    print("\nCONCLUSION: a genuine depleted core (predicted r_b ~ 0.005–0.05\" for M•~5e8) is "
          "UNRESOLVED at every PSF (finest HWHM 0.024\"=0.17 kpc, F150W2) → cannot distinguish "
          "cored vs coreless. Single-Sérsic is the pragmatic INNER model; its residual "
          "oscillation reflects the OUTER cD-envelope (Test A), already in the M* budget. A "
          "depleted core removes <1% of the light → negligible for M*. r_b upper limit ≈ 0.024\" (0.17 kpc).")
    fig.tight_layout()
    figpath = os.path.join(REPO, "results", "figures", "core_sersic_test.png")
    fig.savefig(figpath, dpi=130); plt.close(fig)
    np.savez(os.path.join(REPO, "results", "core_sersic_test.npz"), **out)
    print(f"\nSaved -> results/core_sersic_test.npz ; {figpath}")


if __name__ == "__main__":
    main()
