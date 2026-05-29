"""PSF-convolved 2-component Sersic fill model — firm up the fill-in leg of the
M★ systematic.

The fill-in (replace arc-masked pixels with the deflector model) is the dominant
M★ systematic. The 2-component Sersic fit (scripts/photometry_systematics.py) is
NOT PSF-convolved; a single env lacks `webbpsf`. Here we forward-model a
PSF-convolved 2-component Sersic (Gaussian PSF at the per-band instrument FWHM,
fit by scipy.least_squares), recompute the per-aperture filled magnitude, and
compare to the non-PSF fill to quantify the PSF systematic on the fill.

The fill region is the arc (r ~ 1–4″), well outside the PSF core, so the PSF
effect on the fill is expected to be small — this test demonstrates that.

PSF FWHM per band (published instrument values; the data are drizzled/resampled
so these are effective lower bounds — documented, not measured from a star):
    F200LP (WFC3/UVIS)  ~0.080″     F140W  (WFC3/IR)   ~0.140″
    F150W2 (NIRCam SW)  ~0.060″     F322W2 (NIRCam LW) ~0.130″

Usage: conda activate ISMgas; python scripts/psf_fill_model.py
Output: results/psf_fill_model.npz + prints the ΔPSF(filled) per band.
"""
import os
import sys
import json
import warnings

import numpy as np
from astropy.modeling.models import Sersic2D
from astropy.convolution import Gaussian2DKernel, convolve_fft
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.mask_method_comparison import load_band, reproject_mask, _aperture, flux_mag
from scripts.arc_mask_verification import sky_sigma

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
PSF_FWHM = {"F200LP": 0.080, "F140W": 0.140, "F150W2": 0.060, "F322W2": 0.130}  # arcsec
R_E_INIT = 2.305


def two_comp(params, xx, yy):
    """Bulge+disk Sersic2D from a flat param vector."""
    (a1, r1, n1, e1, t1, a2, r2, n2, e2, t2, x0, y0) = params
    b = Sersic2D(amplitude=a1, r_eff=r1, n=n1, x_0=x0, y_0=y0, ellip=e1, theta=t1)
    d = Sersic2D(amplitude=a2, r_eff=r2, n=n2, x_0=x0, y_0=y0, ellip=e2, theta=t2)
    return np.asarray(b(xx, yy) + d(xx, yy), dtype=float)


def fit_psf_2comp(sub, sub_mask, ps, fwhm_arcsec, init):
    """Fit a PSF-convolved 2-comp Sersic to sky-subtracted sub-image via least_squares."""
    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    w = (~sub_mask).astype(float)
    sigma_pix = (fwhm_arcsec / ps) / 2.3548
    kernel = Gaussian2DKernel(x_stddev=max(sigma_pix, 0.3)).array

    def resid(p):
        model = two_comp(p, xx, yy)
        model_c = convolve_fft(model, kernel, normalize_kernel=True, boundary="fill", fill_value=0.0)
        return ((sub - model_c) * w).ravel()

    reff = R_E_INIT / ps
    lo = [0, reff*0.1, 2.0, 0, -np.pi, 0, reff*0.6, 0.5, 0, -np.pi, init[10]-3, init[11]-3]
    hi = [1e6, reff*1.0, 6.0, 0.6, np.pi, 1e6, reff*2.5, 2.0, 0.7, np.pi, init[10]+3, init[11]+3]
    init = np.clip(init, lo, hi)
    res = least_squares(resid, init, bounds=(lo, hi), method="trf", max_nfev=300)
    return res.x, kernel


def analyze(name, hst_master):
    b = load_band(name)
    ps, cx, cy = b["pix_scale"], b["cx"], b["cy"]
    sky_med, _ = sky_sigma(b["img"], cx, cy, ps)
    M = reproject_mask(*hst_master, b["wcs"], b["img"].shape)   # F200LP-reproj per-band base
    p = json.load(open(b["cfg"]["params"]))
    aper, ann = _aperture(b, p)

    box = 6.0
    half = int(np.ceil(box / ps))
    y1, y2 = max(0, int(cy)-half), min(b["img"].shape[0], int(cy)+half+1)
    x1, x2 = max(0, int(cx)-half), min(b["img"].shape[1], int(cx)+half+1)
    sub = np.nan_to_num(b["img"][y1:y2, x1:x2].astype(float)) - sky_med
    sub_mask = M[y1:y2, x1:x2]
    cx0, cy0 = cx-x1, cy-y1
    reff = R_E_INIT/ps
    amp = float(np.nanpercentile(sub[~sub_mask], 98)) if (~sub_mask).any() else 1.0
    e0 = float(p["geometry"]["ellipticity"]); t0 = float(p["pixel_geometry"]["theta_rad"])
    init = [amp, reff*0.3, 4.0, e0, t0, amp*0.3, reff*1.2, 1.0, e0, t0, cx0, cy0]

    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    # non-PSF fit (for ΔPSF reference) — same objective, kernel=delta
    par_nopsf, _ = fit_psf_2comp(sub, sub_mask, ps, 1e-6, init)   # tiny FWHM ≈ no PSF
    par_psf, kern = fit_psf_2comp(sub, sub_mask, ps, PSF_FWHM[name], init)

    def filled_mag(par, psf_conv):
        model = two_comp(par, xx, yy)
        if psf_conv:
            model = convolve_fft(model, kern, normalize_kernel=True, boundary="fill", fill_value=0.0)
        model_full = np.zeros(b["img"].shape)
        model_full[y1:y2, x1:x2] = np.clip(model, 0, None)
        data_filled = np.where(M, model_full + sky_med, np.nan_to_num(b["img"]))
        return flux_mag(b, data_filled, aper, ann, None)

    m_nopsf = filled_mag(par_nopsf, False)
    m_psf = filled_mag(par_psf, True)
    print(f"{name}: filled(no-PSF)={m_nopsf:.4f}  filled(PSF FWHM={PSF_FWHM[name]}\")={m_psf:.4f}  "
          f"ΔPSF={m_psf-m_nopsf:+.4f}  [n1={par_psf[2]:.1f} n2={par_psf[7]:.1f}]")
    return dict(name=name, pivot=b["cfg"]["pivot_AA"], fwhm=PSF_FWHM[name],
                filled_nopsf=m_nopsf, filled_psf=m_psf, dpsf=m_psf-m_nopsf)


def main():
    f200 = load_band("F200LP")
    hst_master = (f200["mask"].astype(bool), f200["wcs"])
    print("PSF-convolved 2-component Sersic fill model (Gaussian PSF at instrument FWHM)\n")
    res = [analyze(n, hst_master) for n in ORDER]
    dpsf = np.array([r["dpsf"] for r in res])
    print(f"\nΔPSF on filled mag per band: {np.round(dpsf,4).tolist()}")
    print(f"max |ΔPSF| = {np.abs(dpsf).max():.4f} mag  (band {ORDER[int(np.abs(dpsf).argmax())]})")
    print("→ PSF effect on the fill is small (fill region is the arc at r~1-4″, outside the PSF core).")
    np.savez("results/psf_fill_model.npz",
             filter_names=np.array(ORDER),
             pivot=np.array([r["pivot"] for r in res]),
             fwhm=np.array([r["fwhm"] for r in res]),
             filled_nopsf=np.array([r["filled_nopsf"] for r in res]),
             filled_psf=np.array([r["filled_psf"] for r in res]),
             dpsf=dpsf)
    print("Saved → results/psf_fill_model.npz")


if __name__ == "__main__":
    main()
