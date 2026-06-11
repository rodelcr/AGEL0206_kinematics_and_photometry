"""Sersic2D fits per band → analytic total flux → AB mag.

For nb08: replace the empirical-aperture photometry inputs with band-specific
2D Sersic-model total fluxes, so we can quantify the aperture systematic on
the Bagpipes-derived stellar mass.

Bands: HST F140W, HST F200LP, JWST F150W2, JWST F322W2
Output: results/sersic_total_photometry.npz with per-band:
    fit params (n, r_eff, ellip, theta, amplitude, x_0, y_0)
    total flux (analytic + numeric cross-check)
    AB magnitude
    aperture (current nb02) AB mag for comparison

Usage
-----
    conda activate ISMgas
    python scripts/sersic_total_photometry.py
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Sersic2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.centroids import centroid_2dg
from scipy.special import gamma as Gamma

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

# Paths
HST_F140W = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits"
HST_F140W_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits"
HST_F200LP = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits"
HST_F200LP_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits"
JWST_F150W2 = "../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f150w2_i2d.fits"
JWST_F150W2_MASK = "../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f150w2_i2d_mask.fits"
JWST_F322W2 = "../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f322w2_i2d.fits"
JWST_F322W2_MASK = "../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f322w2_i2d_mask.fits"

RA_DEFL, DEC_DEFL = 31.55611, -1.23817
R_E_INIT = 2.305  # paper headline R_e in arcsec — initial Sersic2D guess
FIG_DIR = "results/figures"


def b_n(n):
    """Ciotti & Bertin 1999 approximation to b_n in the Sersic profile.

    Defined such that γ(2n, b_n) = Γ(2n)/2 (i.e., r_eff = half-light radius).
    """
    return 2 * n - 1.0 / 3.0 + 0.009876 / n


def sersic_total_flux_analytic(amplitude, r_eff_pix, n, ellip):
    """Total flux integrated to infinity from Sersic2D parameters.

    Sersic2D's `amplitude` is I(r=r_eff). The analytic 2D total integral is

        F_tot = 2π · n · Γ(2n) · b_n^(-2n) · I_e · r_eff^2 · (1 - ellip) · exp(b_n)

    in units of (image data unit) × (pixel area in pixel² units). The factor
    (1 - ellip) accounts for the elliptical aperture (semi-minor/semi-major).

    Reference: Graham & Driver (2005) eq. 4.
    """
    bn = b_n(n)
    return (2 * np.pi * n * Gamma(2 * n) * np.exp(bn) / bn ** (2 * n)
            * amplitude * r_eff_pix ** 2 * (1 - ellip))


def sersic_total_flux_numeric(fit, box_pix=400, n_grid=2001):
    """Numeric integral on a wide grid, for cross-check with analytic."""
    r_eff_pix = fit.r_eff.value
    half = max(box_pix, int(20 * r_eff_pix))
    yy, xx = np.mgrid[-half:half+1, -half:half+1]
    # Center at origin for a clean integral (Sersic2D evaluated about its own x_0,y_0)
    xx_eval = xx + fit.x_0.value
    yy_eval = yy + fit.y_0.value
    return float(np.sum(fit(xx_eval, yy_eval)))


def moment_seed(sub, weights, cx, cy, pix_scale, rmax_arcsec):
    """Flux-weighted second-moment estimate of (ellip, theta) for a Sérsic-fit seed.

    Standard SExtractor/GALFIT/statmorph-style pre-estimation: it puts the LevMar fit
    in the elliptical basin (a low-ellip init otherwise rails into the spurious circular
    minimum). The moment ellipticity is biased low (the round bright core dominates the
    flux weighting), so it is a SEED only — the Sérsic fit refines it. Returns
    (ellip0 ∈ [0,0.7], theta0 rad). Falls back to (0.2, 0.0) if the window is empty."""
    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    r = np.hypot(xx - cx, yy - cy) * pix_scale
    f = np.clip(sub, 0, None) * weights * (r < rmax_arcsec)
    tot = f.sum()
    if tot <= 0:
        return 0.2, 0.0
    Mxx = (f * (xx - cx) ** 2).sum() / tot
    Myy = (f * (yy - cy) ** 2).sum() / tot
    Mxy = (f * (xx - cx) * (yy - cy)).sum() / tot
    theta = 0.5 * np.arctan2(2 * Mxy, Mxx - Myy)
    s, d = Mxx + Myy, np.hypot(Mxx - Myy, 2 * Mxy)
    ba = np.sqrt(max((s - d) / 2, 1e-6) / max((s + d) / 2, 1e-6))
    return float(np.clip(1 - ba, 0.0, 0.7)), float(theta)


def fit_sersic2d(img, mask, center, r_eff_init_arcsec, box_arcsec, pix_scale):
    """Sersic2D fit with sky subtraction and flux-weighted least squares.

    Steps:
      1. Cut a box of `box_arcsec` around `center`.
      2. Estimate sky background from the outer annulus (R > 2 R_e_init within box,
         excluding mask). Subtract from sub-image.
      3. Initial amplitude from sky-subtracted center value, scaled assuming n=2 Sersic.
      4. Fit Sersic2D with weights = (~mask) × sqrt(data_clip), so bright deflector
         pixels dominate over sky pixels.
    """
    xc, yc = center
    half = int(np.ceil(box_arcsec / pix_scale))
    y1 = max(0, int(yc) - half)
    y2 = min(img.shape[0], int(yc) + half + 1)
    x1 = max(0, int(xc) - half)
    x2 = min(img.shape[1], int(xc) + half + 1)
    sub_raw = img[y1:y2, x1:x2].astype(float)
    sub_mask = mask[y1:y2, x1:x2]
    yy, xx = np.mgrid[:sub_raw.shape[0], :sub_raw.shape[1]]
    cx_box = xc - x1
    cy_box = yc - y1
    r_box = np.hypot(xx - cx_box, yy - cy_box) * pix_scale  # in arcsec

    # Sky from outer annulus (R > 2 × R_e_init), excluding mask + handling NaNs
    sky_ring = (r_box > 2.0 * r_eff_init_arcsec) & ~sub_mask & np.isfinite(sub_raw)
    if sky_ring.sum() > 0:
        sky_bkg = float(np.median(sub_raw[sky_ring]))
    else:
        sky_bkg = 0.0
    sub = np.nan_to_num(sub_raw - sky_bkg)

    # Initial amplitude from sky-subtracted central peak, scaled for n=2 Sersic.
    # For a Sersic with index n, I(r=r_eff) / I(0) = exp(-b_n) ≈ exp(-3.67) ≈ 0.0254 at n=2.
    cy_int, cx_int = int(round(cy_box)), int(round(cx_box))
    cy_int = np.clip(cy_int, 0, sub.shape[0]-1)
    cx_int = np.clip(cx_int, 0, sub.shape[1]-1)
    # 3x3 median around center to suppress noise
    cy0, cy1 = max(0, cy_int-1), min(sub.shape[0], cy_int+2)
    cx0, cx1 = max(0, cx_int-1), min(sub.shape[1], cx_int+2)
    peak = float(np.median(sub[cy0:cy1, cx0:cx1]))
    if peak <= 0:
        peak = max(float(np.percentile(sub[~sub_mask & np.isfinite(sub_raw)], 99.5)), 1e-6)
    # Physical r_eff for AGEL0206 deflector is ~1-2.3" (R_e = 2.305 = HST CoG total).
    # Sersic-fit r_eff often comes out smaller than the CoG R_e because Sersic fits
    # are biased low for n>1 profiles. Bounds: 0.3" to 1.3 × R_e_init.
    r_eff_pix = r_eff_init_arcsec / pix_scale
    r_eff_min_pix = 0.3 / pix_scale
    r_eff_max_pix = r_eff_pix * 2.0  # up to 2 × R_e_init = 4.6"
    amp_init = peak * 0.05  # roughly I(r_eff) for n≈3 (b_n exp ≈ 0.05)

    _bounds = {"n": (0.5, 6.0),
               "r_eff": (r_eff_min_pix, r_eff_max_pix),
               "ellip": (0.0, 0.7),
               "amplitude": (peak * 1e-3, peak * 1.5),
               "x_0": (cx_box - 3, cx_box + 3),
               "y_0": (cy_box - 3, cy_box + 3)}

    # Plain unmasked weight, on the SKY-SUBTRACTED data. Sky is now ≈0 so the
    # fitter doesn't drown the bright deflector signal in sky pixels.
    weights = (~sub_mask).astype(float)

    # MULTI-START — the single-Sérsic χ² has a spurious CIRCULAR local minimum (ellip→0)
    # alongside the true elliptical one; a low-ellip init rails into the circular trap
    # (b/a=1 for all bands), under-reporting the deflector shape AND biasing the total
    # flux via the (1-ellip) factor. Seed from flux-weighted 2nd moments (SExtractor/GALFIT
    # style), plus a boosted-ellip start (escape the circular basin) and a circular fallback
    # (faint bands may genuinely prefer it); keep the lowest weighted-residual fit.
    # (2026-06-11; was a single ellip=0.15 start.)
    e_mom, th_mom = moment_seed(sub, weights, cx_box, cy_box, pix_scale, 1.5 * r_eff_init_arcsec)
    starts = ([(e_mom, th_mom)]                                   # moment seed
              + [(0.4, t) for t in np.linspace(0, np.pi, 6, endpoint=False)]  # PA grid
              + [(0.0, 0.0)])                                     # circular fallback
    best_fit, best_rss = None, np.inf
    for e0, th0 in starts:
        sersic_init = Sersic2D(amplitude=amp_init, r_eff=r_eff_pix * 0.6, n=3.0,
                               x_0=cx_box, y_0=cy_box, ellip=min(e0, 0.69), theta=th0,
                               bounds=_bounds)
        cand = LevMarLSQFitter()(sersic_init, xx, yy, sub, weights=weights, maxiter=1500)
        rss = float(np.sum((weights * (sub - cand(xx, yy))) ** 2))
        if rss < best_rss:
            best_fit, best_rss = cand, rss
    fit = best_fit

    # Residual statistics on unmasked pixels (sky-subtracted scale)
    model = fit(xx, yy)
    residual = sub - model
    if (~sub_mask).any():
        denom = np.maximum(np.abs(sub[~sub_mask]), peak * 0.01)
        res_pct = np.percentile(np.abs(residual[~sub_mask]) / denom, [50, 95])
    else:
        res_pct = (None, None)

    return fit, (x1, y1), dict(sub=sub, sub_raw=sub_raw, sub_mask=sub_mask,
                                model=model, residual=residual,
                                res_med_pct=res_pct[0], res_p95_pct=res_pct[1],
                                sky_bkg=sky_bkg, amp_init=amp_init, peak=peak)


def find_center_2dg(img, mask, wcs, ra0, dec0, box_arcsec, pix_scale):
    """Find centroid via 2D Gaussian fit, with mask. Falls back to WCS guess."""
    xc0, yc0 = wcs.world_to_pixel_values(ra0, dec0)
    xc0, yc0 = float(xc0), float(yc0)
    half = int(np.ceil(box_arcsec / pix_scale))
    y1 = max(0, int(yc0) - half)
    y2 = min(img.shape[0], int(yc0) + half + 1)
    x1 = max(0, int(xc0) - half)
    x2 = min(img.shape[1], int(xc0) + half + 1)
    sub = np.nan_to_num(img[y1:y2, x1:x2])
    sub_m = mask[y1:y2, x1:x2] if mask is not None else None
    try:
        cx_rel, cy_rel = centroid_2dg(sub, mask=sub_m)
        if np.isfinite(cx_rel) and np.isfinite(cy_rel):
            return float(x1 + cx_rel), float(y1 + cy_rel), "2dg"
    except Exception:
        pass
    return xc0, yc0, "wcs_fallback"


# --------------------------------------------------------------------------
# Per-band loaders + flux→AB-mag converters
# --------------------------------------------------------------------------
def load_hst(image_path, mask_path, photflam, photplam):
    """HST WFC3 cutout — drizzled, units e/s.

    PHOTFLAM and PHOTPLAM are passed in because the cutouts have stripped
    the calibration metadata. Values from CLAUDE.md (also Pivot wavelengths
    confirmed against nb02 photometry):
        F200LP: PHOTFLAM=5.1851e-20, PHOTPLAM=4923.48 (ZP_AB=27.344)
        F140W:  PHOTFLAM=1.4737e-20, PHOTPLAM=13922.9 (ZP_AB=26.452)
    """
    with fits.open(image_path) as h:
        img = h[0].data.astype(float)
        hdr = h[0].header
    wcs = WCS(hdr)
    mask = fits.getdata(mask_path).astype(bool)
    pix_scale = float(np.abs(proj_plane_pixel_scales(wcs)[0])) * 3600
    # AB ZP for count rate (e/s):
    zp_st = -2.5 * np.log10(photflam) - 21.10
    zp_ab = zp_st - 5.0 * np.log10(photplam) + 18.6921
    return dict(img=img, mask=mask, wcs=wcs, pix_scale=pix_scale,
                zp_ab=zp_ab, photflam=photflam, photplam=photplam,
                unit="e/s", flux_to_jy=None)


def load_jwst(image_path, mask_path):
    """JWST NIRCam i2d — units MJy/sr."""
    with fits.open(image_path) as h:
        img = h["SCI"].data.astype(float)
        hdr_sci = h["SCI"].header
        hdr_pri = h[0].header
    wcs = WCS(hdr_sci)
    mask = fits.getdata(mask_path).astype(bool)
    if mask.shape != img.shape:
        raise RuntimeError(f"JWST mask shape {mask.shape} != image shape {img.shape}")
    pix_scale = float(np.abs(proj_plane_pixel_scales(wcs)[0])) * 3600
    pixar_sr = float(hdr_sci["PIXAR_SR"])
    # AB ZP per pixel (data unit MJy/sr × PIXAR_SR steradians = MJy/pix):
    # F_Jy = pixel_value_MJy_per_sr * PIXAR_SR * 1e6
    # mag = -2.5 log10(F_Jy / 3631) = -2.5 log10(pix_val) - 2.5 log10(PIXAR_SR * 1e6 / 3631)
    #     = -2.5 log10(pix_val) - 2.5 log10(PIXAR_SR) - 6.1003
    zp_ab = -6.10 - 2.5 * np.log10(pixar_sr)
    return dict(img=img, mask=mask, wcs=wcs, pix_scale=pix_scale,
                zp_ab=zp_ab, pixar_sr=pixar_sr,
                unit="MJy/sr", flux_to_jy=pixar_sr * 1e6)


def total_flux_to_ab(total_flux_per_pix, band):
    """Convert total summed flux (in image's pixel native unit, summed) to AB mag.

    For HST e/s: mag = -2.5 log10(sum_e_per_s) + ZP_AB
    For JWST MJy/sr: this needs care. The Sersic integral is in MJy/sr × pix² because
        F_tot = ∫ I(x,y) dxdy where (x,y) are in pix and I is per-pix (= MJy/sr).
        The actual flux in Jy is sum(MJy/sr × pix²) × (pix_in_sr). Since dxdy integral
        in pixels gives "pix² units" and each pixel covers PIXAR_SR steradians, the
        result needs × PIXAR_SR * 1e6 for Jy.
        AB = -2.5 log10(F_Jy / 3631) = -2.5 log10(F_Jy) + 8.9
        So AB = -2.5 log10(sum × PIXAR_SR × 1e6) + 8.9
              = -2.5 log10(sum) - 2.5 log10(PIXAR_SR) - 15 + 8.9
              = -2.5 log10(sum) + zp_ab_perpix
    """
    return -2.5 * np.log10(total_flux_per_pix) + band["zp_ab"]


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    print("=" * 75)
    print("Sersic2D total photometry — 4 bands")
    print("=" * 75)

    # PHOTFLAM, PHOTPLAM from CLAUDE.md (stripped from the cutouts; use canonical values)
    bands = {
        "F200LP": dict(
            loader=lambda: load_hst(HST_F200LP, HST_F200LP_MASK,
                                    photflam=5.1851e-20, photplam=4923.48),
            pivot_AA=4972.0, ap_mag=22.6126, ap_err=0.0172,
        ),
        "F140W":  dict(
            loader=lambda: load_hst(HST_F140W, HST_F140W_MASK,
                                    photflam=1.4737e-20, photplam=13922.9),
            pivot_AA=13923.0, ap_mag=19.1335, ap_err=0.0044,
        ),
        "F150W2": dict(
            loader=lambda: load_jwst(JWST_F150W2, JWST_F150W2_MASK),
            pivot_AA=16720.0, ap_mag=18.9425, ap_err=0.0005,
        ),
        "F322W2": dict(
            loader=lambda: load_jwst(JWST_F322W2, JWST_F322W2_MASK),
            pivot_AA=32470.0, ap_mag=18.6042, ap_err=0.0003,
        ),
    }

    results = {}
    for name, cfg in bands.items():
        print(f"\n--- {name} ---")
        b = cfg["loader"]()
        print(f"  unit={b['unit']}, pix_scale={b['pix_scale']:.4f}\"/pix, "
              f"shape={b['img'].shape}, mask sum={b['mask'].sum()}")
        print(f"  AB ZP per pixel = {b['zp_ab']:.4f}")

        # Find centroid in this band's pixel grid
        cx, cy, method = find_center_2dg(b["img"], b["mask"], b["wcs"],
                                         RA_DEFL, DEC_DEFL,
                                         box_arcsec=3.0, pix_scale=b["pix_scale"])
        print(f"  centroid pix=({cx:.2f}, {cy:.2f}) [{method}]")

        # Box of ~4" radius (≈ 1.7 R_e on each side) — tight enough that the
        # deflector dominates the fit, but large enough to constrain the wings.
        box_arcsec = 4.0
        fit, off, stats = fit_sersic2d(b["img"], b["mask"], (cx, cy),
                                       r_eff_init_arcsec=R_E_INIT,
                                       box_arcsec=box_arcsec,
                                       pix_scale=b["pix_scale"])
        n = float(fit.n.value)
        r_eff_pix = float(fit.r_eff.value)
        r_eff_arcsec = r_eff_pix * b["pix_scale"]
        amp = float(fit.amplitude.value)
        ellip = float(fit.ellip.value)
        theta_deg = float(np.degrees(fit.theta.value))
        print(f"  sky_bkg={stats['sky_bkg']:.4e}  peak(sky-sub)={stats['peak']:.4e}  amp_init={stats['amp_init']:.4e}")
        print(f"  Sersic2D: n={n:.3f}  r_eff={r_eff_arcsec:.3f}\"  ellip={ellip:.3f}  "
              f"theta={theta_deg:.1f}°  amp={amp:.4e}")
        print(f"  residual median={stats['res_med_pct']*100:.1f}% / p95={stats['res_p95_pct']*100:.1f}%")

        # Total flux: analytic + numeric cross-check
        F_an = sersic_total_flux_analytic(amp, r_eff_pix, n, ellip)
        F_nm = sersic_total_flux_numeric(fit)
        ratio = F_nm / F_an if F_an > 0 else float("nan")
        print(f"  Total flux: analytic={F_an:.4e}  numeric={F_nm:.4e}  numeric/analytic={ratio:.4f}")

        mag_sersic = total_flux_to_ab(F_an, b)
        mag_sersic_num = total_flux_to_ab(F_nm, b)
        delta_mag = mag_sersic - cfg["ap_mag"]
        print(f"  AB_mag(Sersic, analytic) = {mag_sersic:.4f}")
        print(f"  AB_mag(Sersic, numeric)  = {mag_sersic_num:.4f}")
        print(f"  AB_mag(aperture)         = {cfg['ap_mag']:.4f} ± {cfg['ap_err']:.4f}")
        print(f"  Δmag (sersic − aperture) = {delta_mag:+.4f}")

        results[name] = dict(
            n=n, r_eff_arcsec=r_eff_arcsec, r_eff_pix=r_eff_pix,
            amp=amp, ellip=ellip, theta_deg=theta_deg,
            cx=cx, cy=cy, x_off=off[0], y_off=off[1],
            box_arcsec=box_arcsec, pix_scale=b["pix_scale"],
            unit=b["unit"], zp_ab=b["zp_ab"],
            F_total_analytic=F_an, F_total_numeric=F_nm,
            mag_sersic_analytic=mag_sersic, mag_sersic_numeric=mag_sersic_num,
            mag_aperture=cfg["ap_mag"], err_aperture=cfg["ap_err"],
            delta_mag=delta_mag, pivot_AA=cfg["pivot_AA"],
            res_med_pct=float(stats["res_med_pct"]),
            res_p95_pct=float(stats["res_p95_pct"]),
            sub_data=stats["sub"], sub_mask=stats["sub_mask"],
            sub_model=stats["model"], sub_residual=stats["residual"],
        )

    # Summary table
    print("\n" + "=" * 75)
    print("SUMMARY")
    print("=" * 75)
    print(f"{'band':<8} {'n':>6} {'r_eff':>9} {'ellip':>7} "
          f"{'mag_ap':>9} {'mag_ser':>9} {'Δmag':>8} {'res p95':>8}")
    print("-" * 75)
    for name, r in results.items():
        print(f"{name:<8} {r['n']:>6.2f} {r['r_eff_arcsec']:>8.3f}\" {r['ellip']:>7.3f} "
              f"{r['mag_aperture']:>9.4f} {r['mag_sersic_analytic']:>9.4f} "
              f"{r['delta_mag']:>+8.4f} {r['res_p95_pct']*100:>7.1f}%")
    print()

    # Save results — store fit params + delta_mags + photometry vectors for nb08 to load
    # Keep arrays small (subimages) but exclude full-image arrays
    save = {
        "filter_names": np.array(list(results.keys())),
        "pivot_wavelengths_AA": np.array([r["pivot_AA"] for r in results.values()]),
        "mag_aperture": np.array([r["mag_aperture"] for r in results.values()]),
        "mag_aperture_err": np.array([r["err_aperture"] for r in results.values()]),
        "mag_sersic_analytic": np.array([r["mag_sersic_analytic"] for r in results.values()]),
        "mag_sersic_numeric": np.array([r["mag_sersic_numeric"] for r in results.values()]),
        "delta_mag": np.array([r["delta_mag"] for r in results.values()]),
        "sersic_n": np.array([r["n"] for r in results.values()]),
        "sersic_r_eff_arcsec": np.array([r["r_eff_arcsec"] for r in results.values()]),
        "sersic_ellip": np.array([r["ellip"] for r in results.values()]),
        "sersic_theta_deg": np.array([r["theta_deg"] for r in results.values()]),
        "sersic_amplitude": np.array([r["amp"] for r in results.values()]),
        "F_total_analytic": np.array([r["F_total_analytic"] for r in results.values()]),
        "F_total_numeric": np.array([r["F_total_numeric"] for r in results.values()]),
        "residual_median_pct": np.array([r["res_med_pct"] for r in results.values()]),
        "residual_p95_pct": np.array([r["res_p95_pct"] for r in results.values()]),
        "zp_ab": np.array([r["zp_ab"] for r in results.values()]),
        "pix_scale_arcsec": np.array([r["pix_scale"] for r in results.values()]),
    }
    out_npz = "results/sersic_total_photometry.npz"
    np.savez(out_npz, **save)
    print(f"Saved → {out_npz}")

    # 4-band diagnostic figure: per-band data/mask/model/residual/filled
    fig, axes = plt.subplots(4, 5, figsize=(20, 14))
    for row, (name, r) in enumerate(results.items()):
        sub = r["sub_data"]
        sub_mask = r["sub_mask"]
        model = r["sub_model"]
        residual = r["sub_residual"]
        filled = np.where(sub_mask, model, np.nan_to_num(sub))
        if (~sub_mask).any():
            vmax = float(np.nanpercentile(np.abs(sub[~sub_mask]), 99))
        else:
            vmax = float(np.nanpercentile(np.abs(np.nan_to_num(sub)), 99))
        if vmax <= 0:
            vmax = 1.0

        for col, (a, t, cmap, vrange) in enumerate(zip(
            [sub, sub_mask.astype(float), model, residual, filled],
            ["data", "mask", "Sersic model", "residual", "model-filled"],
            ["magma", "Reds", "magma", "RdBu_r", "magma"],
            [(0, vmax), (0, 1), (0, vmax), (-vmax/4, vmax/4), (0, vmax)],
        )):
            ax = axes[row, col]
            ax.imshow(a, origin="lower", cmap=cmap, vmin=vrange[0], vmax=vrange[1])
            ax.set_title(f"{name} — {t}", fontsize=10)
            if col == 0:
                ax.set_ylabel(
                    f"{name}\nn={r['n']:.2f}, r_eff={r['r_eff_arcsec']:.2f}\"\n"
                    f"Δmag={r['delta_mag']:+.3f}", fontsize=9
                )
            ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("§2 — Sersic2D fits per band (rows: F200LP, F140W, F150W2, F322W2)",
                 fontsize=13)
    plt.tight_layout()
    out_png = os.path.join(FIG_DIR, "nb08_sersic_fits.png")
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved → {out_png}")

    return save, results


if __name__ == "__main__":
    main()
