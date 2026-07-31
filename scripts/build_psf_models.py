"""Build PSF models for the 4 photometry bands (2026-06-17, task #10).

Env has no webbpsf/stpsf/grizli and TinyTim is not compiled, so PSFs are built
WITHOUT a synthetic generator, two ways:

  - F140W (WFC3/IR): use the STScI empirical PSF library PSFSTD_WFC3IR_F140W.fits
    (4x-oversampled, 9 detector positions) — the gold-standard ready-made PSF.
    We take the central detector position and block-average to the native cutout
    pixel scale (0.08"/pix). (Field-star EPSF available as a cross-check: 102 stars
    in the full drz.)
  - F200LP (WFC3/UVIS), F150W2 + F322W2 (NIRCam): build an EMPIRICAL effective PSF
    (photutils EPSFBuilder) from isolated, unsaturated field stars. HST stars come
    from the FULL-FRAME drc/drz (the 25" lens cutout has none); JWST from the i2d
    mosaic. The PSF is a property of the detector+drizzle, so a PSF built on the
    full frame applies to the lens cutout carved from the same product.

Saves results/psf_models/<band>_psf.npz with the normalized PSF (sum=1), its pixel
scale, oversampling, and measured FWHM. These feed the PSF-convolved Sérsic refit
(firms up the M* aperture-correction) and the core-Sérsic test (#9).

Run from repo root:  python -m scripts.build_psf_models
"""
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")

from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils.detection import DAOStarFinder
from photutils.psf import extract_stars, EPSFBuilder

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IFU = os.path.join(REPO, "..", "velocity_dispersion_from_IFU")
PSFSTD_F140W = os.path.join(REPO, "..", "202509_CSWA164_modeling", "data_CSWA164",
                            "PSFSTD_WFC3IR_F140W.fits")
OUTDIR = os.path.join(REPO, "results", "psf_models")

# (full-frame image, native pix scale "/pix, nominal FWHM ", DAO fwhm in pix)
EMP = {
    "F200LP": (os.path.join(IFU, "AGEL020613-011417A_F200LP_WFC3_drc_sci.fits"), 0.05, 0.075, 1.5),
    "F150W2": (os.path.join(IFU, "jw05594-o101_t103_nircam_clear-f150w2_i2d.fits"), 0.0307, 0.05, 1.6),
    "F322W2": (os.path.join(IFU, "jw05594-o101_t103_nircam_clear-f322w2_i2d.fits"), 0.063, 0.13, 2.1),
}


def _radial_fwhm(psf, ps, oversamp):
    """FWHM from the half-max isophote AREA (robust for a winged PSF; the 2nd
    radial moment is NOT — it scales with the cutout size). FWHM = 2*sqrt(A/pi)
    where A = number of pixels above half the peak."""
    p = np.clip(psf, 0, None)
    half = 0.5 * p.max()
    area = float((p >= half).sum())          # EPSF-grid pixels above half-max
    fwhm_grid = 2.0 * np.sqrt(area / np.pi)  # equivalent-circle diameter
    return fwhm_grid * ps / oversamp         # arcsec (oversamp grid -> native -> arcsec)


def _get_jwst_sci(path):
    h = fits.open(path)
    for ext in ("SCI", 1, 0):
        try:
            d = h[ext].data
            if d is not None and d.ndim == 2:
                h.close(); return np.asarray(d, float)
        except Exception:
            continue
    h.close(); raise RuntimeError("no 2D SCI")


def build_empirical(name, oversampling=2, n_stars=25, size=25, smoothing_kernel="quartic"):
    path, ps, fwhm_as, fwhm_pix = EMP[name]
    data = _get_jwst_sci(path) if path.endswith("i2d.fits") else np.asarray(fits.getdata(path), float)
    fin = np.isfinite(data) & (data != 0)
    mean, med, std = sigma_clipped_stats(data[fin], sigma=3.0, maxiters=5)
    finder = DAOStarFinder(fwhm=fwhm_pix, threshold=15 * std,
                           sharplo=0.4, sharphi=1.0, roundlo=-0.4, roundhi=0.4)
    t = finder(np.where(fin, data - med, 0.0))
    x = np.array(t["xcentroid"]); y = np.array(t["ycentroid"])
    pk = np.array(t["peak"]); snr = pk / std
    ny, nx = data.shape
    edge = size; excl = max(2.5 * size, 18)
    sat = np.nanpercentile(pk, 99.5)
    rows = []
    for i in np.argsort(-snr):
        if not (20 <= snr[i] <= 2e4) or pk[i] >= sat:
            continue
        if x[i] < edge or x[i] > nx - edge or y[i] < edge or y[i] > ny - edge:
            continue
        d = np.hypot(x - x[i], y - y[i]); d[i] = np.inf
        if d.min() < excl:
            continue
        rows.append((x[i], y[i]))
        if len(rows) >= n_stars:
            break
    clean = np.where(fin, data - med, 0.0)
    clean = np.nan_to_num(clean, nan=0.0, posinf=0.0, neginf=0.0)
    nd = NDData(data=clean)
    # keep only stars whose full cutout is finite & non-degenerate (drops drc gap-edge stars)
    good_rows = []
    h = size // 2
    for (xx0, yy0) in rows:
        xi, yi = int(round(xx0)), int(round(yy0))
        sub = clean[yi - h:yi + h + 1, xi - h:xi + h + 1]
        if sub.shape == (size, size) and np.all(np.isfinite(sub)) and sub.max() > 0:
            good_rows.append((xx0, yy0))
    stars_tbl = Table(rows=good_rows, names=("x", "y"))
    stars = extract_stars(nd, stars_tbl, size=size)
    builder = EPSFBuilder(oversampling=oversampling, maxiters=12,
                          progress_bar=False, recentering_maxiters=10,
                          smoothing_kernel=smoothing_kernel)
    epsf, fitted = builder(stars)
    psf = np.asarray(epsf.data, float)
    psf = np.clip(psf, 0, None); psf /= psf.sum()
    fw = _radial_fwhm(psf, ps, oversampling)
    return dict(psf=psf, pix_scale=ps, oversampling=oversampling, fwhm_as=fw,
                n_stars_used=len(good_rows), source="empirical_EPSF_field_stars",
                fwhm_nominal_as=fwhm_as, image=os.path.basename(path))


def build_f140w_psfstd():
    cube = np.asarray(fits.getdata(PSFSTD_F140W), float)   # (9,101,101), 4x oversampled
    psf = cube[4]                                          # central detector position
    psf = np.clip(psf, 0, None); psf /= psf.sum()
    ps_native = 0.08                                       # F140W cutout drizzle scale
    fw = _radial_fwhm(psf, ps_native, oversamp=4)
    return dict(psf=psf, pix_scale=ps_native, oversampling=4, fwhm_as=fw,
                n_stars_used=0, source="STScI_PSFSTD_WFC3IR_F140W_centralpos",
                fwhm_nominal_as=0.13, image=os.path.basename(PSFSTD_F140W))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    results = {}
    print("Building F140W from STScI PSFSTD library...")
    results["F140W"] = build_f140w_psfstd()
    # per-band builder settings; F200LP (UVIS) is critically undersampled (FWHM~1.5pix)
    # → quartic smoothing produces NaNs (use quadratic). oversampling=2 recovers the
    # correct core FWHM (0.085" ≈ UVIS nominal); ov=4 over-sharpens it, ov=1 is too coarse.
    EMP_OPTS = {"F200LP": dict(oversampling=2, smoothing_kernel="quadratic"),
                "F150W2": dict(oversampling=2, smoothing_kernel="quartic"),
                "F322W2": dict(oversampling=2, smoothing_kernel="quartic")}
    for n in EMP:
        print(f"Building {n} empirical EPSF from field stars...")
        try:
            results[n] = build_empirical(n, **EMP_OPTS[n])
        except Exception as e:
            print(f"  {n}: FAILED ({type(e).__name__}: {e})")
            results[n] = None

    print(f"\n{'band':8s} {'source':40s} {'Nstars':>6s} {'PSF FWHM(\")':>11s} {'(nominal)':>10s}")
    print("-" * 80)
    for n in ["F200LP", "F140W", "F150W2", "F322W2"]:
        r = results.get(n)
        if r is None:
            print(f"{n:8s} FAILED"); continue
        np.savez(os.path.join(OUTDIR, f"{n}_psf.npz"),
                 psf=r["psf"], pix_scale=r["pix_scale"], oversampling=r["oversampling"],
                 fwhm_as=r["fwhm_as"], n_stars_used=r["n_stars_used"],
                 source=r["source"], image=r["image"])
        print(f"{n:8s} {r['source']:40s} {r['n_stars_used']:6d} "
              f"{r['fwhm_as']:11.3f} {r['fwhm_nominal_as']:10.3f}")
    print(f"\nSaved PSF kernels -> {OUTDIR}/<band>_psf.npz")


if __name__ == "__main__":
    main()
