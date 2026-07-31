"""PSF feasibility — field-star census for an EMPIRICAL PSF (2026-06-17, task #10).

The env (ISMgas) has NO webbpsf / stpsf / grizli / TinyTim, so a synthetic PSF
is off the table without changing the env (forbidden). The only route is an
EMPIRICAL PSF built from field stars (photutils EPSFBuilder). This script asks
the prerequisite question: ARE there usable, isolated, unsaturated point sources
in each band's image, and how many?

For each band:
  - estimate background + noise (sigma-clipped),
  - DAOStarFinder at the band's nominal FWHM,
  - keep star-like detections (sharpness/roundness cuts, S/N range that avoids
    saturation and noise), isolated (no neighbour within an exclusion radius),
    and away from the deflector + image edges,
  - measure each candidate's empirical FWHM (radial 2nd moment / Gaussian-ish),
  - report counts and the candidate table.

Output: results/psf_star_census.npz + a per-band printed summary. No PSF is built
here — this only decides feasibility and which bands can support an empirical PSF.

Run from repo root:  python -m scripts.psf_star_census
"""
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")

from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

from scripts.mask_method_comparison import load_band

# nominal PSF FWHM (arcsec) per band — drives DAOStarFinder + isolation radius.
# HST: WFC3 UVIS ~0.07", IR ~0.13"; JWST NIRCam SW ~0.05" (F150W2), LW ~0.13" (F322W2).
FWHM_AS = {"F200LP": 0.075, "F140W": 0.13, "F150W2": 0.05, "F322W2": 0.13}
BANDS = ["F200LP", "F140W", "F150W2", "F322W2"]


def measure_fwhm(cut, fwhm_pix_guess):
    """Crude empirical FWHM via the 2nd radial moment of a background-subtracted,
    flux-weighted cutout (FWHM = 2.355*sigma for a Gaussian). Robust enough to
    flag a real point source vs a cosmic ray / extended blob."""
    ny, nx = cut.shape
    yy, xx = np.mgrid[:ny, :nx]
    w = np.clip(cut - np.median(cut), 0, None)
    if w.sum() <= 0:
        return np.nan
    cx = (w * xx).sum() / w.sum(); cy = (w * yy).sum() / w.sum()
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    sig = np.sqrt((w * r2).sum() / w.sum() / 2.0)   # 2D: var = 2*sigma^2
    return 2.3548 * sig


def census_band(name):
    b = load_band(name)
    img = np.asarray(b["img"], float)
    ps = b["pix_scale"]
    fwhm_pix = FWHM_AS[name] / ps
    cxd, cyd = b["cx"], b["cy"]

    finite = np.isfinite(img)
    mean, med, std = sigma_clipped_stats(img[finite], sigma=3.0, maxiters=5)
    finder = DAOStarFinder(fwhm=fwhm_pix, threshold=8.0 * std,
                           sharplo=0.3, sharphi=1.2, roundlo=-0.6, roundhi=0.6)
    tbl = finder(np.where(finite, img - med, 0.0))
    if tbl is None or len(tbl) == 0:
        return dict(name=name, ps=ps, fwhm_pix=fwhm_pix, std=std, n_found=0,
                    n_usable=0, cands=[])

    xs = np.array(tbl["xcentroid"]); ys = np.array(tbl["ycentroid"])
    flux = np.array(tbl["flux"]); peak = np.array(tbl["peak"])
    sharp = np.array(tbl["sharpness"]); rnd = np.array(tbl["roundness2"])
    snr = peak / std

    # saturation guard: drop the brightest 0.5% peaks (likely saturated/non-linear)
    sat = np.nanpercentile(peak, 99.5)
    ny, nx = img.shape
    edge = 5 * fwhm_pix
    excl = max(8 * fwhm_pix, 12)   # isolation radius (pix)
    rdefl = np.hypot(xs - cxd, ys - cyd)

    cands = []
    for i in range(len(xs)):
        if not (10 <= snr[i] <= 1e4):           # bright but not saturated
            continue
        if peak[i] >= sat:
            continue
        if xs[i] < edge or xs[i] > nx - edge or ys[i] < edge or ys[i] > ny - edge:
            continue
        if rdefl[i] < 6.0 / ps:                 # >6" from the deflector
            continue
        d = np.hypot(xs - xs[i], ys - ys[i]); d[i] = np.inf
        if d.min() < excl:                      # isolation
            continue
        # measure FWHM on a cutout
        half = int(round(4 * fwhm_pix))
        x0, x1 = int(xs[i]) - half, int(xs[i]) + half + 1
        y0, y1 = int(ys[i]) - half, int(ys[i]) + half + 1
        cut = img[y0:y1, x0:x1] - med
        fw = measure_fwhm(cut, fwhm_pix) * ps   # arcsec
        # accept if measured FWHM within 0.6-1.8x nominal (rejects CRs/blobs)
        if not (0.6 * FWHM_AS[name] <= fw <= 1.8 * FWHM_AS[name]):
            continue
        cands.append(dict(x=float(xs[i]), y=float(ys[i]), snr=float(snr[i]),
                          peak=float(peak[i]), sharp=float(sharp[i]),
                          rnd=float(rnd[i]), fwhm_as=float(fw),
                          sep_defl_as=float(rdefl[i] * ps)))
    cands.sort(key=lambda c: -c["snr"])
    return dict(name=name, ps=ps, fwhm_pix=fwhm_pix, std=float(std),
                n_found=len(xs), n_usable=len(cands), cands=cands)


def main():
    out = {}
    print(f"{'band':8s} {'FOV':>10s} {'found':>6s} {'usable':>6s} {'median FWHM(\")':>14s}  verdict")
    print("-" * 70)
    for n in BANDS:
        r = census_band(n)
        out[f"{n}_n_found"] = r["n_found"]
        out[f"{n}_n_usable"] = r["n_usable"]
        if r["cands"]:
            fwhms = np.array([c["fwhm_as"] for c in r["cands"]])
            out[f"{n}_fwhm_med"] = float(np.median(fwhms))
            out[f"{n}_cand_x"] = np.array([c["x"] for c in r["cands"]])
            out[f"{n}_cand_y"] = np.array([c["y"] for c in r["cands"]])
            out[f"{n}_cand_snr"] = np.array([c["snr"] for c in r["cands"]])
            out[f"{n}_cand_fwhm_as"] = fwhms
            med_fw = np.median(fwhms)
        else:
            out[f"{n}_fwhm_med"] = np.nan; med_fw = np.nan
        b = load_band(n)
        fov = b["img"].shape[0] * r["ps"]
        nu = r["n_usable"]
        verdict = ("EPSF OK (>=5)" if nu >= 5 else
                   "marginal (1-4)" if nu >= 1 else "NO stars")
        print(f"{n:8s} {fov:7.0f}\"   {r['n_found']:6d} {nu:6d} {med_fw:14.3f}  {verdict}")
    np.savez("results/psf_star_census.npz", **out)
    print("\nSaved -> results/psf_star_census.npz")
    print("\nTop candidates per band (x, y, S/N, FWHM\"):")
    for n in BANDS:
        r = census_band(n)
        if r["cands"]:
            print(f"  {n}:")
            for c in r["cands"][:6]:
                print(f"    ({c['x']:7.1f},{c['y']:7.1f})  S/N={c['snr']:7.1f}  "
                      f"FWHM={c['fwhm_as']:.3f}\"  sep_defl={c['sep_defl_as']:.1f}\"")


if __name__ == "__main__":
    main()
