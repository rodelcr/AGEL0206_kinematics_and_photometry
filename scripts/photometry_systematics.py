"""Multi-systematic arc-mask photometry study (all 4 bands).

Mask philosophy (per user 2026-05-29):
  - F200LP is the most accurate locator of the BLUE z=1.302 source light → it
    defines WHERE the arc is (reprojected to every band).
  - The deeper/sharper IR bands reveal the source's FURTHER EXTENT that HST misses
    → grow the F200LP arc footprint into each band's single-Sersic-residual
    source pixels that are spatially CONTIGUOUS with the arc (region growing, so we
    cannot re-grab the symmetric core pedestal or unrelated field companions).

Deflector model = single-component Sersic2D (per user 2026-06-10) — the deflector is
an elliptical, so a single Sersic (free n, init n=4) is the physically-motivated
light model. NOTE: a single Sersic can under-fit the bright extended JWST galaxy body
(positive residual → grow_extension may over-mask in F150W2/F322W2); inspect the
residual panels and per-band mask sizes after running. PSF convolution unavailable in
ISMgas.

Two mask definitions compared:
  per_band : F200LP-reproj ∪ (this band's IR source extension)
  global   : union of ALL per-band masks, reprojected onto each band

Two flux estimates per (band, mask):
  raw    : photutils aperture, masked pixels excluded (under-counts under-arc light)
  filled : masked pixels replaced by the 2-comp Sersic model (recovers under-arc light)

Output: results/photometry_systematics.npz (4 mag vectors + masks + per-band params)
        results/figures/photsys_{band}.png

Usage: conda activate ISMgas; python scripts/photometry_systematics.py
"""
import os
import sys
import json
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.modeling.models import Sersic2D
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.ndimage import label, shift as ndshift
from scipy.signal import fftconvolve

warnings.filterwarnings("ignore")
REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.mask_method_comparison import (load_band, reproject_mask, _aperture,
                                            flux_mag)
from scripts.arc_mask_verification import sky_sigma, core_exclude, reproject_intensity

FIG_DIR = "results/figures"
ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
R_E_INIT = 2.305
K_EXT = 3.0          # HST source-extension residual threshold (sigma)
K_EXT_JWST = 2.0     # deeper JWST data: lower threshold to capture diffuse arc emission
CORE_R = 0.4         # never mask inside this radius
ARC_RMAX = 5.0       # don't grow the arc beyond this radius (arcsec)
DCOLOR = 0.5         # HST color gate: IR-extension must be >=this many mag bluer
                     # (F200LP-F140W) than the deflector core.
MORPH_FRAC = 1.0     # JWST morphology guard: a pixel joins only if the source EXCESS
                     # exceeds MORPH_FRAC x the single-Sersic deflector model (i.e. it
                     # is source-dominated / OUTSIDE the Sersic-dominated region). Since
                     # the single Sersic captures 66-87% of the deflector flux, this
                     # rejects the deflector body while admitting the diffuse arc.
JWST_BANDS = ("F150W2", "F322W2")


def register_shift(band_img, reproj_img, cx, cy, ps, box_arcsec=3.0):
    """Integer-pixel (dy, dx) that aligns `reproj_img` to `band_img` near the deflector
    via arc cross-correlation.

    Corrects the sub-arcsec HST<->JWST/IR registration residual that a raw WCS
    reprojection leaves behind (the JWST i2d GWCS differs from the FITS-header WCS that
    astropy reads by ~0.1-0.18"; HST UVIS vs IR also differ by ~0.18"). Measured directly
    on the arc light, so it registers exactly what we mask. Returns (0, 0) if either
    cutout is flat."""
    half = int(np.ceil(box_arcsec / ps))
    y1, y2 = max(0, int(cy - half)), int(cy + half)
    x1, x2 = max(0, int(cx - half)), int(cx + half)
    A = np.nan_to_num(band_img[y1:y2, x1:x2]).astype(float)
    B = np.nan_to_num(reproj_img[y1:y2, x1:x2]).astype(float)
    A = A - np.median(A); B = B - np.median(B)
    if A.std() <= 0 or B.std() <= 0 or A.shape != B.shape:
        return 0, 0
    A /= A.std(); B /= B.std()
    cc = fftconvolve(A, B[::-1, ::-1], mode="same")
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    return int(iy - A.shape[0] // 2), int(ix - A.shape[1] // 2)


def reg_shift_for(src, b):
    """Cross-corr shift (dy,dx) registering a reprojected `src` band onto target band b."""
    rep_img = reproject_intensity(src["img"], src["wcs"], src["pix_scale"],
                                  b["wcs"], b["img"].shape)
    return register_shift(b["img"], rep_img, b["cx"], b["cy"], b["pix_scale"])


def arc_color_on_grid(b, f200, f140, sh200=(0, 0), sh140=(0, 0)):
    """F200LP-F140W AB color on band b's grid (blue-minus-red; arc more negative).

    Reproject both HST bands' surface brightness onto b's WCS, **register each to b's
    image** by its cross-corr shift (sh200/sh140), and form the color. Constant ZP
    offsets are dropped: grow_extension thresholds on the contrast vs the deflector
    core, where the offset cancels."""
    sb200 = reproject_intensity(f200["img"], f200["wcs"], f200["pix_scale"],
                                b["wcs"], b["img"].shape)
    sb140 = reproject_intensity(f140["img"], f140["wcs"], f140["pix_scale"],
                                b["wcs"], b["img"].shape)
    if sh200 != (0, 0):
        sb200 = ndshift(np.nan_to_num(sb200), sh200, order=1, mode="constant", cval=np.nan)
    if sh140 != (0, 0):
        sb140 = ndshift(np.nan_to_num(sb140), sh140, order=1, mode="constant", cval=np.nan)
    good = np.isfinite(sb200) & np.isfinite(sb140) & (sb200 > 0) & (sb140 > 0)
    color = np.full(b["img"].shape, np.nan)
    color[good] = -2.5 * np.log10(sb200[good] / sb140[good])  # arc bluer -> more negative
    return color


def fit_sersic(b, mask, box_arcsec=6.0):
    """Single-component Sersic2D deflector model, masked pixels excluded.

    The deflector is an elliptical, so a single Sersic is the physically-motivated
    light model (per user 2026-06-10). Index free (0.5-8), initialised n=3.

    CRITICAL: Sersic2D's `amplitude` is I(r=r_eff), NOT the central peak. With
    I(0)/I(r_eff) = exp(b_n) (~2150 at n=4), initialising amplitude at the peak
    over-predicts the centre by ~2000x and LevMar collapses the amplitude to ~0
    (the whole galaxy then shows as positive residual). So initialise
    amp ~ peak*exp(-b_n) and BOUND amplitude tightly to the peak, exactly as the
    proven fitter in sersic_total_photometry.fit_sersic2d. Returns
    (model_full sky-sub, sky_med).
    """
    ps, xc, yc = b["pix_scale"], b["cx"], b["cy"]
    half = int(np.ceil(box_arcsec / ps))
    y1, y2 = max(0, int(yc) - half), min(b["img"].shape[0], int(yc) + half + 1)
    x1, x2 = max(0, int(xc) - half), min(b["img"].shape[1], int(xc) + half + 1)
    sub = np.nan_to_num(b["img"][y1:y2, x1:x2].astype(float))
    sm = mask[y1:y2, x1:x2]
    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    cx0, cy0 = xc - x1, yc - y1
    r = np.hypot(xx - cx0, yy - cy0) * ps
    ring = (r > 2 * R_E_INIT) & ~sm
    sky_med = float(np.median(sub[ring])) if ring.any() else 0.0
    sub = sub - sky_med
    w = (~sm).astype(float)
    reff = R_E_INIT / ps
    # robust central peak (3x3 median at centre, sky-subtracted) sets the amplitude scale
    ci, cj = int(np.clip(round(cy0), 0, sub.shape[0] - 1)), int(np.clip(round(cx0), 0, sub.shape[1] - 1))
    peak = float(np.median(sub[max(0, ci - 1):ci + 2, max(0, cj - 1):cj + 2]))
    if peak <= 0:
        peak = max(float(np.nanpercentile(sub[~sm], 99.5)) if (~sm).any() else 1.0, 1e-6)
    amp_init = peak * 0.05    # ~ I(r_eff) for n~3 (exp(-b_n) ~ 0.05)
    cb = dict(x_0=(cx0 - 3, cx0 + 3), y_0=(cy0 - 3, cy0 + 3))
    _bnds = {"n": (0.5, 8.), "r_eff": (reff * 0.1, reff * 2.5),
             "ellip": (0., 0.7), "amplitude": (peak * 1e-3, peak * 1.5), **cb}
    # MULTI-START with moment seeding (SExtractor/GALFIT style): the single-Sérsic χ² has
    # a spurious circular minimum (ellip→0) that a low-ellip init rails into; seed (ellip,
    # theta) from flux-weighted 2nd moments + a boosted start + a circular fallback, keep
    # the lowest weighted-residual fit, so we recover the true deflector shape (b/a~0.80-0.86)
    # — matters for the (1-ellip) total-flux factor. (2026-06-11; was ellip=0.2 single start.)
    from scripts.sersic_total_photometry import moment_seed as _moment_seed
    e_mom, th_mom = _moment_seed(sub, w, cx0, cy0, ps, 1.5 * R_E_INIT)
    starts = ([(e_mom, th_mom)]
              + [(0.4, t) for t in np.linspace(0, np.pi, 6, endpoint=False)]
              + [(0.0, 0.0)])
    fit, _rss = None, np.inf
    for e0, th0 in starts:
        sersic = Sersic2D(amplitude=amp_init, r_eff=reff * 0.6, n=3.0, x_0=cx0, y_0=cy0,
                          ellip=min(e0, 0.69), theta=th0, bounds=_bnds)
        cand = LevMarLSQFitter()(sersic, xx, yy, sub, weights=w, maxiter=1500)
        rss = float(np.sum((w * (sub - cand(xx, yy))) ** 2))
        if rss < _rss:
            fit, _rss = cand, rss
    yyf, xxf = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    model_full = np.clip(np.asarray(fit(xxf - x1, yyf - y1), dtype=float), 0, None)
    return model_full, sky_med


def grow_extension(M_f200, resid, sky_sig, b, color=None, model=None,
                   mode="color", k=K_EXT):
    """Grow the F200LP arc seed into contiguous source-extension pixels.

    Two gating modes for the "is this pixel source, not deflector body?" test:

    - mode="color" (HST F200LP/F140W): the pixel must be **bluer than the deflector
      core** by >=DCOLOR in F200LP-F140W. The blue z=1.302 source vs red passive
      deflector contrast. Robust where HST has good S/N.

    - mode="morph" (deep JWST F150W2/F322W2): the pixel must be **source-dominated** —
      its excess over the single-Sersic deflector model exceeds MORPH_FRAC x the model
      (`resid > MORPH_FRAC * model`), i.e. it is OUTSIDE the Sersic-dominated region.
      No HST color is required, so the deep JWST data can capture the **diffuse** arc
      emission that HST is too shallow to color-confirm. Because the single Sersic
      already accounts for the majority of the deflector flux, requiring the excess to
      beat the model rejects the deflector body (model-dominated) while admitting the
      diffuse source (excess-dominated). Threshold k=K_EXT_JWST (lower, deeper data).

    Both modes additionally require a >k*sky_sig residual, r in (CORE_R, ARC_RMAX), and
    spatial contiguity with the bright arc seed (tangential arc-connection)."""
    ps, cx, cy = b["pix_scale"], b["cx"], b["cy"]
    yy, xx = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    r = np.hypot(xx - cx, yy - cy) * ps
    if mode == "morph":
        gate = resid > (MORPH_FRAC * np.nan_to_num(model))   # source-dominated
    else:
        core_ann = (r > 0.4) & (r < 0.8) & np.isfinite(color)
        defl_color = float(np.nanmedian(color[core_ann])) if core_ann.any() else np.nan
        gate = np.isfinite(color) & (color < (defl_color - DCOLOR))
    sig = (resid > k * sky_sig) & (r > CORE_R) & (r < ARC_RMAX) & gate
    seed = M_f200 & (r < ARC_RMAX)
    lbl, n = label(sig | seed)
    keep_labels = set(np.unique(lbl[seed & (lbl > 0)]))
    grown = np.isin(lbl, list(keep_labels)) if keep_labels else seed.copy()
    # preserve the full F200LP footprint; add only IR-grown pixels within the arc radius
    return M_f200 | (grown & (r < ARC_RMAX))


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    bands = {n: load_band(n) for n in ORDER}
    f200 = bands["F200LP"]
    f140 = bands["F140W"]
    src200 = (f200["mask"].astype(bool), f200["wcs"])

    per_band, models, skys = {}, {}, {}
    for n in ORDER:
        b = bands[n]
        # cross-corr registration: align the reprojected F200LP (and F140W) arc to THIS
        # band, correcting the ~0.1-0.18" HST<->JWST/IR WCS residual before masking.
        sh200 = reg_shift_for(f200, b)
        sh140 = reg_shift_for(f140, b)
        M_f200 = reproject_mask(*src200, b["wcs"], b["img"].shape)
        if sh200 != (0, 0):
            M_f200 = ndshift(M_f200.astype(float), sh200, order=0, mode="constant") > 0.5
        ddx, ddy = sh200[1] * b["pix_scale"], sh200[0] * b["pix_scale"]
        print(f"  {n}: F200LP arc-registration shift = ({ddx:+.3f}, {ddy:+.3f})\"")
        model, sky_med = fit_sersic(b, M_f200)
        resid = np.nan_to_num(b["img"]) - sky_med - model
        _, sky_sig = sky_sigma(b["img"], b["cx"], b["cy"], b["pix_scale"])
        if n in JWST_BANDS:
            # deep JWST: morphology guard (source-dominated) captures diffuse arc
            M_band = grow_extension(M_f200, resid, sky_sig, b, model=model,
                                    mode="morph", k=K_EXT_JWST)
        else:
            # HST: color gate (bluer than deflector core)
            color = arc_color_on_grid(b, f200, f140, sh200, sh140)
            M_band = grow_extension(M_f200, resid, sky_sig, b, color=color,
                                    mode="color", k=K_EXT)
        per_band[n] = M_band
        models[n] = (model, sky_med)
        skys[n] = sky_sig
        ext_px = int((M_band & ~M_f200).sum())
        print(f"{n}: F200-reproj={int(M_f200.sum())}px  +IR-extension={ext_px}px  "
              f"-> per-band={int(M_band.sum())}px")

    # global mask = union of all per-band masks reprojected onto each band
    glob = {}
    for t in ORDER:
        bt = bands[t]
        g = np.zeros(bt["img"].shape, bool)
        for s in ORDER:
            g |= reproject_mask(per_band[s], bands[s]["wcs"], bt["wcs"], bt["img"].shape)
        glob[t] = g
        print(f"{t}: global-union={int(g.sum())}px (vs per-band {int(per_band[t].sum())}px)")

    # photometry: raw + filled, for per-band and global masks
    out = {"filter_names": np.array(ORDER),
           "pivot": np.array([bands[n]["cfg"]["pivot_AA"] for n in ORDER]),
           "headline": np.array([bands[n]["cfg"]["ap_mag"] for n in ORDER])}
    vectors = {k: [] for k in ("perband_raw", "perband_filled", "global_raw", "global_filled")}
    print(f"\n{'band':<8}{'mask':<9}{'raw':>9}{'filled':>9}{'fill':>8}")
    for n in ORDER:
        b = bands[n]
        p = json.load(open(b["cfg"]["params"]))
        aper, ann = _aperture(b, p)
        model, sky_med = models[n]
        for mlabel, M in (("perband", per_band[n]), ("global", glob[n])):
            raw = flux_mag(b, b["img"], aper, ann, M)
            filled_img = np.where(M, model + sky_med, np.nan_to_num(b["img"]))
            filled = flux_mag(b, filled_img, aper, ann, None)
            vectors[f"{mlabel}_raw"].append(raw)
            vectors[f"{mlabel}_filled"].append(filled)
            print(f"{n:<8}{mlabel:<9}{raw:>9.4f}{filled:>9.4f}{filled-raw:>+8.4f}")

        # figure
        ps = b["pix_scale"]; half = int(np.ceil(6.0 / ps))
        sl = (slice(max(0, int(b["cy"] - half)), int(b["cy"] + half)),
              slice(max(0, int(b["cx"] - half)), int(b["cx"] + half)))
        flux = np.nan_to_num(b["img"]) - sky_med
        vmax = float(np.nanpercentile(np.abs(flux[sl]), 99))
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        ax[0].imshow(flux[sl], origin="lower", cmap="magma", vmin=0, vmax=vmax)
        M_f200 = reproject_mask(*src200, b["wcs"], b["img"].shape)
        ax[0].contour(M_f200[sl].astype(float), levels=[.5], colors="cyan", linewidths=.8)
        ax[0].contour(per_band[n][sl].astype(float), levels=[.5], colors="lime", linewidths=.8)
        ax[0].contour(glob[n][sl].astype(float), levels=[.5], colors="red", linewidths=.6)
        ax[0].set_title(f"{n}: cyan=F200reproj lime=per-band(+IR) red=global", fontsize=9)
        ax[1].imshow((model)[sl], origin="lower", cmap="magma", vmin=0, vmax=vmax)
        ax[1].set_title(f"{n} single-Sersic model", fontsize=9)
        rv = 6 * skys[n]
        ax[2].imshow((flux - model)[sl], origin="lower", cmap="RdBu_r", vmin=-rv, vmax=rv)
        ax[2].set_title(f"{n} residual (data-model)", fontsize=9)
        for a in ax:
            a.set_xticks([]); a.set_yticks([])
        plt.tight_layout()
        plt.savefig(f"{FIG_DIR}/photsys_{n}.png", dpi=140, bbox_inches="tight"); plt.close()

    for k, v in vectors.items():
        out[k] = np.array(v)
    for n in ORDER:
        out[f"{n}_perband_mask"] = per_band[n]
        out[f"{n}_global_mask"] = glob[n]
    np.savez("results/photometry_systematics.npz", **out)

    print("\n" + "=" * 60)
    print("MAG VECTORS (F200LP, F140W, F150W2, F322W2)")
    for k in ("perband_raw", "perband_filled", "global_raw", "global_filled"):
        print(f"  {k:<16} {np.round(out[k],4).tolist()}")
    print("Saved → results/photometry_systematics.npz")


if __name__ == "__main__":
    main()
