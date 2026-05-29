"""Multi-systematic arc-mask photometry study (all 4 bands).

Mask philosophy (per user 2026-05-29):
  - F200LP is the most accurate locator of the BLUE z=1.302 source light → it
    defines WHERE the arc is (reprojected to every band).
  - The deeper/sharper IR bands reveal the source's FURTHER EXTENT that HST misses
    → grow the F200LP arc footprint into each band's 2-component-Sersic-residual
    source pixels that are spatially CONTIGUOUS with the arc (region growing, so we
    cannot re-grab the symmetric core pedestal or unrelated field companions).

Deflector model = 2-component (bulge+disk) Sersic2D — single Sersic under-fits the
bright extended JWST galaxy (median fractional residual +1.0 → over-masks); 2-comp
cuts the galaxy-body RMS ~3.5x (validated). PSF convolution unavailable in ISMgas.

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
from scipy.ndimage import label

warnings.filterwarnings("ignore")
REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.mask_method_comparison import (load_band, reproject_mask, _aperture,
                                            flux_mag)
from scripts.arc_mask_verification import sky_sigma, core_exclude

FIG_DIR = "results/figures"
ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
R_E_INIT = 2.305
K_EXT = 3.0          # source-extension residual threshold (sigma)
CORE_R = 0.4         # never mask inside this radius
ARC_RMAX = 5.0       # don't grow the arc beyond this radius (arcsec)


def fit_2comp(b, mask, box_arcsec=6.0):
    """Bulge+disk Sersic2D, masked pixels excluded. Returns (model_full sky-sub, sky_med)."""
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
    amp = float(np.nanpercentile(sub[~sm], 98)) if (~sm).any() else 1.0
    cb = dict(x_0=(cx0 - 3, cx0 + 3), y_0=(cy0 - 3, cy0 + 3))
    bulge = Sersic2D(amplitude=amp, r_eff=reff * 0.3, n=4.0, x_0=cx0, y_0=cy0, ellip=0.2,
                     theta=0, bounds={"n": (2., 6.), "r_eff": (reff * 0.1, reff * 1.0),
                     "ellip": (0., 0.6), "amplitude": (1e-6, 1e6), **cb})
    disk = Sersic2D(amplitude=amp * 0.3, r_eff=reff * 1.2, n=1.0, x_0=cx0, y_0=cy0, ellip=0.2,
                    theta=0, bounds={"n": (0.5, 2.), "r_eff": (reff * 0.6, reff * 2.5),
                    "ellip": (0., 0.7), "amplitude": (1e-6, 1e6), **cb})
    fit = LevMarLSQFitter()(bulge + disk, xx, yy, sub, weights=w, maxiter=1500)
    yyf, xxf = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    model_full = np.clip(np.asarray(fit(xxf - x1, yyf - y1), dtype=float), 0, None)
    return model_full, sky_med


def grow_extension(M_f200, resid, sky_sig, b):
    """Grow the F200LP arc into contiguous IR-significant positive residual."""
    ps, cx, cy = b["pix_scale"], b["cx"], b["cy"]
    yy, xx = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    r = np.hypot(xx - cx, yy - cy) * ps
    sig = (resid > K_EXT * sky_sig) & (r > CORE_R) & (r < ARC_RMAX)
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
    src200 = (f200["mask"].astype(bool), f200["wcs"])

    per_band, models, skys = {}, {}, {}
    for n in ORDER:
        b = bands[n]
        M_f200 = reproject_mask(*src200, b["wcs"], b["img"].shape)
        model, sky_med = fit_2comp(b, M_f200)
        resid = np.nan_to_num(b["img"]) - sky_med - model
        _, sky_sig = sky_sigma(b["img"], b["cx"], b["cy"], b["pix_scale"])
        M_band = grow_extension(M_f200, resid, sky_sig, b)
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
        ax[1].set_title(f"{n} 2-comp model", fontsize=9)
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
