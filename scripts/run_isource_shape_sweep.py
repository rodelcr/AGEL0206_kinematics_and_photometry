"""§6cum I(r)-SHAPE sweep at R<R_e — fixed F200-raw mask, varies I(r) shape only.

Distinct from `run_isource_sweep.py` (which mixed mask choice and I-source
choice). Here we hold the spaxel selection fixed at the headline F200-raw
arc-mask choice, and only vary the SHAPE of the I-weight map.

Ten I-shapes:
    1. IFU_band            (HEADLINE — 6500-7500 Å observed-band IFU sum)
    2. IFU_wl              (full IFU white-light)
    3. F140W_raw           (HST F140W reprojected onto IFU grid, no HST mask)
    4. F200LP_raw          (HST F200LP reprojected onto IFU grid, no HST mask)
    5. F140W_arcmasked     (F140W with arc zeroed via F200LP-mask reprojected → IFU)
    6. F200LP_arcmasked    (F200LP with own arc-mask zeroed → reprojected to IFU)
    7. F140W_1Dcog         (1D CoG annular-mean from masked F140W → interp1d to spaxels)
    8. F200LP_1Dcog        (1D CoG annular-mean from masked F200LP → interp1d to spaxels)
    9. F140W_Sersic2D      (2D Sersic2D fit to F140W with mask weights → IFU grid)
   10. F200LP_Sersic2D     (2D Sersic2D fit to F200LP with mask weights → IFU grid)

For each shape: §6cum cumulative I-weighted ppxf at R<R_E with
    sel = (r_spax < R_E) & ~arc_spax_mask    # F200 raw mask, headline
    w = I_shape[sel]
    flux = sum(cube[:, sel] * (w / w.sum())[None, :], axis=1)
ppxf bootstrap × 3 SPS × 15 polynomial degrees × N=500.

Reuses:
    scripts.run_isource_sweep.load_setup       — cube, center, mask, R_E
    notebooks/07c §4b make_1d_Imap             — ported here verbatim
    notebooks/07c §4c fit_sersic2d / sersic_to_ifu — ported here verbatim

Cache: results/annular_bootstrap_07c_ishape/{shape}_{sps}.npz
"""

import argparse
import os
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Sersic2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402
from scripts.bootstrap_ppxf_parallel import (  # noqa: E402
    run_bootstrap_single_degree_parallel,
)
from scripts.measure_Re import measure_Re_from_profile  # noqa: E402
from scripts.run_isource_sweep import load_setup  # noqa: E402

CACHE_DIR = "results/annular_bootstrap_07c_ishape"
SPS_LIBS = ["fsps", "emiles", "xsl"]
DEGREES = np.arange(15, 30)
WINDOW = 75
Z_SYSTEMIC = 0.67564
BOOT_SEED_BASE = 42 + 70000  # offset so bootstrap draws are independent


# --------------------------------------------------------------------------
# Build the 8 I-shape maps (all on the IFU spaxel grid, ny × nx)
# --------------------------------------------------------------------------
def build_all_shapes(state):
    cube = state["cube"]
    img_f140 = state["img_f140"]
    mask_f140 = state["mask_f140"]
    wcs_f140 = state["wcs_f140"]
    pix_f140 = state["pix_f140"]
    img_f200 = state["img_f200"]
    mask_f200 = state["mask_f200"]
    wcs_f200 = state["wcs_f200"]
    pix_f200 = state["pix_f200"]
    wcs_ifu = state["wcs_ifu"]
    ny, nx = state["ny"], state["nx"]
    r_spax = state["r_spax"]
    R_E = state["R_E"]
    lam = state["lam"]

    # 1. IFU_band — sum cube over 6500-7500 Å
    band_mask = (lam >= 6500.0) & (lam <= 7500.0)
    I_band = np.clip(np.sum(cube[band_mask, :, :], axis=0).astype(float), 0, None)

    # 2. IFU_wl — sum cube over all wavelengths
    I_wl = np.clip(np.sum(cube, axis=0).astype(float), 0, None)

    # 3, 4. HST raw reprojected (no HST mask applied to image)
    yy_, xx_ = np.mgrid[0:ny, 0:nx]
    ra_g, dec_g = wcs_ifu.pixel_to_world_values(xx_.ravel(), yy_.ravel())

    def reproject(img, wcs_hst):
        xh, yh = wcs_hst.world_to_pixel_values(ra_g, dec_g)
        out = map_coordinates(img, [yh, xh], order=1, mode="constant", cval=0.0).reshape(ny, nx)
        return np.clip(out, 0, None)

    I_F140W_raw = reproject(img_f140, wcs_f140)
    I_F200LP_raw = reproject(img_f200, wcs_f200)

    # 5, 6. Direct arc-masked images (deflector-only) reprojected to IFU
    # F140W: the F140W _mask.fits covers the deflector core (NOT arc-only),
    #        so we use the F200LP mask reprojected onto F140W's pixel grid
    #        as the arc-only mask before zeroing.
    F140W_yy, F140W_xx = np.mgrid[:img_f140.shape[0], :img_f140.shape[1]]
    ra_f140grid, dec_f140grid = wcs_f140.pixel_to_world_values(
        F140W_xx.ravel(), F140W_yy.ravel()
    )
    xh200, yh200 = wcs_f200.world_to_pixel_values(ra_f140grid, dec_f140grid)
    arc_on_f140 = map_coordinates(
        mask_f200.astype(float), [yh200, xh200],
        order=0, mode="constant", cval=0.0,
    ).reshape(img_f140.shape).astype(bool)
    img_f140_arcmasked = np.where(arc_on_f140, 0.0, img_f140)
    I_F140W_arcmasked = reproject(img_f140_arcmasked, wcs_f140)

    # F200LP: native mask IS the arc-only mask, so just zero it out
    img_f200_arcmasked = np.where(mask_f200, 0.0, img_f200)
    I_F200LP_arcmasked = reproject(img_f200_arcmasked, wcs_f200)

    # 7, 8. 1D CoG annular-mean → interp1d to spaxels
    # (matches nb07c §4b make_1d_Imap exactly)
    def cog_profile(img, center, pix_scale, mask=None, r_max=6.0, r_step=0.08):
        xc, yc = center
        ny_, nx_ = img.shape
        yy_a, xx_a = np.mgrid[:ny_, :nx_]
        r = np.hypot(xx_a - xc, yy_a - yc) * pix_scale
        edges = np.arange(0, r_max + r_step, r_step)
        rmid, Imean = [], []
        for j in range(len(edges) - 1):
            ann = (r >= edges[j]) & (r < edges[j+1])
            if mask is not None:
                ann = ann & ~mask
            if ann.sum() == 0:
                continue
            Imean.append(float(np.mean(img[ann])))
            rmid.append(0.5 * (edges[j] + edges[j+1]))
        return np.array(rmid), np.array(Imean)

    # F140W center, F200LP center are not in state — recompute from HST_mean RA/Dec
    RA_C, DEC_C = state["RA_C"], state["DEC_C"]
    xc_f140, yc_f140 = wcs_f140.world_to_pixel_values(RA_C, DEC_C)
    xc_f140, yc_f140 = float(xc_f140), float(yc_f140)
    xc_f200, yc_f200 = wcs_f200.world_to_pixel_values(RA_C, DEC_C)
    xc_f200, yc_f200 = float(xc_f200), float(yc_f200)

    r_f140, I_f140cog = cog_profile(img_f140, (xc_f140, yc_f140), pix_f140, mask=mask_f140)
    r_f200, I_f200cog = cog_profile(img_f200, (xc_f200, yc_f200), pix_f200, mask=mask_f200)

    def make_1d_Imap(r_profile, I_profile, r_spax_map):
        I = np.clip(I_profile, 0, None)
        interp = interp1d(r_profile, I, bounds_error=False,
                          fill_value=(I[0], 0.0), kind="linear")
        return np.clip(interp(r_spax_map), 0, None)

    I_F140W_1Dcog = make_1d_Imap(r_f140, I_f140cog, r_spax)
    I_F200LP_1Dcog = make_1d_Imap(r_f200, I_f200cog, r_spax)

    # 9, 10. 2D Sersic2D fits with mask weights → reprojected to IFU
    def fit_sersic2d(img, mask, center, r_eff_init, box_arcsec, pix_scale):
        xc, yc = center
        half = int(np.ceil(box_arcsec / pix_scale))
        y1 = max(0, int(yc) - half)
        y2 = min(img.shape[0], int(yc) + half + 1)
        x1 = max(0, int(xc) - half)
        x2 = min(img.shape[1], int(xc) + half + 1)
        sub = img[y1:y2, x1:x2]
        sub_mask = mask[y1:y2, x1:x2]
        yy_s, xx_s = np.mgrid[:sub.shape[0], :sub.shape[1]]
        r_eff_pix = r_eff_init / pix_scale
        amp_init = float(np.nanpercentile(sub[~sub_mask], 98)) if (~sub_mask).any() else 1.0
        sersic_init = Sersic2D(
            amplitude=amp_init, r_eff=r_eff_pix, n=2.0,
            x_0=xc - x1, y_0=yc - y1, ellip=0.2, theta=0.0,
            bounds={"n": (0.3, 8.0), "r_eff": (r_eff_pix*0.3, r_eff_pix*3),
                    "ellip": (0.0, 0.95), "amplitude": (1e-4, 1e4)},
        )
        weights = (~sub_mask).astype(float)
        try:
            fit = LevMarLSQFitter()(sersic_init, xx_s, yy_s, np.nan_to_num(sub),
                                     weights=weights, maxiter=300)
        except Exception as e:
            print(f"  Sersic fit failed for shape {center}: {e}")
            fit = sersic_init
        return fit, (x1, y1)

    fit_f140, off_f140 = fit_sersic2d(img_f140, mask_f140, (xc_f140, yc_f140), R_E, 6.0, pix_f140)
    fit_f200, off_f200 = fit_sersic2d(img_f200, mask_f200, (xc_f200, yc_f200), R_E, 6.0, pix_f200)
    print(f"  Sersic F140W: n={fit_f140.n.value:.2f}  r_eff={fit_f140.r_eff.value*pix_f140:.3f}\"  "
          f"ellip={fit_f140.ellip.value:.2f}")
    print(f"  Sersic F200LP: n={fit_f200.n.value:.2f}  r_eff={fit_f200.r_eff.value*pix_f200:.3f}\"  "
          f"ellip={fit_f200.ellip.value:.2f}")

    def sersic_to_ifu(fit, offset_xy, wcs_hst):
        x1, y1 = offset_xy
        xh, yh = wcs_hst.world_to_pixel_values(ra_g, dec_g)
        model = fit(xh - x1, yh - y1).reshape(ny, nx)
        return np.clip(model, 0, None)

    I_F140W_Sersic2D = sersic_to_ifu(fit_f140, off_f140, wcs_f140)
    I_F200LP_Sersic2D = sersic_to_ifu(fit_f200, off_f200, wcs_f200)

    return {
        "IFU_band":          I_band,
        "IFU_wl":            I_wl,
        "F140W_raw":         I_F140W_raw,
        "F200LP_raw":        I_F200LP_raw,
        "F140W_arcmasked":   I_F140W_arcmasked,
        "F200LP_arcmasked":  I_F200LP_arcmasked,
        "F140W_1Dcog":       I_F140W_1Dcog,
        "F200LP_1Dcog":      I_F200LP_1Dcog,
        "F140W_Sersic2D":    I_F140W_Sersic2D,
        "F200LP_Sersic2D":   I_F200LP_Sersic2D,
    }


# --------------------------------------------------------------------------
# §6cum bootstrap (same as nb07c, fixed at F200-raw mask)
# --------------------------------------------------------------------------
def run_one_shape(shape_name, I_map, sps, state, n_bootstrap, n_jobs, seed_offset, log_prefix=""):
    cache_fn = os.path.join(CACHE_DIR, f"{shape_name}_{sps}.npz")
    if os.path.exists(cache_fn):
        d = np.load(cache_fn, allow_pickle=True)
        if "sig_boot" in d.files and d["sig_boot"].shape == (len(DEGREES), n_bootstrap):
            print(f"{log_prefix}[skip] cached {cache_fn}")
            return cache_fn

    cube = state["cube"]
    r_spax = state["r_spax"]
    arc_spax_mask = state["arc_spax_mask"]
    R_E = state["R_E"]
    band_mask = state["band_mask"]
    noise_sky = state["noise_sky"]

    # Headline selection: F200-raw mask
    sel = (r_spax < R_E) & ~arc_spax_mask
    w = I_map[sel]
    if w.sum() <= 0:
        print(f"{log_prefix}FAIL: I-weight sum non-positive for {shape_name}")
        return None
    w_norm = w / w.sum()
    flux = np.sum(cube[:, sel] * w_norm[None, :], axis=1)
    n_kept = int(sel.sum())
    sn = float(np.median(flux[band_mask]) / np.median(noise_sky[band_mask]))
    print(f"{log_prefix}[{shape_name:<18s} | {sps:6s}] N_kept={n_kept:3d} S/N_band={sn:5.1f}  I_sum={I_map[sel].sum():.2e}")

    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise_sky.copy(), state["hdr"], sps_name=sps, z=Z_SYSTEMIC, verbose=False,
    )
    n_pix = len(inputs["galaxy"])
    n_deg = len(DEGREES)
    bf = np.zeros((n_deg, n_pix))
    V_o = np.zeros(n_deg)
    sig_o = np.zeros(n_deg)
    V_b = np.full((n_deg, n_bootstrap), np.nan)
    sig_b = np.full((n_deg, n_bootstrap), np.nan)

    t0 = time.time()
    for i, deg in enumerate(DEGREES):
        pp = ppxf(
            inputs["sps"].templates, inputs["galaxy"], inputs["noise"],
            inputs["velscale"], inputs["start"],
            goodpixels=inputs["goodpixels"], plot=False, moments=2,
            trig=False, degree=int(deg), mdegree=0,
            lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"], quiet=True,
        )
        bf[i] = pp.bestfit
        V_o[i] = pp.sol[0]
        sig_o[i] = pp.sol[1]
        rb = run_bootstrap_single_degree_parallel(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_bootstrap, window=WINDOW,
            seed=BOOT_SEED_BASE + seed_offset + i, n_jobs=n_jobs,
        )
        V_b[i] = rb["V_samples"]
        sig_b[i] = rb["sigma_samples"]

    elapsed = time.time() - t0
    out = dict(
        V_orig=V_o, sig_orig=sig_o, V_boot=V_b, sig_boot=sig_b,
        best_fit=bf, galaxy=inputs["galaxy"], noise=inputs["noise"],
        lam_gal_rest=inputs["lam_gal_rest"], goodpixels=inputs["goodpixels"],
        degrees=np.asarray(DEGREES), r_max=R_E, sps=sps,
        ishape=shape_name, n_kept=n_kept, sn_band=sn, n_bootstrap=n_bootstrap,
    )
    os.makedirs(CACHE_DIR, exist_ok=True)
    np.savez(cache_fn, **out)
    print(f"{log_prefix}  → σ_orig={sig_o.min():.0f}-{sig_o.max():.0f} | sig_boot p50={np.median(sig_b):.1f} | {elapsed:.1f}s | saved {cache_fn}")
    return cache_fn


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="Smoke test: 1 shape × 1 SPS × N=10")
    parser.add_argument("--n-bootstrap", type=int, default=250)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--only-shape", default=None)
    args = parser.parse_args()

    state = load_setup()
    print()
    print("Building 8 I-shape maps...")
    I_shapes = build_all_shapes(state)
    print()
    print("Sanity check (each I-map at R<R_e should have non-zero sum):")
    sel = (state["r_spax"] < state["R_E"]) & ~state["arc_spax_mask"]
    for name, I in I_shapes.items():
        s = I[sel].sum()
        peak = np.unravel_index(np.argmax(I), I.shape)
        print(f"  {name:<18s}  sum_in_aperture={s:.3e}  peak_pix={peak}  "
              f"% non-zero in aperture = {100 * (I[sel] > 0).sum()/sel.sum():.1f}%")
    print()

    if args.smoke:
        run_one_shape("IFU_band", I_shapes["IFU_band"], "fsps", state,
                      n_bootstrap=10, n_jobs=args.n_jobs, seed_offset=0,
                      log_prefix="[smoke] ")
        return

    shapes = list(I_shapes.keys())
    if args.only_shape is not None:
        shapes = [s for s in shapes if s == args.only_shape]

    total = len(shapes) * len(SPS_LIBS)
    counter = 0
    t_total = time.time()
    print(f"=== Running {len(shapes)} I-shapes × {len(SPS_LIBS)} SPS × N={args.n_bootstrap} ===\n")
    for s_idx, shape_name in enumerate(shapes):
        for sps_idx, sps in enumerate(SPS_LIBS):
            counter += 1
            seed_offset = 100 * s_idx + 10 * sps_idx
            prefix = f"[{counter}/{total}] "
            try:
                run_one_shape(shape_name, I_shapes[shape_name], sps, state,
                              n_bootstrap=args.n_bootstrap, n_jobs=args.n_jobs,
                              seed_offset=seed_offset, log_prefix=prefix)
            except Exception as e:
                print(f"{prefix}FAILED on {shape_name}/{sps}: {type(e).__name__}: {e}")
    print(f"\n=== DONE in {(time.time()-t_total)/60:.1f} min ===")


if __name__ == "__main__":
    main()
