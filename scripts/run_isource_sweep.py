"""§6cum I(r) sweep × {arcmask, nomask} at R<R_e.

Reuses nb06's build_I_weight_map machinery (8 (source, mask_strategy) combos)
and applies the §6cum cumulative I-weighted ppxf recipe from nb07c. Produces
per-config caches in results/annular_bootstrap_07c_isource/.

Pre-flight: this script accepts --smoke (1 config × 1 SPS × N=10 ~30 s) and a
full-sweep mode at N=500.

Usage
-----
    conda activate ISMgas
    cd /Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry
    python scripts/run_isource_sweep.py --smoke
    python scripts/run_isource_sweep.py --n-bootstrap 500
"""

import argparse
import os
import sys
import time
from itertools import product

import numpy as np
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Sersic2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.centroids import centroid_2dg
from scipy.ndimage import gaussian_filter, map_coordinates

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.bootstrap_ppxf import (  # noqa: E402
    setup_ppxf_inputs_from_spectrum,
    NOISE_SLICE,
)
from scripts.bootstrap_ppxf_parallel import (  # noqa: E402
    run_bootstrap_single_degree_parallel,
)
from scripts.measure_Re import measure_Re_from_profile  # noqa: E402


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
IFU_FILE = "Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits"
HST_F140W = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits"
HST_F140W_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits"
HST_F200LP = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits"
HST_F200LP_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits"
CACHE_DIR = "results/annular_bootstrap_07c_isource"

Z_SYSTEMIC = 0.67564
RA_DEFL, DEC_DEFL = 31.55611, -1.23817
SPS_LIBS = ["fsps", "emiles", "xsl"]
DEGREES = np.arange(15, 30)
KCWI_SEEING_FWHM = 1.27   # arcsec
CONTAM_THRESHOLD = 0.15
WINDOW = 75
BOOT_SEED_BASE = 42 + 50000  # offset from §6cum (1000) and §6cum-nomask (9000)

# Configurations: (Isrc, maskstrat) — drawn from nb06 viz_combos
NB06_VIZ_COMBOS = [
    ("IFU_band", "unmasked"),
    ("IFU_band", "15pct_psf"),
    ("IFU_wl", "unmasked"),
    ("IFU_wl", "15pct_psf"),
    ("F140W", "unmasked"),
    ("F140W", "arc_only_ifu"),
    ("F200LP", "unmasked"),
    ("F200LP", "arc_only_ifu"),
]

# selstrat options applied to the spaxel selection mask:
#   'arcmask' = (r_spax < R_e) & ~arc_spax_mask   (F200LP raw mask, nb07c convention)
#   'nomask'  = (r_spax < R_e)                     (admit all spaxels in aperture)
SEL_STRATEGIES = ["arcmask", "nomask"]


def cache_path(isrc, maskstrat, selstrat, sps):
    return os.path.join(
        CACHE_DIR, f"{isrc}_{maskstrat}_{selstrat}_{sps}.npz"
    )


# --------------------------------------------------------------------------
# Setup: load cube, find HST_mean center, build maps
# --------------------------------------------------------------------------
def load_setup():
    print("[setup] Loading IFU cube...")
    with fits.open(IFU_FILE) as h:
        hdr = h[0].header
        cube = h[0].data.astype(float)
    wcs_ifu = WCS(hdr, naxis=2)
    ny, nx = cube.shape[1], cube.shape[2]
    crval, cdelt, crpix = hdr["CRVAL3"], hdr["CD3_3"], hdr.get("CRPIX3", 1.0)
    lam = crval + cdelt * (np.arange(hdr["NAXIS3"]) + 1 - crpix)
    pix_scale_ifu = float(np.abs(proj_plane_pixel_scales(wcs_ifu)[0])) * 3600
    band_mask = (lam >= 6500) & (lam <= 7500)

    print(f"[setup] cube shape={cube.shape}  pix={pix_scale_ifu:.3f}\"/spax")

    # HST images and masks
    print("[setup] Loading HST F140W and F200LP ...")
    with fits.open(HST_F140W) as h:
        img_f140 = h[0].data.astype(float)
        wcs_f140 = WCS(h[0].header)
    mask_f140 = fits.getdata(HST_F140W_MASK).astype(bool)
    pix_scale_f140 = float(np.abs(proj_plane_pixel_scales(wcs_f140)[0])) * 3600

    with fits.open(HST_F200LP) as h:
        img_f200 = h[0].data.astype(float)
        wcs_f200 = WCS(h[0].header)
    mask_f200 = fits.getdata(HST_F200LP_MASK).astype(bool)
    pix_scale_f200 = float(np.abs(proj_plane_pixel_scales(wcs_f200)[0])) * 3600

    # HST_mean center (matches nb07c headline)
    def find_center(img, mask, wcs, ra0, dec0, box_arcsec=3.0, pix_scale=None):
        xc0, yc0 = wcs.world_to_pixel_values(ra0, dec0)
        xc0, yc0 = float(xc0), float(yc0)
        half = int(np.ceil(box_arcsec / pix_scale))
        y1 = max(0, int(yc0) - half)
        y2 = min(img.shape[0], int(yc0) + half + 1)
        x1 = max(0, int(xc0) - half)
        x2 = min(img.shape[1], int(xc0) + half + 1)
        sub = np.nan_to_num(img[y1:y2, x1:x2])
        sub_m = mask[y1:y2, x1:x2]
        cx_rel, cy_rel = centroid_2dg(sub, mask=sub_m)
        return float(x1 + cx_rel), float(y1 + cy_rel)

    xc_f140, yc_f140 = find_center(img_f140, mask_f140, wcs_f140, RA_DEFL, DEC_DEFL, pix_scale=pix_scale_f140)
    xc_f200, yc_f200 = find_center(img_f200, mask_f200, wcs_f200, RA_DEFL, DEC_DEFL, pix_scale=pix_scale_f200)
    ra_f140, dec_f140 = wcs_f140.pixel_to_world_values(xc_f140, yc_f140)
    ra_f200, dec_f200 = wcs_f200.pixel_to_world_values(xc_f200, yc_f200)
    RA_C = (float(ra_f140) + float(ra_f200)) / 2
    DEC_C = (float(dec_f140) + float(dec_f200)) / 2
    cx_ifu, cy_ifu = wcs_ifu.world_to_pixel_values(RA_C, DEC_C)
    cx_ifu, cy_ifu = float(cx_ifu), float(cy_ifu)
    print(f"[setup] HST_mean center → IFU sub-pixel ({cx_ifu:.3f}, {cy_ifu:.3f})")

    # r_spax from HST_mean
    yy, xx = np.mgrid[:ny, :nx]
    ra_s, dec_s = wcs_ifu.pixel_to_world_values(xx.ravel(), yy.ravel())
    dra = (ra_s.reshape(ny, nx) - RA_C) * np.cos(np.radians(DEC_C)) * 3600
    ddec = (dec_s.reshape(ny, nx) - DEC_C) * 3600
    r_spax = np.sqrt(dra**2 + ddec**2)

    # ifu_band map
    ifu_band = np.sum(cube[band_mask, :, :], axis=0)

    # arc_spax_mask (raw F200LP reprojected, order=0 — nb07c convention)
    yy_, xx_ = np.mgrid[:ny, :nx]
    ra_g, dec_g = wcs_ifu.pixel_to_world_values(xx_.ravel(), yy_.ravel())
    xh, yh = wcs_f200.world_to_pixel_values(ra_g, dec_g)
    arc_spax_mask = map_coordinates(
        mask_f200.astype(float), [yh, xh], order=0, mode="constant", cval=0.0
    ).reshape(ny, nx).astype(bool)
    print(
        f"[setup] arc_spax_mask (F200 raw): {arc_spax_mask.sum()} / {ny*nx} "
        f"({100*arc_spax_mask.mean():.2f}%)"
    )

    # spaxel_masked (PSF-aware contamination, nb06 convention) — uses F140W mask
    psf_sigma_hst_pix = (KCWI_SEEING_FWHM / pix_scale_f140) / 2.355
    f140_mask_conv = gaussian_filter(mask_f140.astype(float), sigma=psf_sigma_hst_pix)
    xh140, yh140 = wcs_f140.world_to_pixel_values(ra_g, dec_g)
    spaxel_contam = map_coordinates(
        f140_mask_conv, [yh140, xh140], order=1, mode="constant", cval=0.0
    ).reshape(ny, nx)
    spaxel_masked = spaxel_contam > CONTAM_THRESHOLD
    print(
        f"[setup] spaxel_masked (PSF-aware F140W): {spaxel_masked.sum()} / {ny*nx} "
        f"({100*spaxel_masked.mean():.2f}%)"
    )

    # R_e — mean of F140W and F200LP masked CoGs (paper headline)
    KPC_PER_ARCSEC = 7.04

    def cog(img, center, pix_scale, mask=None, r_max=6.0, r_step=0.08):
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
        rmid = np.array(rmid)
        Imean = np.array(Imean)
        Re_, _, _ = measure_Re_from_profile(rmid, Imean)
        return rmid, Imean, Re_

    r_f140, I_f140, Re_f140 = cog(img_f140, (xc_f140, yc_f140), pix_scale_f140, mask=mask_f140)
    r_f200, I_f200, Re_f200 = cog(img_f200, (xc_f200, yc_f200), pix_scale_f200, mask=mask_f200)
    R_E = 0.5 * (Re_f140 + Re_f200)
    print(f"[setup] R_e = mean({Re_f140:.3f}, {Re_f200:.3f}) = {R_E:.3f}\"")

    # Sersic2D fits (will be lazy-loaded inside build_I_weight_map; precompute here)
    def fit_sersic2d(img, mask, center, r_eff_init, box_arcsec, pix_scale):
        xc, yc = center
        half = int(np.ceil(box_arcsec / pix_scale))
        y1 = max(0, int(yc) - half)
        y2 = min(img.shape[0], int(yc) + half + 1)
        x1 = max(0, int(xc) - half)
        x2 = min(img.shape[1], int(xc) + half + 1)
        sub = img[y1:y2, x1:x2]
        sub_mask = mask[y1:y2, x1:x2]
        yy_, xx_ = np.mgrid[:sub.shape[0], :sub.shape[1]]
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
            fit = LevMarLSQFitter()(sersic_init, xx_, yy_, np.nan_to_num(sub),
                                     weights=weights, maxiter=300)
        except Exception as e:
            print(f"  Sersic fit failed: {e}")
            fit = sersic_init
        return fit, (x1, y1)

    # Noise spectrum on native grid
    noise_sky = np.std(cube[:, NOISE_SLICE[0], NOISE_SLICE[1]], axis=(1, 2))

    state = dict(
        hdr=hdr, cube=cube, lam=lam, band_mask=band_mask, ifu_band=ifu_band,
        wcs_ifu=wcs_ifu, ny=ny, nx=nx,
        img_f140=img_f140, mask_f140=mask_f140, wcs_f140=wcs_f140, pix_f140=pix_scale_f140,
        img_f200=img_f200, mask_f200=mask_f200, wcs_f200=wcs_f200, pix_f200=pix_scale_f200,
        cx_ifu=cx_ifu, cy_ifu=cy_ifu, RA_C=RA_C, DEC_C=DEC_C,
        r_spax=r_spax, arc_spax_mask=arc_spax_mask, spaxel_masked=spaxel_masked,
        R_E=R_E, noise_sky=noise_sky,
    )
    return state


# --------------------------------------------------------------------------
# build_I_weight_map (verbatim from nb06 cell 31, with state passed in)
# --------------------------------------------------------------------------
def build_I_weight_map(source, mask_strategy, state):
    """Return a (ny, nx) intensity map on the IFU spaxel grid."""
    cube = state["cube"]
    spaxel_masked = state["spaxel_masked"]
    lam_native = state["lam"]
    ny, nx = state["ny"], state["nx"]
    wcs_ifu_2d = state["wcs_ifu"]

    if source == "IFU_wl":
        I = np.sum(cube, axis=0).astype(float)
        if mask_strategy == "15pct_psf":
            I = np.where(spaxel_masked, 0.0, I)
        elif mask_strategy != "unmasked":
            raise ValueError(f"IFU_wl mask_strategy must be unmasked or 15pct_psf")
        return np.clip(I, 0, None)

    if source == "IFU_band":
        band_mask = (lam_native >= 6500.0) & (lam_native <= 7500.0)
        I = np.sum(cube[band_mask, :, :], axis=0).astype(float)
        if mask_strategy == "15pct_psf":
            I = np.where(spaxel_masked, 0.0, I)
        elif mask_strategy != "unmasked":
            raise ValueError(f"IFU_band mask_strategy must be unmasked or 15pct_psf")
        return np.clip(I, 0, None)

    if source == "F140W":
        img = state["img_f140"].copy()
        wcs_img = state["wcs_f140"]
        msk = state["mask_f140"]
    elif source == "F200LP":
        img = state["img_f200"].copy()
        wcs_img = state["wcs_f200"]
        msk = state["mask_f200"]
    else:
        raise ValueError(f"unknown source {source!r}")

    if mask_strategy == "hst_mask_excl":
        img = np.where(msk, 0.0, img)
    elif mask_strategy not in ("unmasked", "arc_only_ifu"):
        raise ValueError(f"unknown mask_strategy {mask_strategy!r}")

    # Reproject onto IFU grid
    yy_, xx_ = np.mgrid[0:ny, 0:nx]
    ra_, dec_ = wcs_ifu_2d.pixel_to_world_values(xx_.ravel(), yy_.ravel())
    x_im, y_im = wcs_img.world_to_pixel_values(ra_, dec_)
    I_ifu = map_coordinates(img, [y_im, x_im], order=1, mode="constant", cval=0.0).reshape(ny, nx)
    I_ifu = np.clip(I_ifu, 0, None)

    if mask_strategy == "arc_only_ifu":
        I_ifu = np.where(spaxel_masked, 0.0, I_ifu)

    return I_ifu


# --------------------------------------------------------------------------
# §6cum bootstrap for one (Isrc, maskstrat, selstrat, sps)
# --------------------------------------------------------------------------
def build_aperture_spectrum(I_map, state, selstrat):
    """Build the I-weighted aperture spectrum at R<R_e for the given selstrat."""
    cube = state["cube"]
    r_spax = state["r_spax"]
    arc_spax_mask = state["arc_spax_mask"]
    R_E = state["R_E"]
    band_mask = state["band_mask"]
    noise_sky = state["noise_sky"]

    if selstrat == "arcmask":
        sel = (r_spax < R_E) & ~arc_spax_mask
    elif selstrat == "nomask":
        sel = (r_spax < R_E)
    else:
        raise ValueError(f"unknown selstrat {selstrat!r}")

    n_kept = int(sel.sum())
    if n_kept < 1:
        raise RuntimeError("no spaxels in aperture")

    w = I_map[sel]
    w_sum = float(w.sum())
    if w_sum <= 0:
        raise RuntimeError(f"I-weight sum non-positive ({w_sum})")
    w_norm = w / w_sum
    flux = np.sum(cube[:, sel] * w_norm[None, :], axis=1)
    sn = float(np.median(flux[band_mask]) / np.median(noise_sky[band_mask]))
    return flux, noise_sky.copy(), n_kept, sn


def run_one_config(isrc, maskstrat, selstrat, sps, state, n_bootstrap, n_jobs, seed_offset, log_prefix=""):
    cache_fn = cache_path(isrc, maskstrat, selstrat, sps)
    if os.path.exists(cache_fn):
        d = np.load(cache_fn, allow_pickle=True)
        if "sig_boot" in d.files and d["sig_boot"].shape == (len(DEGREES), n_bootstrap):
            print(f"{log_prefix}[skip] cached {cache_fn}")
            return cache_fn, dict(d)

    I_map = build_I_weight_map(isrc, maskstrat, state)
    flux, noise, n_kept, sn = build_aperture_spectrum(I_map, state, selstrat)
    print(f"{log_prefix}[{isrc:8s} | {maskstrat:13s} | {selstrat:7s} | {sps:6s}] N_kept={n_kept:3d} S/N_band={sn:5.1f}")

    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, state["hdr"], sps_name=sps, z=Z_SYSTEMIC, verbose=False
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
        degrees=np.asarray(DEGREES), r_max=state["R_E"], sps=sps,
        isource=isrc, maskstrat=maskstrat, selstrat=selstrat,
        n_kept=n_kept, sn_band=sn, n_bootstrap=n_bootstrap,
    )
    os.makedirs(CACHE_DIR, exist_ok=True)
    np.savez(cache_fn, **out)
    print(f"{log_prefix}  → σ_orig={sig_o.min():.0f}-{sig_o.max():.0f} | sig_boot p50={np.median(sig_b):.1f} | {elapsed:.1f}s | saved {cache_fn}")
    return cache_fn, out


# --------------------------------------------------------------------------
# Build the matrix of configs to actually run
# --------------------------------------------------------------------------
def build_config_matrix():
    """Return list of (isrc, maskstrat, selstrat) tuples to run.

    Skip pairs where the I-map already zero-masks the arc AND selstrat=arcmask
    is exactly redundant (these would only differ if the I-map zeros differ
    from arc_spax_mask, which they generally do — keep them).
    """
    configs = []
    for isrc, maskstrat in NB06_VIZ_COMBOS:
        for selstrat in SEL_STRATEGIES:
            configs.append((isrc, maskstrat, selstrat))
    return configs


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="Smoke test: 1 config × 1 SPS × N=10")
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--only-isrc", default=None,
                        help="Restrict to a single I-source (e.g. IFU_band)")
    args = parser.parse_args()

    state = load_setup()

    if args.smoke:
        # 1 config × 1 SPS × N=10
        cfg = ("IFU_band", "unmasked", "arcmask")
        sps = "fsps"
        run_one_config(*cfg, sps=sps, state=state, n_bootstrap=10, n_jobs=args.n_jobs,
                       seed_offset=0, log_prefix="[smoke] ")
        return

    configs = build_config_matrix()
    if args.only_isrc is not None:
        configs = [c for c in configs if c[0] == args.only_isrc]

    print(f"\n=== Running {len(configs)} configs × {len(SPS_LIBS)} SPS × N={args.n_bootstrap} ===\n")

    total = len(configs) * len(SPS_LIBS)
    counter = 0
    t_total = time.time()
    for cfg_idx, (isrc, maskstrat, selstrat) in enumerate(configs):
        for s_idx, sps in enumerate(SPS_LIBS):
            counter += 1
            seed_offset = 100 * cfg_idx + 10 * s_idx
            prefix = f"[{counter}/{total}] "
            try:
                run_one_config(
                    isrc, maskstrat, selstrat, sps, state,
                    n_bootstrap=args.n_bootstrap, n_jobs=args.n_jobs,
                    seed_offset=seed_offset, log_prefix=prefix,
                )
            except Exception as e:
                print(f"{prefix}FAILED on {(isrc, maskstrat, selstrat, sps)}: {type(e).__name__}: {e}")
    print(f"\n=== DONE in {(time.time()-t_total)/60:.1f} min ===")


if __name__ == "__main__":
    main()
