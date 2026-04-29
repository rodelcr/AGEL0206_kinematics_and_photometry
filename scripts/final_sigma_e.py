"""Final paper σ_e production pipeline (notebook 09).

Best-practice settings baked in:
  * Frame-aware ppxf: FSPS galaxy in vacuum, EMILES + XSL galaxy in air.
    Removes the −90 km/s FSPS V_sys offset diagnosed Apr 2026 (vac/air mismatch).
  * §6cum cumulative I-weighted aperture ppxf at three apertures:
      R < R_e/8, R < R_e/2, R < R_e   (Kormendy & Ho 2013 / Greene+2020 σ_e).
  * I-weight = IFU 6500–7500 Å white-light band (deflector-dominant).
  * F200LP arc mask reprojected onto the IFU grid (nearest-neighbor) and
    hard-applied to drop arc spaxels.
  * HST-mean center: F140W + F200LP centroid_2dg averaged in world coordinates,
    then propagated to IFU sub-pixel via WCS.
  * R_e = 0.5 × (R_e_F140W + R_e_F200LP) from masked curve-of-growth = 2.305".
  * Per-SPS bootstrap (FSPS/EMILES/XSL) over additive polynomial degrees
    15–29 → SPS-pooled posterior at each aperture.
  * z_systemic = 0.67564 (line-fit, notebook 04).

Output: results/final_sigma_e_paper.npz  + per-(aperture, SPS) caches in
        results/final_sigma_e_paper/

CLI:
    python scripts/final_sigma_e.py --n_bootstrap 50          # smoke (~6 min)
    python scripts/final_sigma_e.py --n_bootstrap 500         # production (~1 h)
    python scripts/final_sigma_e.py --n_bootstrap 50 --force  # rebuild caches
"""

import argparse
import os
import sys
from pathlib import Path
from time import perf_counter as clock

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.cosmology import FlatLambdaCDM
from photutils.centroids import centroid_2dg
from scipy.ndimage import map_coordinates
from ppxf.ppxf import ppxf

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from scripts.bootstrap_ppxf import (   # noqa: E402
    setup_ppxf_inputs_from_spectrum,
    run_bootstrap_single_degree,
    SPS_NATIVE_FRAME,
    NOISE_SLICE,
    C_KMS,
)
from scripts.measure_Re import measure_Re_from_profile  # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
IFU_FILE = "Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits"
HST_F140W = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits"
HST_F140W_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits"
HST_F200LP = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits"
HST_F200LP_MASK = "../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits"

CACHE_DIR = REPO / "results" / "final_sigma_e_paper"
RESULTS_NPZ = REPO / "results" / "final_sigma_e_paper.npz"

Z_SYSTEMIC = 0.67564
SPS_LIBS = ("fsps", "emiles", "xsl")
DEGREES = np.arange(15, 30)
BOOT_SEED = 42
WINDOW = 75
RA_DEFL, DEC_DEFL = 31.55611, -1.23817
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
KPC_PER_ARCSEC = COSMO.kpc_proper_per_arcmin(Z_SYSTEMIC).value / 60
# R<R_e/8 dropped (seeing-limited: 3 spaxels at 0.29" is below FWHM/2 = 0.64").
APERTURE_FRACS = (1.0 / 2, 1.0)
APERTURE_LABELS = ("Re_2", "Re")


# ─────────────────────────────────────────────────────────────────────────────
# Setup helpers
# ─────────────────────────────────────────────────────────────────────────────
def find_center(img, mask, wcs, ra0, dec0, box_arcsec=3.0, pix_scale=None):
    xc0, yc0 = wcs.world_to_pixel_values(ra0, dec0)
    xc0, yc0 = float(xc0), float(yc0)
    half = int(np.ceil(box_arcsec / pix_scale))
    y1 = max(0, int(yc0) - half); y2 = min(img.shape[0], int(yc0) + half + 1)
    x1 = max(0, int(xc0) - half); x2 = min(img.shape[1], int(xc0) + half + 1)
    sub = np.nan_to_num(img[y1:y2, x1:x2]); sub_m = mask[y1:y2, x1:x2]
    cx_rel, cy_rel = centroid_2dg(sub, mask=sub_m)
    if not (np.isfinite(cx_rel) and np.isfinite(cy_rel)):
        return xc0, yc0
    return float(x1 + cx_rel), float(y1 + cy_rel)


def cahk_g_line_depth_map(cube, lam_obs, z):
    """Per-spaxel total Ca H/K + G-band absorption *flux* deficit (stellar I-map).

    For each absorption feature the continuum is estimated as the median of
    two side windows; the metric is
        sum_pix (continuum - flux)   over the line window
    summed across all three lines and clipped at zero. Units are flux ×
    pixels — i.e., the "missing flux" filled in by the absorption. Unlike
    a normalised line depth (cont-flux)/cont, this metric scales with
    continuum brightness and so falls to ~0 at large radii where the
    deflector does not contribute, giving a stable curve-of-growth.

    Rest-frame line / continuum windows (Å):
        Ca K   line 3925-3942   cont 3895-3915, 3955-3975
        Ca H   line 3960-3976   cont 3940-3955, 3985-4005
        G-band line 4297-4313   cont 4275-4290, 4320-4340
    """
    lam_rest = lam_obs / (1.0 + z)
    bands = (
        {"line": (3925, 3942), "cont": [(3895, 3915), (3955, 3975)]},
        {"line": (3960, 3976), "cont": [(3940, 3955), (3985, 4005)]},
        {"line": (4297, 4313), "cont": [(4275, 4290), (4320, 4340)]},
    )
    n_lam, ny, nx = cube.shape
    out = np.zeros((ny, nx), dtype=float)
    for w in bands:
        line_mask = (lam_rest >= w["line"][0]) & (lam_rest <= w["line"][1])
        cont_mask = np.zeros_like(line_mask)
        for c0, c1 in w["cont"]:
            cont_mask |= (lam_rest >= c0) & (lam_rest <= c1)
        if line_mask.sum() == 0 or cont_mask.sum() == 0:
            continue
        cont = np.median(cube[cont_mask], axis=0)
        line_mean = np.mean(cube[line_mask], axis=0)
        # integrated absorption flux deficit, in flux × pixel-equivalent units
        out += np.nan_to_num(cont - line_mean) * float(line_mask.sum())
    return np.clip(out, 0.0, None)


def curve_of_growth(image, center, pix_scale, mask=None, r_max=6.0, r_step=0.08):
    xc, yc = center
    ny_, nx_ = image.shape
    yy, xx = np.mgrid[:ny_, :nx_]
    r = np.hypot(xx - xc, yy - yc) * pix_scale
    edges = np.arange(0, r_max + r_step, r_step)
    r_mid, I_mean = [], []
    for j in range(len(edges) - 1):
        ann = (r >= edges[j]) & (r < edges[j + 1])
        if mask is not None:
            ann = ann & ~mask
        if ann.sum() == 0:
            continue
        I_mean.append(float(np.mean(image[ann])))
        r_mid.append(0.5 * (edges[j] + edges[j + 1]))
    r_mid = np.asarray(r_mid); I_mean = np.asarray(I_mean)
    Re, _, _ = measure_Re_from_profile(r_mid, I_mean)
    return Re


def load_setup():
    """Load cube + HST images + masks, compute center, R_e, r_spax map, arc mask."""
    print("─" * 70)
    print("§1  Load IFU cube + HST images + arc masks")
    print("─" * 70)
    with fits.open(IFU_FILE) as h:
        hdr = h[0].header
        cube = h[0].data.astype(float)
    wcs_ifu = WCS(hdr, naxis=2)
    ny, nx = cube.shape[1], cube.shape[2]
    crval = hdr["CRVAL3"]; cdelt = hdr["CD3_3"]; crpix = hdr.get("CRPIX3", 1.0)
    lam = crval + cdelt * (np.arange(hdr["NAXIS3"]) + 1 - crpix)
    pix_scale_ifu = float(np.abs(proj_plane_pixel_scales(wcs_ifu)[0])) * 3600
    band_mask = (lam >= 6500) & (lam <= 7500)
    ifu_band = np.sum(cube[band_mask, :, :], axis=0)
    print(f"  IFU cube      shape={cube.shape}  pix={pix_scale_ifu:.3f}\"/spax  "
          f"λ={lam[0]:.0f}-{lam[-1]:.0f} Å")

    with fits.open(HST_F140W) as h:
        img_f140 = h[0].data.astype(float)
    wcs_f140 = WCS(fits.getheader(HST_F140W))
    mask_f140 = fits.getdata(HST_F140W_MASK).astype(bool)
    pix_scale_f140 = float(np.abs(proj_plane_pixel_scales(wcs_f140)[0])) * 3600

    with fits.open(HST_F200LP) as h:
        img_f200 = h[0].data.astype(float)
    wcs_f200 = WCS(fits.getheader(HST_F200LP))
    mask_f200 = fits.getdata(HST_F200LP_MASK).astype(bool)
    pix_scale_f200 = float(np.abs(proj_plane_pixel_scales(wcs_f200)[0])) * 3600
    print(f"  F140W+mask    pix={pix_scale_f140:.3f}\"/pix  mask={mask_f140.sum()} pix")
    print(f"  F200LP+mask   pix={pix_scale_f200:.3f}\"/pix  mask={mask_f200.sum()} pix (arc-tuned)")

    print("\n§2  HST-mean center via centroid_2dg on F140W + F200LP")
    xc_140, yc_140 = find_center(img_f140, mask_f140, wcs_f140, RA_DEFL, DEC_DEFL,
                                  pix_scale=pix_scale_f140)
    xc_200, yc_200 = find_center(img_f200, mask_f200, wcs_f200, RA_DEFL, DEC_DEFL,
                                  pix_scale=pix_scale_f200)
    ra_140, dec_140 = wcs_f140.pixel_to_world_values(xc_140, yc_140)
    ra_200, dec_200 = wcs_f200.pixel_to_world_values(xc_200, yc_200)
    d_centroid = np.hypot(
        (ra_140 - ra_200) * 3600 * np.cos(np.radians((dec_140 + dec_200) / 2)),
        (dec_140 - dec_200) * 3600,
    )
    print(f"  F140W centroid  (RA,Dec) = ({float(ra_140):.5f}, {float(dec_140):.5f})")
    print(f"  F200LP centroid (RA,Dec) = ({float(ra_200):.5f}, {float(dec_200):.5f})")
    print(f"  |Δcentroid| = {d_centroid:.3f}\"  ({'PASS' if d_centroid < 0.2 else 'FAIL'})")

    ra_c = (float(ra_140) + float(ra_200)) / 2
    dec_c = (float(dec_140) + float(dec_200)) / 2
    cx_ifu, cy_ifu = wcs_ifu.world_to_pixel_values(ra_c, dec_c)
    cx_ifu, cy_ifu = float(cx_ifu), float(cy_ifu)
    print(f"  HST-mean → IFU sub-pixel ({cx_ifu:.3f}, {cy_ifu:.3f})")

    print("\n§3  R_e from masked curve-of-growth (F140W + F200LP)")
    Re_140 = curve_of_growth(img_f140, (xc_140, yc_140), pix_scale_f140, mask=mask_f140)
    Re_200 = curve_of_growth(img_f200, (xc_200, yc_200), pix_scale_f200, mask=mask_f200)
    R_E = 0.5 * (Re_140 + Re_200)
    print(f"  R_e(F140W masked)  = {Re_140:.3f}\"")
    print(f"  R_e(F200LP masked) = {Re_200:.3f}\"")
    print(f"  R_e (headline mean) = {R_E:.3f}\"  = {R_E*KPC_PER_ARCSEC:.2f} kpc")

    # R_e from a Ca H+K + G-band absorption-depth I-map (stellar-only). Built on
    # the IFU cube directly, so the centre is the IFU sub-pixel HST-mean point.
    cahk_map = cahk_g_line_depth_map(cube, lam, Z_SYSTEMIC)
    Re_cahk = curve_of_growth(cahk_map, (cx_ifu, cy_ifu), pix_scale_ifu)
    print(f"  R_e(CaHK + G-band depth) = {Re_cahk:.3f}\"  "
          f"({(Re_cahk - R_E):+.3f}\" vs headline mean)")

    print("\n§4  r_spax map + reprojected F200LP arc mask")
    yy_, xx_ = np.mgrid[:ny, :nx]
    ra_s, dec_s = wcs_ifu.pixel_to_world_values(xx_.ravel(), yy_.ravel())
    dra = (ra_s.reshape(ny, nx) - ra_c) * np.cos(np.radians(dec_c)) * 3600
    ddec = (dec_s.reshape(ny, nx) - dec_c) * 3600
    r_spax = np.sqrt(dra**2 + ddec**2)
    xh, yh = wcs_f200.world_to_pixel_values(ra_s, dec_s)
    arc_spax_mask = map_coordinates(
        mask_f200.astype(float), [yh, xh], order=0, mode="constant", cval=0.0
    ).reshape(ny, nx).astype(bool)
    print(f"  arc_spax_mask: {arc_spax_mask.sum()} of {ny*nx} spaxels flagged "
          f"({100*arc_spax_mask.mean():.2f}%)")

    noise_sky = np.std(cube[:, NOISE_SLICE[0], NOISE_SLICE[1]], axis=(1, 2))

    return dict(
        hdr=hdr, cube=cube, lam=lam, ifu_band=ifu_band, noise_sky=noise_sky,
        ny=ny, nx=nx, pix_scale_ifu=pix_scale_ifu, band_mask=band_mask,
        cx_ifu=cx_ifu, cy_ifu=cy_ifu,
        ra_center=ra_c, dec_center=dec_c, d_centroid=float(d_centroid),
        Re_140=Re_140, Re_200=Re_200, R_E=R_E,
        Re_cahk=float(Re_cahk), cahk_map=cahk_map,
        r_spax=r_spax, arc_spax_mask=arc_spax_mask,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Aperture extraction
# ─────────────────────────────────────────────────────────────────────────────
def extract_aperture_spectrum(state, r_max, mask_weight=None, mask_on=None):
    """I-weighted spectrum of all spaxels with r < r_max.

    The `mask_weight` parameter generalizes mask handling:
      0.0 → hard mask (drop arc spaxels via zero weight) — DEFAULT, headline behavior
      1.0 → no mask (use all spaxels at full I-weight) — sensitivity check
      0.5 → soft mask (arc spaxels kept but down-weighted by 0.5) — interpolation track
      arbitrary float in [0, 1] → continuous weighting

    Backward compat: if `mask_on` is given (legacy boolean), translate to
    mask_weight = 0.0 (mask_on=True) or 1.0 (mask_on=False).
    """
    if mask_weight is None:
        mask_weight = 0.0 if (mask_on is None or mask_on) else 1.0

    sel = state["r_spax"] < r_max
    I = np.clip(state["ifu_band"], 0, None)
    I_eff = I.copy()
    I_eff[state["arc_spax_mask"]] *= mask_weight
    w = I_eff[sel]
    n_active = int((w > 0).sum())
    w_norm = w / max(w.sum(), 1e-30)
    flux = np.sum(state["cube"][:, sel] * w_norm[None, :], axis=1)
    sn_band = float(np.median(flux[state["band_mask"]])
                    / np.median(state["noise_sky"][state["band_mask"]]))
    return flux, state["noise_sky"].copy(), n_active, sn_band


# ─────────────────────────────────────────────────────────────────────────────
# Per-(aperture, SPS) ppxf + bootstrap
# ─────────────────────────────────────────────────────────────────────────────
def _mask_suffix(mask_weight):
    """Cache filename suffix for a given mask_weight."""
    if mask_weight == 0.0:
        return ""
    if mask_weight == 1.0:
        return "_nomask"
    s = f"{mask_weight:.2f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"_softmask_w{s}"


def run_aperture_sps(state, label, r_max, sps_name, n_bootstrap, force,
                     mask_on=True, mask_weight=None, fallback_n=None):
    """Run §6cum at one (label, sps, mask) config with bootstrap.

    Cache lookup chain: N={n_bootstrap}, then N={fallback_n} (if specified).
    The actual N used is recorded in the returned dict's `n_bootstrap` key.
    """
    if mask_weight is None:
        mask_weight = 0.0 if mask_on else 1.0
    suffix = _mask_suffix(mask_weight)
    chain = [n_bootstrap]
    if fallback_n is not None and fallback_n != n_bootstrap:
        chain.append(fallback_n)
    if not force:
        for n_try in chain:
            cache = CACHE_DIR / f"{label}_{sps_name}_N{n_try}{suffix}.npz"
            if cache.exists():
                d = dict(np.load(cache, allow_pickle=True))
                for k in list(d.keys()):
                    if d[k].shape == ():
                        d[k] = d[k].item()
                return d
    cache = CACHE_DIR / f"{label}_{sps_name}_N{n_bootstrap}{suffix}.npz"

    flux, noise, n_kept, sn_band = extract_aperture_spectrum(
        state, r_max, mask_weight=mask_weight,
    )
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, state["hdr"], sps_name=sps_name, z=Z_SYSTEMIC,
        verbose=False, frame_galaxy="auto",
    )
    n_pix = len(inputs["galaxy"])
    n_deg = len(DEGREES)
    V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
    V_boot = np.full((n_deg, n_bootstrap), np.nan)
    sig_boot = np.full((n_deg, n_bootstrap), np.nan)
    bf_arr = np.zeros((n_deg, n_pix))
    t0 = clock()
    for i, deg in enumerate(DEGREES):
        pp = ppxf(
            inputs["sps"].templates, inputs["galaxy"], inputs["noise"],
            inputs["velscale"], inputs["start"],
            goodpixels=inputs["goodpixels"], plot=False, moments=2, trig=False,
            degree=int(deg), mdegree=0,
            lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"], quiet=True,
        )
        V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
        bf_arr[i] = pp.bestfit
        try:
            label_idx = APERTURE_LABELS.index(label)
        except ValueError:
            # Ad-hoc labels (e.g. R_e systematic variants) get a stable
            # offset above the APERTURE_LABELS range.
            label_idx = 100 + (abs(hash(label)) % 1000)
        rb = run_bootstrap_single_degree(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_bootstrap, window=WINDOW,
            seed=BOOT_SEED + 1000 * label_idx + i,
        )
        V_boot[i] = rb["V_samples"]
        sig_boot[i] = rb["sigma_samples"]
    elapsed = clock() - t0

    out = dict(
        label=label, sps_name=sps_name, r_max=float(r_max),
        mask_on=bool(mask_weight == 0.0), mask_weight=float(mask_weight),
        n_spax=int(n_kept), sn_band=float(sn_band),
        degrees=np.asarray(DEGREES),
        V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
        V_boot=V_boot, sig_boot=sig_boot,
        frame_galaxy=inputs["frame_galaxy"],
        velscale=float(inputs["velscale"]),
        n_pix=int(n_pix), n_bootstrap=int(n_bootstrap),
        elapsed_s=float(elapsed),
        galaxy=inputs["galaxy"], noise=inputs["noise"],
        lam_gal_rest=inputs["lam_gal_rest"], best_fit=bf_arr,
    )
    np.savez(cache, **out)
    tag = (
        "masked" if mask_weight == 0.0
        else "nomask" if mask_weight == 1.0
        else f"soft_w{mask_weight}"
    )
    print(f"    [{tag:11s}] {label}/{sps_name:6s}  σ={sig_orig.min():.0f}-{sig_orig.max():.0f} "
          f"V={V_orig.mean():+.1f}  frame={inputs['frame_galaxy']:7s}  "
          f"{elapsed:.0f}s  ({n_kept} spax, S/N={sn_band:.1f})")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# SPS pooling + error budget
# ─────────────────────────────────────────────────────────────────────────────
def pool_sps(per_sps):
    """Pool the per-SPS bootstrap distributions (concatenate, all degrees + all SPS)."""
    pooled_sig, pooled_V = [], []
    per_sps_summary = {}
    for sps, d in per_sps.items():
        # Subtract per-SPS V_sys (median across all bootstrap samples) to remove
        # constant frame-residual offsets; this matches the §6cum-aware Gültekin
        # pipeline elsewhere in the project.
        V = d["V_boot"]
        sig = d["sig_boot"]
        V_sys = float(np.nanmedian(V))
        V_offset = V - V_sys
        pooled_V.append(V_offset.ravel())
        pooled_sig.append(sig.ravel())
        per_sps_summary[sps] = dict(
            sigma_p16=float(np.nanpercentile(sig, 16)),
            sigma_p50=float(np.nanpercentile(sig, 50)),
            sigma_p84=float(np.nanpercentile(sig, 84)),
            V_sys=V_sys,
            frame=d["frame_galaxy"],
        )
    sigma_pool = np.concatenate(pooled_sig)
    V_pool = np.concatenate(pooled_V)
    sigma_pool = sigma_pool[np.isfinite(sigma_pool)]
    V_pool = V_pool[np.isfinite(V_pool)]
    return dict(
        sigma_p16=float(np.percentile(sigma_pool, 16)),
        sigma_p50=float(np.percentile(sigma_pool, 50)),
        sigma_p84=float(np.percentile(sigma_pool, 84)),
        V_p16=float(np.percentile(V_pool, 16)),
        V_p50=float(np.percentile(V_pool, 50)),
        V_p84=float(np.percentile(V_pool, 84)),
        sigma_samples=sigma_pool,
        per_sps=per_sps_summary,
    )


def error_budget(pool_Re):
    """Quadrature-combine known systematics with statistical pooled error.

    Stat:        bootstrap pooled posterior 1σ (asymmetric → use mean of ±)
    I-shape:     13 km/s — from nb07c I-shape sweep (analyze_isource_shape_sweep.py)
    Mask:        16 km/s — F200LP raw vs no-mask (results/annular_bootstrap_07c_nomask)
    Frame fix:    5 km/s — σ shift between air and vacuum galaxy (max across SPS)
    Centering:    4 km/s — 5-center sweep (NOTES_centering_investigation_2026-04-27.md)
    """
    stat = 0.5 * ((pool_Re["sigma_p84"] - pool_Re["sigma_p50"])
                  + (pool_Re["sigma_p50"] - pool_Re["sigma_p16"]))
    I_shape = 13.0
    mask = 16.0
    frame = 5.0
    centering = 4.0
    total = float(np.sqrt(stat**2 + I_shape**2 + mask**2 + frame**2 + centering**2))
    return dict(
        stat=float(stat), I_shape=I_shape, mask=mask, frame=frame, centering=centering,
        total=total,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────────────────────
def main(n_bootstrap, force):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n{'═'*70}\nFINAL σ_e PIPELINE — N_BOOTSTRAP={n_bootstrap}, force={force}\n{'═'*70}")
    print(f"  SPS_NATIVE_FRAME = {SPS_NATIVE_FRAME}")
    print(f"  z_systemic = {Z_SYSTEMIC}    DEGREES = {DEGREES[0]}..{DEGREES[-1]}")
    print(f"  SPS_LIBS = {SPS_LIBS}        APERTURES = {APERTURE_LABELS}")
    print(f"  cache → {CACHE_DIR.relative_to(REPO)}")

    state = load_setup()

    print(f"\n§5  Three tracks × {len(APERTURE_LABELS)} apertures × {len(SPS_LIBS)} SPS libraries\n" + "─"*70)
    print(f"      Track A — F200LP arc mask ON   (mask_weight=0.0, headline)")
    print(f"      Track B — arc mask OFF         (mask_weight=1.0, sensitivity)")
    print(f"      Track C — soft mask            (mask_weight=0.5, interpolation)")
    tracks = {}
    track_specs = (
        ("masked",         0.0),
        ("nomask",         1.0),
        ("softmask_w0p5",  0.5),
    )
    for track_label, mask_weight in track_specs:
        print(f"\n  ─── Track: {track_label} (mask_weight={mask_weight}) ───")
        apertures = {}
        for frac, label in zip(APERTURE_FRACS, APERTURE_LABELS):
            r_max = frac * state["R_E"]
            print(f"  {track_label}/{label}:  R_max = {r_max:.3f}\"  ({frac:.3f} R_e)")
            per_sps = {}
            for sps in SPS_LIBS:
                per_sps[sps] = run_aperture_sps(
                    state, label, r_max, sps,
                    n_bootstrap=n_bootstrap, force=force,
                    mask_weight=mask_weight,
                    fallback_n=50,
                )
            apertures[label] = dict(
                r_max=float(r_max), frac=float(frac),
                per_sps=per_sps,
                pool=pool_sps(per_sps),
            )
        tracks[track_label] = apertures

    pool_Re_masked = tracks["masked"]["Re"]["pool"]
    pool_Re_nomask = tracks["nomask"]["Re"]["pool"]
    pool_Re_soft   = tracks["softmask_w0p5"]["Re"]["pool"]
    budget = error_budget(pool_Re_masked)
    delta_mask = pool_Re_nomask["sigma_p50"] - pool_Re_masked["sigma_p50"]
    delta_soft = pool_Re_soft["sigma_p50"]   - pool_Re_masked["sigma_p50"]

    print(f"\n§6  Headline: σ_e(<R_e) — three tracks compared")
    print("─" * 70)
    def _show(name, p):
        print(f"  Track {name}:  σ_e = {p['sigma_p50']:.2f} "
              f"(+{p['sigma_p84']-p['sigma_p50']:.2f} "
              f"−{p['sigma_p50']-p['sigma_p16']:.2f}) km/s")
    _show("A (masked,    w=0.0)", pool_Re_masked)
    _show("B (nomask,    w=1.0)", pool_Re_nomask)
    _show("C (softmask,  w=0.5)", pool_Re_soft)
    print(f"  Δ(nomask − masked)   = {delta_mask:+.2f} km/s")
    print(f"  Δ(softmask − masked) = {delta_soft:+.2f} km/s")
    if abs(delta_mask) > 1e-6:
        print(f"  Linearity check: Δ_soft / Δ_nomask = "
              f"{delta_soft/delta_mask:.3f}  (expect 0.5 if linear in arc weight)")
    print(f"\n  Error budget (masked headline):")
    for k in ("stat", "I_shape", "mask", "frame", "centering"):
        print(f"    σ_{k:9s} = {budget[k]:.1f} km/s")
    print(f"  Quadrature total: ±{budget['total']:.1f} km/s")
    print(f"\n  HEADLINE: σ_e(<R_e) = {pool_Re_masked['sigma_p50']:.0f} ± {budget['total']:.0f} km/s")
    apertures = tracks["masked"]  # downstream save uses the headline track

    # ── §7  R_e source systematic (single track: masked, w=0.0) ─────────────
    # Re-run §6cum at three alternative R_e definitions so we can quote the
    # σ_e sensitivity to the choice of R_e itself. Always uses the headline
    # mask.
    print(f"\n§7  R_e systematic — 3 alternative R_e definitions (masked track)")
    print("─" * 70)
    Re_variants = (
        ("Re_F140W",  state["Re_140"]),
        ("Re_F200LP", state["Re_200"]),
        ("Re_CaHK",   state["Re_cahk"]),
    )
    re_sys = {}
    for var_label, r_max in Re_variants:
        print(f"  {var_label}:  R_max = {r_max:.3f}\"  "
              f"({(r_max - state['R_E'])*100/state['R_E']:+.1f}% vs headline mean)")
        per_sps = {}
        for sps in SPS_LIBS:
            per_sps[sps] = run_aperture_sps(
                state, var_label, float(r_max), sps,
                n_bootstrap=n_bootstrap, force=force,
                mask_weight=0.0,
                fallback_n=50,
            )
        re_sys[var_label] = dict(
            r_max=float(r_max), per_sps=per_sps, pool=pool_sps(per_sps),
        )

    print(f"\n  Comparison vs headline mean R_e (Track A, masked):")
    print(f"    {'source':<14} {'R_e [\"]':>9} {'σ_e [km/s]':>12} {'Δσ_e':>8}")
    hl_sigma = pool_Re_masked["sigma_p50"]
    print(f"    {'mean (paper)':<14} {state['R_E']:>9.3f} {hl_sigma:>12.2f} {0.0:>+8.2f}")
    for var_label, d in re_sys.items():
        s = d["pool"]["sigma_p50"]
        print(f"    {var_label:<14} {d['r_max']:>9.3f} {s:>12.2f} {s - hl_sigma:>+8.2f}")

    # Pack everything into a single npz for the notebook to load
    # Per-(label, sps, track) record of which N each fit was loaded at.
    # Important for the paper: the masked headline is N=500 throughout, but
    # the no-mask sensitivity may fall back to N=50 if a slow fit was killed.
    n_used = {}
    for tname, tdict in tracks.items():
        for lab in APERTURE_LABELS:
            for s in SPS_LIBS:
                n_used[f"{tname}/{lab}/{s}"] = int(tdict[lab]["per_sps"][s]["n_bootstrap"])
    print("\n  N_bootstrap actually used per (track / aperture / SPS):")
    for k, v in n_used.items():
        print(f"    {k:<26} N={v}")
    save = dict(
        z_systemic=Z_SYSTEMIC,
        ra_center=state["ra_center"], dec_center=state["dec_center"],
        d_centroid_arcsec=state["d_centroid"],
        cx_ifu=state["cx_ifu"], cy_ifu=state["cy_ifu"],
        Re_140=state["Re_140"], Re_200=state["Re_200"], R_E=state["R_E"],
        Re_cahk=state["Re_cahk"],
        kpc_per_arcsec=KPC_PER_ARCSEC,
        n_bootstrap=int(n_bootstrap),
        n_bootstrap_used=np.array(list(n_used.values())),
        n_bootstrap_used_keys=np.array(list(n_used.keys())),
        degrees=np.asarray(DEGREES),
        sps_libs=np.array(SPS_LIBS),
        sps_native_frame=np.array([SPS_NATIVE_FRAME[s] for s in SPS_LIBS]),
        aperture_labels=np.array(APERTURE_LABELS),
        aperture_fracs=np.array(APERTURE_FRACS),
        aperture_r_max=np.array([apertures[l]["r_max"] for l in APERTURE_LABELS]),
        # Per-aperture pooled posterior summary
        sigma_p16=np.array([apertures[l]["pool"]["sigma_p16"] for l in APERTURE_LABELS]),
        sigma_p50=np.array([apertures[l]["pool"]["sigma_p50"] for l in APERTURE_LABELS]),
        sigma_p84=np.array([apertures[l]["pool"]["sigma_p84"] for l in APERTURE_LABELS]),
        V_p16=np.array([apertures[l]["pool"]["V_p16"] for l in APERTURE_LABELS]),
        V_p50=np.array([apertures[l]["pool"]["V_p50"] for l in APERTURE_LABELS]),
        V_p84=np.array([apertures[l]["pool"]["V_p84"] for l in APERTURE_LABELS]),
        # Per-aperture per-SPS V_sys + sigma summaries (masked track)
        per_sps_summary=np.array(
            [[apertures[l]["pool"]["per_sps"][s]["sigma_p50"] for s in SPS_LIBS]
             for l in APERTURE_LABELS]
        ),
        per_sps_V_sys=np.array(
            [[apertures[l]["pool"]["per_sps"][s]["V_sys"] for s in SPS_LIBS]
             for l in APERTURE_LABELS]
        ),
        per_sps_frame=np.array(
            [[apertures[l]["pool"]["per_sps"][s]["frame"] for s in SPS_LIBS]
             for l in APERTURE_LABELS]
        ),
        # Pooled posterior samples for the headline aperture (masked track) — for figures
        sigma_samples_Re=apertures["Re"]["pool"]["sigma_samples"],
        sigma_samples_Re_2=apertures["Re_2"]["pool"]["sigma_samples"],
        # No-mask track summaries (sensitivity, mask_weight=1.0)
        nomask_sigma_p16=np.array([tracks["nomask"][l]["pool"]["sigma_p16"] for l in APERTURE_LABELS]),
        nomask_sigma_p50=np.array([tracks["nomask"][l]["pool"]["sigma_p50"] for l in APERTURE_LABELS]),
        nomask_sigma_p84=np.array([tracks["nomask"][l]["pool"]["sigma_p84"] for l in APERTURE_LABELS]),
        nomask_sigma_samples_Re=tracks["nomask"]["Re"]["pool"]["sigma_samples"],
        nomask_sigma_samples_Re_2=tracks["nomask"]["Re_2"]["pool"]["sigma_samples"],
        nomask_per_sps_summary=np.array(
            [[tracks["nomask"][l]["pool"]["per_sps"][s]["sigma_p50"] for s in SPS_LIBS]
             for l in APERTURE_LABELS]
        ),
        delta_mask_kms=float(delta_mask),
        # Soft-mask track summaries (interpolation, mask_weight=0.5)
        soft_sigma_p16=np.array([tracks["softmask_w0p5"][l]["pool"]["sigma_p16"] for l in APERTURE_LABELS]),
        soft_sigma_p50=np.array([tracks["softmask_w0p5"][l]["pool"]["sigma_p50"] for l in APERTURE_LABELS]),
        soft_sigma_p84=np.array([tracks["softmask_w0p5"][l]["pool"]["sigma_p84"] for l in APERTURE_LABELS]),
        soft_sigma_samples_Re=tracks["softmask_w0p5"]["Re"]["pool"]["sigma_samples"],
        soft_sigma_samples_Re_2=tracks["softmask_w0p5"]["Re_2"]["pool"]["sigma_samples"],
        soft_per_sps_summary=np.array(
            [[tracks["softmask_w0p5"][l]["pool"]["per_sps"][s]["sigma_p50"] for s in SPS_LIBS]
             for l in APERTURE_LABELS]
        ),
        delta_soft_kms=float(delta_soft),
        # Error budget (masked track is headline)
        budget_stat=budget["stat"], budget_I_shape=budget["I_shape"],
        budget_mask=budget["mask"], budget_frame=budget["frame"],
        budget_centering=budget["centering"], budget_total=budget["total"],
        # §7 R_e systematic: σ_e at each of three alt R_e definitions
        re_sys_labels=np.array([k for k, _ in Re_variants]),
        re_sys_r_max=np.array([re_sys[k]["r_max"] for k, _ in Re_variants]),
        re_sys_sigma_p16=np.array([re_sys[k]["pool"]["sigma_p16"] for k, _ in Re_variants]),
        re_sys_sigma_p50=np.array([re_sys[k]["pool"]["sigma_p50"] for k, _ in Re_variants]),
        re_sys_sigma_p84=np.array([re_sys[k]["pool"]["sigma_p84"] for k, _ in Re_variants]),
        re_sys_per_sps_sigma=np.array(
            [[re_sys[k]["pool"]["per_sps"][s]["sigma_p50"] for s in SPS_LIBS]
             for k, _ in Re_variants]
        ),
    )
    np.savez(RESULTS_NPZ, **save)
    print(f"\nSaved → {RESULTS_NPZ.relative_to(REPO)}")
    return save


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n_bootstrap", type=int, default=50,
                   help="bootstrap iterations per (aperture, SPS, degree); 50=smoke, 500=production")
    p.add_argument("--force", action="store_true",
                   help="ignore caches and rebuild from scratch")
    args = p.parse_args()
    main(args.n_bootstrap, args.force)
