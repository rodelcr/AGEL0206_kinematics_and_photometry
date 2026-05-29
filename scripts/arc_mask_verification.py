"""Verify the expert-painted F200LP arc-light mask against two independent,
objective, physically-motivated selectors, and characterise the masking-vs-noise
tradeoff via an S/N-regime sweep.

Motivation
----------
The adopted arc masks (`*_mask.fits`) were hand-painted pixel-by-pixel. A referee
cannot reproduce a hand mask. This script shows the hand mask is recoverable from
two objective criteria, and quantifies *how much* arc masking is actually needed
relative to the image noise (the S/N sweep):

  (A) COLOR  (F200LP - F140W): the lensed source is a blue z=1.302 star-forming
      galaxy; the deflector is a red z=0.676 elliptical. Arc-dominated pixels are
      *bluer* than the deflector core and sit away from it.
  (B) SERSIC-RESIDUAL: subtract the fitted 2D Sersic deflector model (reused from
      scripts.sersic_total_photometry) from the image; positive excess flags the
      non-smooth (arc / companion) light.

Over-masking guard (cf. nb07d: dilating the F200 mask inflates sigma — noise-driven):
every candidate mask requires S/N above sky, only positive Sersic residual is
flagged, and we report the masked-flux fraction attributable to the smooth Sersic
deflector model (rising fraction = eating deflector light = over-masking onset).

Cross-grid alignment uses scipy.ndimage.map_coordinates + WCS (the repo's own
no-dependency reprojection, as in measure_Re.build_psf_contamination_map) because
the `reproject` package is not installed in the ISMgas env.

Usage
-----
    conda activate ISMgas
    python scripts/arc_mask_verification.py            # full run, save npz + figs
    python scripts/arc_mask_verification.py --quick    # skip Re (faster smoke)

Outputs
-------
    results/arc_mask_verification.npz
    results/figures/arcmask_{F200LP,F140W}_maps.png
    results/figures/arcmask_snr_sweep.png
"""

import os
import sys
import json
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.modeling.models import Sersic2D
from scipy.ndimage import map_coordinates

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

# Reuse tested code: image loaders, centroid, Sersic fit
from scripts.sersic_total_photometry import load_hst, find_center_2dg

VDI = "../velocity_dispersion_from_IFU"
RA_DEFL, DEC_DEFL = 31.55611, -1.23817
R_E_INIT = 2.305
FIG_DIR = "results/figures"
PARAM_DIR = "example_outputs"

# Canonical PHOTFLAM/PHOTPLAM (stripped from cutouts; from CLAUDE.md / sersic_total_photometry)
BANDS = {
    "F200LP": dict(
        image=f"{VDI}/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits",
        mask=f"{VDI}/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits",
        params=f"{PARAM_DIR}/AGEL020613-011417A_F200LP_WFC3_cutout_L3_params.json",
        photflam=5.1851e-20, photplam=4923.48, pivot_AA=4972.0, ap_mag=22.6126,
    ),
    "F140W": dict(
        image=f"{VDI}/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits",
        mask=f"{VDI}/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits",
        params=f"{PARAM_DIR}/AGEL020613-011417A_F140W_WFC3_cutout_L3_params.json",
        photflam=1.4737e-20, photplam=13922.9, pivot_AA=13923.0, ap_mag=19.1335,
    ),
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def load_band(name):
    """Image, WCS, pix scale ("/pix), AB ZP (per count-rate), expert mask."""
    cfg = BANDS[name]
    b = load_hst(cfg["image"], cfg["mask"], cfg["photflam"], cfg["photplam"])
    b["name"] = name
    b["cfg"] = cfg
    cx, cy, method = find_center_2dg(b["img"], b["mask"], b["wcs"],
                                     RA_DEFL, DEC_DEFL, box_arcsec=3.0,
                                     pix_scale=b["pix_scale"])
    b["cx"], b["cy"] = cx, cy
    return b


def sky_sigma(img, cx, cy, ps, r_in=6.0, r_out=9.0):
    """Robust sky sigma (MAD) from an arc-free outer ring."""
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - cx, yy - cy) * ps
    ring = (r > r_in) & (r < r_out) & np.isfinite(img)
    s = img[ring]
    med = np.median(s)
    mad = np.median(np.abs(s - med)) * 1.4826
    return float(med), float(mad)


def reproject_intensity(src_img, src_wcs, src_ps, dst_wcs, dst_shape):
    """Reproject a *surface-brightness* (per-arcsec^2) field from src grid onto
    dst grid via WCS + bilinear map_coordinates. Returns dst-grid SB array."""
    src_sb = np.nan_to_num(src_img) / (src_ps ** 2)  # cts / arcsec^2
    yy, xx = np.mgrid[:dst_shape[0], :dst_shape[1]]
    ra, dec = dst_wcs.pixel_to_world_values(xx.ravel(), yy.ravel())
    xs, ys = src_wcs.world_to_pixel_values(ra, dec)
    out = map_coordinates(src_sb, [ys, xs], order=1, mode="constant", cval=np.nan)
    return out.reshape(dst_shape)


def color_map(primary, other):
    """AB color in a FIXED blue-minus-red convention (m_F200LP - m_F140W), on the
    primary grid, regardless of which band is primary. Arc (blue z=1.302 source)
    is therefore always at *more negative* color than the red deflector core, on
    both grids. Finite only where both bands have positive surface brightness;
    constant ZP offsets are harmless (we threshold on the contrast vs the core)."""
    p_sb = np.nan_to_num(primary["img"]) / (primary["pix_scale"] ** 2)
    o_sb = reproject_intensity(other["img"], other["wcs"], other["pix_scale"],
                               primary["wcs"], primary["img"].shape)
    good = (p_sb > 0) & (o_sb > 0) & np.isfinite(o_sb)
    m_p = -2.5 * np.log10(np.where(good, p_sb, np.nan)) + primary["zp_ab"]
    m_o = -2.5 * np.log10(np.where(good, o_sb, np.nan)) + other["zp_ab"]
    # Order so the bluer (shorter pivot) band is the minuend → arc is more negative.
    if primary["cfg"]["pivot_AA"] <= other["cfg"]["pivot_AA"]:
        c = m_p - m_o   # primary is the blue band
    else:
        c = m_o - m_p   # other (reprojected) is the blue band
    color = np.full(primary["img"].shape, np.nan)
    color[good] = c[good]
    return color, p_sb, o_sb


def _fit_sersic_local(img, mask, center, r_eff_init_arcsec, box_arcsec, ps,
                      ellip_init, theta_init):
    """Multi-init Sersic2D fit with FREE position angle, seeded from the aperture
    geometry. Frees theta (unlike sersic_total_photometry / run_isource_shape_sweep,
    which freeze theta=0 and rail ellip→0 on the round core). Returns (fit,(x1,y1))."""
    from astropy.modeling.fitting import LevMarLSQFitter
    xc, yc = center
    half = int(np.ceil(box_arcsec / ps))
    y1, y2 = max(0, int(yc) - half), min(img.shape[0], int(yc) + half + 1)
    x1, x2 = max(0, int(xc) - half), min(img.shape[1], int(xc) + half + 1)
    sub = np.nan_to_num(img[y1:y2, x1:x2].astype(float))
    sub_mask = mask[y1:y2, x1:x2]
    # sky-subtract from outer ring
    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    r = np.hypot(xx - (xc - x1), yy - (yc - y1)) * ps
    ring = (r > 2 * r_eff_init_arcsec) & ~sub_mask
    sub = sub - (np.median(sub[ring]) if ring.any() else 0.0)
    r_eff_pix = r_eff_init_arcsec / ps
    amp_init = float(np.nanpercentile(sub[~sub_mask], 98)) if (~sub_mask).any() else 1.0
    weights = (~sub_mask).astype(float)
    best, best_chi2 = None, np.inf
    for n_init in (1.5, 2.5, 4.0):
        for e0 in (ellip_init, max(ellip_init - 0.1, 0.05)):
            init = Sersic2D(
                amplitude=amp_init, r_eff=r_eff_pix * 0.7, n=n_init,
                x_0=xc - x1, y_0=yc - y1, ellip=e0, theta=theta_init,
                bounds={"n": (1.0, 6.0), "r_eff": (r_eff_pix * 0.3, r_eff_pix * 2.0),
                        "ellip": (0.0, 0.7), "amplitude": (1e-4, 1e5),
                        "x_0": (xc - x1 - 3, xc - x1 + 3),
                        "y_0": (yc - y1 - 3, yc - y1 + 3)})
            try:
                cand = LevMarLSQFitter()(init, xx, yy, sub, weights=weights, maxiter=800)
                chi2 = float(np.sum(((sub - cand(xx, yy)) * weights) ** 2))
                if chi2 < best_chi2:
                    best_chi2, best = chi2, cand
            except Exception:
                pass
    return best, (x1, y1)


def render_sersic_full(b, ellip_init, theta_init):
    """Fit Sersic2D (box-limited, free PA seeded from aperture geometry) then
    render on the FULL band grid. Returns (model_full, fit)."""
    fit, (x1, y1) = _fit_sersic_local(
        b["img"], b["mask"], (b["cx"], b["cy"]), R_E_INIT, 4.0,
        b["pix_scale"], ellip_init, theta_init)
    yy, xx = np.mgrid[:b["img"].shape[0], :b["img"].shape[1]]
    model_full = np.clip(np.asarray(fit(xx - x1, yy - y1), dtype=float), 0, None)
    return model_full, fit


def core_exclude(shape, cx, cy, ps, r_arcsec):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    return (np.hypot(xx - cx, yy - cy) * ps) < r_arcsec


# ---------------------------------------------------------------------------
# candidate masks
# ---------------------------------------------------------------------------
def color_arc_mask(color, defl_color, dcolor_thresh, flux, sky_sig, snr_thresh,
                   core_excl):
    """Bluer than deflector core by dcolor_thresh AND S/N>snr_thresh AND outside core."""
    bluer = np.isfinite(color) & (color < (defl_color - dcolor_thresh))
    detected = flux > (snr_thresh * sky_sig)
    return bluer & detected & ~core_excl


def sersic_residual_mask(residual, sky_sig, k_sigma, flux, snr_thresh, core_excl):
    """Positive Sersic-residual excess > k_sigma*sky AND S/N>snr_thresh AND outside core."""
    excess = residual > (k_sigma * sky_sig)
    detected = flux > (snr_thresh * sky_sig)
    return excess & detected & ~core_excl


# ---------------------------------------------------------------------------
# comparison + over-masking
# ---------------------------------------------------------------------------
def compare_masks(expert, auto, region=None):
    """IoU, expert-recovery, auto-excess. Optionally restricted to a region."""
    e = expert.astype(bool)
    a = auto.astype(bool)
    if region is not None:
        e = e & region
        a = a & region
    inter = (e & a).sum()
    union = (e | a).sum()
    iou = inter / union if union else 0.0
    recovery = inter / e.sum() if e.sum() else 0.0           # frac of expert recovered
    excess = (a & ~e).sum() / a.sum() if a.sum() else 0.0    # frac of auto not in expert
    return dict(iou=float(iou), recovery=float(recovery), excess=float(excess),
                n_expert=int(e.sum()), n_auto=int(a.sum()), n_inter=int(inter))


def overmasking_metrics(mask, model_full, img, sky_med):
    """Fraction of masked flux attributable to the smooth Sersic deflector model
    vs. genuine excess. Rising model-fraction = eating deflector light."""
    m = mask.astype(bool)
    if m.sum() == 0:
        return dict(n=0, model_flux_frac=0.0, masked_flux=0.0)
    data_sub = np.nan_to_num(img) - sky_med
    masked_data = np.clip(data_sub[m], 0, None).sum()
    masked_model = np.clip(model_full[m], 0, None).sum()
    frac = float(masked_model / masked_data) if masked_data > 0 else float("nan")
    return dict(n=int(m.sum()), model_flux_frac=frac, masked_flux=float(masked_data))


# ---------------------------------------------------------------------------
# photometry: aperture AB mag under a given mask
# ---------------------------------------------------------------------------
def aperture_ab_mag(b, mask):
    """Elliptical-aperture AB mag using the params.json geometry + canonical ZP.
    Reproduces the formula in photometry_masking_HST.run_photometry_math but
    returns the value (and is headless-safe). mask: bool array (True=excluded)."""
    from photutils.aperture import EllipticalAperture, EllipticalAnnulus, ApertureStats
    with open(b["cfg"]["params"]) as f:
        p = json.load(f)
    pg = p["pixel_geometry"]
    xc, yc = b["wcs"].world_to_pixel_values(p["target_ra"], p["target_dec"])
    aper = EllipticalAperture((float(xc), float(yc)), a=pg["a"], b=pg["b"],
                              theta=pg["theta_rad"])
    ann = EllipticalAnnulus((float(xc), float(yc)),
                            a_in=pg["annulus_a_in"], a_out=pg["annulus_a_out"],
                            b_out=pg["annulus_a_out"] * (pg["b"] / pg["a"]),
                            theta=pg["theta_rad"])
    data = np.nan_to_num(b["img"])
    tot = None if mask is None else np.asarray(mask).astype(bool)
    ap_stats = ApertureStats(data, aper, mask=tot)
    ann_stats = ApertureStats(data, ann, mask=tot)
    net = ap_stats.sum - ann_stats.median * ap_stats.sum_aper_area.value
    photflam, photplam = b["cfg"]["photflam"], b["cfg"]["photplam"]
    ab_zp = -2.5 * np.log10(photflam) - 21.10 - 5 * np.log10(photplam) + 18.69
    return float(-2.5 * np.log10(net) + ab_zp) if net > 0 else float("nan")


# ---------------------------------------------------------------------------
# per-band analysis
# ---------------------------------------------------------------------------
SNR_REGIMES = [2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0]
KSIGMA_REGIMES = [2.0, 3.0, 5.0, 8.0]
DCOLOR_THRESH = 0.5     # mag bluer than deflector core to call "arc"
CORE_R_ARCSEC = 0.4     # never mask inside this radius (deflector peak)


def analyze_band(primary_name, other_name, do_re=True):
    print(f"\n{'='*70}\n  {primary_name}  (color vs {other_name})\n{'='*70}")
    primary = load_band(primary_name)
    other = load_band(other_name)
    img, ps, cx, cy = primary["img"], primary["pix_scale"], primary["cx"], primary["cy"]
    sky_med, sky_sig = sky_sigma(img, cx, cy, ps)
    print(f"  pix={ps:.3f}\"  center=({cx:.1f},{cy:.1f})  sky_med={sky_med:.3e}  sky_sig={sky_sig:.3e}")

    color, p_sb, o_sb = color_map(primary, other)
    excl = core_exclude(img.shape, cx, cy, ps, CORE_R_ARCSEC)
    # deflector core color = median within 0.4-0.8" annulus (avoid saturated peak / arc)
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    rr = np.hypot(xx - cx, yy - cy) * ps
    core_ring = (rr > 0.4) & (rr < 0.8) & np.isfinite(color)
    defl_color = float(np.median(color[core_ring]))
    print(f"  deflector core color (m_F200LP-m_F140W, fixed) = {defl_color:.3f}")

    with open(primary["cfg"]["params"]) as f:
        _pg = json.load(f)
    ellip_init = float(_pg["geometry"]["ellipticity"])
    theta_init = float(_pg["pixel_geometry"]["theta_rad"])
    model_full, fit = render_sersic_full(primary, ellip_init, theta_init)
    residual = np.nan_to_num(img) - sky_med - model_full
    print(f"  Sersic2D: n={fit.n.value:.2f} r_eff={fit.r_eff.value*ps:.2f}\" "
          f"ellip={fit.ellip.value:.2f} theta={np.degrees(fit.theta.value)%180:.0f}deg")

    expert = primary["mask"].astype(bool)
    flux = np.nan_to_num(img) - sky_med
    region_ap = rr < (primary_name == "F200LP" and 1.5 or 3.0)  # aperture-ish footprint

    # --- S/N regime sweep: color mask ---
    sweep_color = []
    for snr in SNR_REGIMES:
        m = color_arc_mask(color, defl_color, DCOLOR_THRESH, flux, sky_sig, snr, excl)
        cmp_full = compare_masks(expert, m)
        cmp_ap = compare_masks(expert, m, region_ap)
        om = overmasking_metrics(m, model_full, img, sky_med)
        dmag = aperture_ab_mag(primary, m) - primary["cfg"]["ap_mag"]
        sweep_color.append(dict(snr=snr, **cmp_full,
                                iou_ap=cmp_ap["iou"], model_flux_frac=om["model_flux_frac"],
                                dmag=dmag))
        print(f"   color snr>{snr:>4}: n={cmp_full['n_auto']:>5} IoU={cmp_full['iou']:.2f} "
              f"recov={cmp_full['recovery']:.2f} mdlfrac={om['model_flux_frac']:.2f} dmag={dmag:+.4f}")

    # --- S/N regime sweep: sersic residual mask ---
    sweep_sersic = []
    for k in KSIGMA_REGIMES:
        m = sersic_residual_mask(residual, sky_sig, k, flux, 3.0, excl)
        cmp_full = compare_masks(expert, m)
        om = overmasking_metrics(m, model_full, img, sky_med)
        dmag = aperture_ab_mag(primary, m) - primary["cfg"]["ap_mag"]
        sweep_sersic.append(dict(k_sigma=k, **cmp_full,
                                 model_flux_frac=om["model_flux_frac"], dmag=dmag))
        print(f"   sersic k>{k:>4}: n={cmp_full['n_auto']:>5} IoU={cmp_full['iou']:.2f} "
              f"recov={cmp_full['recovery']:.2f} mdlfrac={om['model_flux_frac']:.2f} dmag={dmag:+.4f}")

    # --- adopted automated masks at a fiducial regime (color snr>5, sersic k>3) ---
    m_color = color_arc_mask(color, defl_color, DCOLOR_THRESH, flux, sky_sig, 5.0, excl)
    m_sersic = sersic_residual_mask(residual, sky_sig, 3.0, flux, 3.0, excl)
    m_union = m_color | m_sersic

    # --- photometry under each mask ---
    masks = {"none": None, "expert": expert, "color": m_color,
             "sersic": m_sersic, "union": m_union}
    phot = {}
    for label, m in masks.items():
        mag = aperture_ab_mag(primary, m)
        rec = dict(ab_mag=mag, dmag=mag - primary["cfg"]["ap_mag"])
        if do_re:
            from scripts.measure_Re import hst_Re
            ov = None if (m is None) else np.asarray(m).astype(bool)
            try:
                re = hst_Re(primary_name, masking="proper",
                            mask_override=ov)["Re"] if ov is not None else \
                     hst_Re(primary_name, masking="none")["Re"]
                rec["Re"] = float(re)
            except Exception as e:
                rec["Re"] = float("nan")
                print(f"   [Re fail {label}: {e}]")
        phot[label] = rec
        print(f"   PHOT {label:>7}: AB={rec['ab_mag']:.4f}  dmag={rec['dmag']:+.4f}"
              + (f"  Re={rec.get('Re', float('nan')):.3f}\"" if do_re else ""))

    return dict(
        name=primary_name, other=other_name, ps=ps, cx=cx, cy=cy,
        sky_med=sky_med, sky_sig=sky_sig, defl_color=defl_color,
        color=color, model_full=model_full, residual=residual,
        expert=expert, m_color=m_color, m_sersic=m_sersic, m_union=m_union,
        sweep_color=sweep_color, sweep_sersic=sweep_sersic, phot=phot,
        sersic_n=float(fit.n.value), sersic_reff_arcsec=float(fit.r_eff.value * ps),
        sersic_ellip=float(fit.ellip.value),
    )


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def fig_band_maps(res):
    name = res["name"]
    img_cx, img_cy = res["cx"], res["cy"]
    fig, ax = plt.subplots(2, 4, figsize=(18, 9))
    color = res["color"]
    vmin, vmax = np.nanpercentile(color, [5, 95])
    ax[0, 0].imshow(color, origin="lower", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    # color is ALWAYS m_F200LP - m_F140W (fixed blue-red convention), on this grid
    ax[0, 0].set_title("color  $m_{F200LP}-m_{F140W}$ (AB SB)\nblue=arc, red=deflector")
    ax[0, 1].imshow(res["model_full"], origin="lower", cmap="magma")
    ax[0, 1].set_title("Sersic model")
    rv = np.nanpercentile(np.abs(res["residual"]), 99)
    ax[0, 2].imshow(res["residual"], origin="lower", cmap="RdBu_r", vmin=-rv, vmax=rv)
    ax[0, 2].set_title("residual (data-model)")
    ax[0, 3].imshow(res["expert"], origin="lower", cmap="Greys")
    ax[0, 3].set_title(f"expert mask (n={res['expert'].sum()})")
    for a, (m, t) in zip(ax[1], [(res["m_color"], "color mask"),
                                  (res["m_sersic"], "sersic mask"),
                                  (res["m_union"], "union"),
                                  (res["expert"], "expert")]):
        a.imshow(m, origin="lower", cmap="Reds")
        a.set_title(f"{t} (n={int(np.asarray(m).sum())})")
    # overlay expert outline on the union panel
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
        a.set_xlim(img_cx - 90, img_cx + 90); a.set_ylim(img_cy - 90, img_cy + 90)
    fig.suptitle(f"Arc-mask verification — {name}", fontsize=14)
    plt.tight_layout()
    out = f"{FIG_DIR}/arcmask_{name}_maps.png"
    plt.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print(f"  saved {out}")


def fig_subtraction_diagnostic(res):
    """Verify the Sersic subtraction is working: data / model / residual, residual
    in S/N units with the candidate-mask contour, an azimuthal radial profile
    (data vs model vs residual), and a residual-pixel histogram (arc-free region
    should be ~N(0,1); a coherent positive tail = the arc the subtraction reveals)."""
    name = res["name"]
    cx, cy, ps, sky_sig = res["cx"], res["cy"], res["ps"], res["sky_sig"]
    model = res["model_full"]
    # data (sky-subtracted) reconstructed from stored arrays: resid = data - sky - model
    data = res["residual"] + model
    resid = res["residual"]
    yy, xx = np.mgrid[:data.shape[0], :data.shape[1]]
    r = np.hypot(xx - cx, yy - cy) * ps

    fig = plt.figure(figsize=(20, 9))
    gs = fig.add_gridspec(2, 4, hspace=0.28, wspace=0.22)
    half = 90 if name == "F200LP" else 60
    sl = (slice(int(cy - half), int(cy + half)), slice(int(cx - half), int(cx + half)))
    vmax = float(np.nanpercentile(np.abs(data[sl]), 99))
    rv = 6 * sky_sig

    def panel(ax, arr, title, cmap, vmin, vmaxx):
        im = ax.imshow(arr[sl], origin="lower", cmap=cmap, vmin=vmin, vmax=vmaxx)
        ax.set_title(title, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046)

    panel(fig.add_subplot(gs[0, 0]), data, f"{name} data (sky-sub)", "magma", 0, vmax)
    panel(fig.add_subplot(gs[0, 1]), model, "Sersic model", "magma", 0, vmax)
    panel(fig.add_subplot(gs[0, 2]), resid, "residual = data − model", "RdBu_r", -rv, rv)
    # residual in S/N units, with the sersic mask contour
    axsn = fig.add_subplot(gs[0, 3])
    im = axsn.imshow((resid / sky_sig)[sl], origin="lower", cmap="RdBu_r", vmin=-6, vmax=6)
    axsn.contour(res["m_sersic"][sl].astype(float), levels=[0.5], colors="lime", linewidths=0.8)
    axsn.set_title("residual / σ_sky  (+sersic mask)", fontsize=11)
    axsn.set_xticks([]); axsn.set_yticks([]); plt.colorbar(im, ax=axsn, fraction=0.046)

    # azimuthal radial profile: data, model, |residual|
    redges = np.arange(0, (half * ps), max(ps, 0.05))
    rc, d_prof, m_prof, res_prof = [], [], [], []
    for j in range(len(redges) - 1):
        ring = (r >= redges[j]) & (r < redges[j + 1])
        if ring.sum() == 0:
            continue
        rc.append(0.5 * (redges[j] + redges[j + 1]))
        d_prof.append(np.median(data[ring])); m_prof.append(np.median(model[ring]))
        res_prof.append(np.median(resid[ring]))
    axp = fig.add_subplot(gs[1, :2])
    axp.plot(rc, d_prof, "k.-", label="data", ms=4)
    axp.plot(rc, m_prof, "r-", label="Sersic model", lw=1.5)
    axp.plot(rc, res_prof, "b-", label="residual (median)", lw=1)
    axp.axhline(0, color="grey", lw=0.5)
    axp.axhline(res["sky_sig"], color="b", ls=":", lw=0.8, label="±σ_sky")
    axp.axhline(-res["sky_sig"], color="b", ls=":", lw=0.8)
    axp.set_yscale("symlog", linthresh=sky_sig)
    axp.set(xlabel="R (arcsec)", ylabel="median I (e/s/pix)",
            title=f"{name} azimuthal profile — model tracks data; residual ~0 except at arc")
    axp.legend(fontsize=9); axp.grid(alpha=0.3)

    # residual histogram: arc-free annulus vs whole inner frame
    axh = fig.add_subplot(gs[1, 2:])
    arc_free = (r > 1.0) & (r < 4.0) & ~res["m_sersic"] & ~res["expert"]
    inner = (r < 4.0)
    bins = np.linspace(-6, 8, 80)
    axh.hist((resid[arc_free] / sky_sig), bins=bins, density=True, alpha=0.6,
             color="grey", label="arc-free region")
    axh.hist((resid[inner] / sky_sig), bins=bins, density=True, histtype="step",
             color="red", lw=1.5, label="all inner pixels")
    g = np.exp(-0.5 * bins ** 2) / np.sqrt(2 * np.pi)
    axh.plot(bins, g, "k--", lw=1, label="N(0,1)")
    axh.set(xlabel="residual / σ_sky", ylabel="density",
            title="arc-free ≈ N(0,1); positive tail = arc revealed")
    axh.legend(fontsize=9)

    out = f"{FIG_DIR}/arcmask_{name}_subtraction.png"
    plt.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print(f"  saved {out}")


def fig_snr_sweep(results):
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    for res in results:
        snr = [s["snr"] for s in res["sweep_color"]]
        ax[0].plot(snr, [s["n_auto"] for s in res["sweep_color"]], "o-", label=res["name"])
        ax[1].plot(snr, [abs(s["dmag"]) for s in res["sweep_color"]], "o-", label=res["name"])
        ax[2].plot(snr, [s["model_flux_frac"] for s in res["sweep_color"]], "o-", label=res["name"])
    ax[0].set(xlabel="S/N threshold", ylabel="masked px (color)", title="mask size vs S/N")
    ax[1].set(xlabel="S/N threshold", ylabel="|Δmag|", title="photometry shift vs S/N")
    ax[1].axhline(0.1, ls="--", c="r", label="10% flux err")
    ax[2].set(xlabel="S/N threshold", ylabel="masked-flux frac from Sersic model",
              title="over-masking onset")
    for a in ax:
        a.legend(); a.grid(alpha=0.3)
    plt.tight_layout()
    out = f"{FIG_DIR}/arcmask_snr_sweep.png"
    plt.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print(f"  saved {out}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip Re (faster smoke)")
    args = ap.parse_args()
    os.makedirs(FIG_DIR, exist_ok=True)

    res_f200 = analyze_band("F200LP", "F140W", do_re=not args.quick)
    res_f140 = analyze_band("F140W", "F200LP", do_re=not args.quick)
    results = [res_f200, res_f140]

    for r in results:
        fig_band_maps(r)
        fig_subtraction_diagnostic(r)
    fig_snr_sweep(results)

    # Save npz (object dtype for the nested dicts/arrays)
    save = {}
    for r in results:
        nm = r["name"]
        save[f"{nm}_color"] = r["color"]
        save[f"{nm}_model_full"] = r["model_full"]
        save[f"{nm}_residual"] = r["residual"]
        save[f"{nm}_expert"] = r["expert"]
        save[f"{nm}_m_color"] = r["m_color"]
        save[f"{nm}_m_sersic"] = r["m_sersic"]
        save[f"{nm}_m_union"] = r["m_union"]
        save[f"{nm}_sweep_color"] = np.array(r["sweep_color"], dtype=object)
        save[f"{nm}_sweep_sersic"] = np.array(r["sweep_sersic"], dtype=object)
        save[f"{nm}_phot"] = np.array([r["phot"]], dtype=object)
        save[f"{nm}_meta"] = np.array([dict(
            ps=r["ps"], sky_med=r["sky_med"], sky_sig=r["sky_sig"],
            defl_color=r["defl_color"], sersic_n=r["sersic_n"],
            sersic_reff_arcsec=r["sersic_reff_arcsec"], sersic_ellip=r["sersic_ellip"],
        )], dtype=object)
    save["config"] = np.array([dict(
        SNR_REGIMES=SNR_REGIMES, KSIGMA_REGIMES=KSIGMA_REGIMES,
        DCOLOR_THRESH=DCOLOR_THRESH, CORE_R_ARCSEC=CORE_R_ARCSEC,
        fiducial="color snr>5, sersic k>3",
    )], dtype=object)
    out = "results/arc_mask_verification.npz"
    np.savez(out, **save)
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    main()
