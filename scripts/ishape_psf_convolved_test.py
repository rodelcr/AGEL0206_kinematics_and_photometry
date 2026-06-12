"""Does PSF-matching the HST I-weight maps to the KCWI seeing shift σ_e? (2026-06-12)

The I-shape systematic (run_isource_shape_sweep.py) uses HST F140W/F200LP as alternative
luminosity-weight maps, reprojected onto the IFU grid at NATIVE HST resolution — they are
NOT convolved to the KCWI ~1.27" seeing (nor are the Sérsic maps). This quantifies that gap:
re-measure σ_e with the HST I-weight (a) raw and (b) Gaussian-convolved to the 1.27" KCWI PSF
before reprojection, holding everything else at the headline (R_e=2.097", arc mask, FSPS).

The relevant number is the raw→PSF-convolved SHIFT (the resolution mismatch); the headline
IFU-white-light weight (already at KCWI resolution) is the resolution-matched anchor.

Output: results/ishape_psf_convolved_test.npz (+ printed table)
Usage: conda activate ISMgas; python scripts/ishape_psf_convolved_test.py
"""
import os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import map_coordinates, gaussian_filter
from ppxf.ppxf import ppxf

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, str(REPO)); os.chdir(REPO)
import scripts.final_sigma_e as fse
import scripts.run_wide_sigma_e as rws
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum

SEEING_FWHM = 1.27
R_E = 2.097
SPS = "fsps"


def reproject(img, wcs_hst, wcs_ifu, ny, nx):
    yy, xx = np.mgrid[0:ny, 0:nx]
    ra, dec = wcs_ifu.pixel_to_world_values(xx.ravel(), yy.ravel())
    xh, yh = wcs_hst.world_to_pixel_values(ra, dec)
    return np.clip(map_coordinates(img, [yh, xh], order=1, mode="constant",
                                   cval=0.0).reshape(ny, nx), 0, None)


def sigma_e_for_weight(state, Imap):
    """Median σ_e (FSPS, 15 deg) using Imap as the I-weight at R_e, arc-masked."""
    st = dict(state); st["ifu_band"] = np.clip(Imap, 0, None)
    flux, noise, n_kept, sn = fse.extract_aperture_spectrum(st, r_max=R_E, mask_weight=0.0)
    n_med = float(np.nanmedian(noise[noise > 0]))
    noise = np.where(noise > 0, noise, n_med * 0.1)
    lo, hi = rws.LAM_OBS_RANGE
    lr = (max(lo, rws.sps_safe_obs_min(SPS)), hi)
    inp = setup_ppxf_inputs_from_spectrum(flux, noise, state["hdr"], sps_name=SPS,
        z=fse.Z_SYSTEMIC, lam_obs_range=lr, lam_fit_range=lr,
        lam_range_temp=rws.SPS_LAM_RANGE_TEMP[SPS], verbose=False,
        frame_galaxy='auto', mask_balmer=False)
    gp = rws._apply_arc_mask(inp['goodpixels'], inp['lam_gal_rest']) if rws.ARC_MASK else inp['goodpixels']
    gp = rws._apply_bad_pixels_mask(gp, inp['lam_gal_rest'])
    sig = []
    for deg in rws.DEGREES:
        pp = ppxf(inp['sps'].templates, inp['galaxy'], inp['noise'], inp['velscale'],
                  inp['start'], goodpixels=gp, plot=False, moments=2, trig=False,
                  degree=int(deg), mdegree=0, lam=inp['lam_gal_rest'],
                  lam_temp=inp['lam_temp'], quiet=True)
        sig.append(pp.sol[1])
    return float(np.median(sig)), int(n_kept), float(sn)


def main():
    fse.IFU_FILE = "raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits"
    st = fse.load_setup()
    wcs_ifu = WCS(st["hdr"], naxis=2); ny, nx = st["ny"], st["nx"]

    maps = {"IFU white-light (headline, KCWI-res)": st["ifu_band"]}
    for band, path in (("F140W", fse.HST_F140W), ("F200LP", fse.HST_F200LP)):
        with fits.open(path) as h:
            img = h[0].data.astype(float); w = WCS(h[0].header)
        ps = float(np.abs(proj_plane_pixel_scales(w)[0]) * 3600)
        sig_pix = (SEEING_FWHM / 2.3548) / ps
        img_psf = gaussian_filter(np.nan_to_num(img), sig_pix, mode="nearest")
        maps[f"{band} raw (HST-res)"] = reproject(img, w, wcs_ifu, ny, nx)
        maps[f"{band} PSF-conv 1.27\""] = reproject(img_psf, w, wcs_ifu, ny, nx)
        print(f"  {band}: pixscale {ps:.3f}\"/pix → σ_psf {sig_pix:.1f} pix")

    print(f"\n  {'I-weight map':36s} {'σ_e (km/s)':>11s} {'Δ vs headline':>14s}")
    print("  " + "-" * 64)
    res = {}
    hl = None
    for name, Imap in maps.items():
        s, nk, sn = sigma_e_for_weight(st, Imap)
        if hl is None:
            hl = s
        res[name] = s
        print(f"  {name:36s} {s:>11.2f} {s-hl:>+13.2f}")

    print("\n  ── PSF-matching shift (raw → PSF-convolved) ──")
    for band in ("F140W", "F200LP"):
        raw = res[f"{band} raw (HST-res)"]; conv = res[f"{band} PSF-conv 1.27\""]
        print(f"    {band}: raw {raw:.2f} → PSF-conv {conv:.2f}  (Δ = {conv-raw:+.2f} km/s)")
    print(f"  (I-shape systematic = ±2.27 km/s; headline σ_e = 267.31)")

    np.savez("results/ishape_psf_convolved_test.npz",
             names=np.array(list(res.keys())), sigma_e=np.array(list(res.values())),
             seeing_fwhm=SEEING_FWHM, R_e=R_E)
    print("\n  → results/ishape_psf_convolved_test.npz")


if __name__ == "__main__":
    main()
