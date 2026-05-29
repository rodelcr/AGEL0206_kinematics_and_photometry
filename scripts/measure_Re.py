#!/usr/bin/env python
"""
Measure the effective radius (R_e) of the AGEL0206 deflector using multiple
data sources and masking strategies.

Sources:
  - IFU white-light (KCWI, 0.30"/spaxel)
  - HST F140W (WFC3/IR, 0.08"/pix)
  - HST F200LP (WFC3/UVIS, 0.05"/pix)   # cutout scale read from WCS (not hard-coded)

Masking strategies:
  - unmasked: no arc removal
  - zeroed: masked pixels set to 0 (biased — deflates annular means)
  - proper: masked pixels excluded via ApertureStats mask parameter
  - PSF-convolved: HST mask convolved with KCWI seeing PSF before applying

Usage:
  python scripts/measure_Re.py                  # run all, save results
  python scripts/measure_Re.py --plot           # also save comparison figure
"""
import numpy as np
import os
import argparse
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter

# ── Constants ──
RA_CENTER, DEC_CENTER = 31.55611, -1.23817
Z_LENS = 0.675
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
KPC_PER_ARCSEC = COSMO.kpc_proper_per_arcmin(Z_LENS).value / 60

# Paths (relative to repo root)
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
IFU_FILE = os.path.join(REPO, 'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits')
VDI = os.path.join(os.path.dirname(REPO), 'velocity_dispersion_from_IFU')
HST_FILES = {
    'F140W': {
        'image': os.path.join(VDI, 'AGEL020613-011417A_F140W_WFC3_cutout_L3.fits'),
        'mask': os.path.join(VDI, 'AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits'),
        'pixscale': 0.08,
    },
    'F200LP': {
        'image': os.path.join(VDI, 'AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits'),
        'mask': os.path.join(VDI, 'AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits'),
        'pixscale': 0.05,  # F200LP cutout is 0.05"/pix (was wrongly 0.08); hst_Re now reads WCS
    },
}


def measure_Re_from_profile(r_arcsec, intensity):
    """
    Half-light radius from a 1D radial intensity profile.
    Integrates I(r) * 2*pi*r dr and finds the radius enclosing 50% of flux.

    Returns: (Re, cumulative_flux_array, total_flux)
    """
    order = np.argsort(r_arcsec)
    r = r_arcsec[order]
    I = intensity[order]

    if r[0] > 0:
        r = np.concatenate([[0], r])
        I = np.concatenate([[I[0]], I])

    integrand = I * 2 * np.pi * r
    cumulative = np.concatenate([[0], cumulative_trapezoid(integrand, r)])
    total = cumulative[-1]

    if total <= 0:
        return np.nan, cumulative, total

    Re = np.interp(total / 2, cumulative, r)
    return Re, cumulative, total


# ─────────────────────────────────────────────
# IFU white-light R_e
# ─────────────────────────────────────────────

def load_ifu_cube():
    """Load IFU cube and return (cube, header, wavelength array)."""
    with fits.open(IFU_FILE) as hdul:
        hdr = hdul[0].header
        cube = np.asarray(hdul[0].data, dtype=float)
    crval = hdr['CRVAL3']
    cdelt = hdr['CD3_3']
    crpix = hdr.get('CRPIX3', 1.0)
    lam = crval + cdelt * (np.arange(hdr['NAXIS3']) + 1 - crpix)
    return cube, hdr, lam


def ifu_white_light_Re(cube, hdr, mask_contamination=None,
                       contamination_threshold=None):
    """
    R_e from the IFU white-light image.

    Parameters
    ----------
    mask_contamination : 2D array or None
        Per-spaxel contamination fraction (0-1). If None, no masking.
    contamination_threshold : float or None
        Spaxels with contamination > threshold are excluded.
        If None but mask_contamination is provided, uses 0.15.

    Returns dict with keys: Re, r_profile, I_profile, method, label
    """
    # Wide sub-image
    y_wide, x_wide = slice(35, 75), slice(30, 70)
    cube_wide = cube[:, y_wide, x_wide]
    wl = np.sum(cube_wide, axis=0)
    ny, nx = wl.shape

    # Center from smoothed peak
    wl_s = gaussian_filter(np.sum(cube, axis=0), sigma=2)
    cy, cx = np.unravel_index(np.argmax(wl_s), wl_s.shape)
    dx = cx - x_wide.start
    dy = cy - y_wide.start

    # Radial distances in arcsec
    wcs_2d = WCS(hdr, naxis=2)
    yy, xx = np.mgrid[y_wide.start:y_wide.stop, x_wide.start:x_wide.stop]
    ra, dec = wcs_2d.pixel_to_world_values(xx.ravel(), yy.ravel())
    ra_c, dec_c = wcs_2d.pixel_to_world_values(cx, cy)
    dra = (ra.reshape(ny, nx) - ra_c) * np.cos(np.radians(dec_c)) * 3600
    ddec = (dec.reshape(ny, nx) - dec_c) * 3600
    r = np.sqrt(dra**2 + ddec**2)

    # Masking
    if mask_contamination is not None:
        contam = mask_contamination[y_wide, x_wide]
        if contamination_threshold is None:
            contamination_threshold = 0.15
        bad = contam > contamination_threshold
        label = f'IFU white-light (masked >{contamination_threshold:.0%})'
    else:
        bad = np.zeros((ny, nx), dtype=bool)
        label = 'IFU white-light (unmasked)'

    # Radial profile in annuli, excluding masked spaxels
    r_edges = np.arange(0, 8, 0.3)
    r_prof, I_prof = [], []
    for j in range(len(r_edges) - 1):
        ann = (r >= r_edges[j]) & (r < r_edges[j + 1]) & ~bad
        if np.sum(ann) > 0:
            r_prof.append((r_edges[j] + r_edges[j + 1]) / 2)
            I_prof.append(np.mean(wl[ann]))

    r_prof = np.array(r_prof)
    I_prof = np.array(I_prof)
    Re, cog, total = measure_Re_from_profile(r_prof, I_prof)

    return {
        'Re': Re, 'r_profile': r_prof, 'I_profile': I_prof,
        'cog': cog, 'total': total,
        'method': 'IFU', 'label': label,
    }


# ─────────────────────────────────────────────
# HST R_e
# ─────────────────────────────────────────────

def hst_Re(band, masking='proper', mask_override=None):
    """
    R_e from an HST image.

    Parameters
    ----------
    band : str
        'F140W' or 'F200LP'
    masking : str
        'none'   — no masking
        'zeroed' — set masked pixels to 0 (old approach)
        'proper' — exclude masked pixels via ApertureStats mask parameter
    mask_override : ndarray or None
        If given, use this boolean mask instead of loading the band's default
        `_mask.fits`. Used by scripts/arc_mask_verification.py to test the
        color- and Sersic-residual-derived masks. Must match the image shape.

    Returns dict with Re, r_profile, I_profile, method, label
    """
    from photutils.aperture import CircularAnnulus, ApertureStats

    info = HST_FILES[band]
    with fits.open(info['image']) as hdul:
        img = hdul[0].data
        wcs = WCS(hdul[0].header)

    x_c, y_c = wcs.world_to_pixel_values(RA_CENTER, DEC_CENTER)
    x_c, y_c = float(x_c), float(y_c)
    # Read the pixel scale from the WCS, not the HST_FILES literal: the F200LP
    # cutout is 0.05"/pix (the dict previously hard-coded 0.08 for both bands,
    # which biased the *absolute* F200LP R_e high by 0.08/0.05 = 1.6x). The
    # headline R_e (scripts/final_sigma_e.py) already reads pix scale from the
    # WCS and is unaffected; this only corrects the measure_Re diagnostic.
    pscale = float(np.abs(proj_plane_pixel_scales(wcs)[0])) * 3600.0

    # Load mask (override takes precedence over the band's default _mask.fits)
    if mask_override is not None:
        mask = np.asarray(mask_override).astype(bool)
        if mask.shape != img.shape:
            raise ValueError(f"mask_override shape {mask.shape} != image {img.shape}")
    else:
        try:
            mask = fits.getdata(info['mask']).astype(bool)
        except FileNotFoundError:
            mask = None

    img_clean = np.nan_to_num(img, nan=0.0)

    # Prepare data and mask for ApertureStats
    if masking == 'none' or mask is None:
        data = img_clean
        ap_mask = None
        label = f'HST {band} (unmasked)'
    elif masking == 'zeroed':
        data = np.where(~mask, img_clean, 0)
        ap_mask = None
        label = f'HST {band} (zeroed)'
    elif masking == 'proper':
        data = img_clean
        ap_mask = mask
        label = f'HST {band} (proper mask)'
    else:
        raise ValueError(f"Unknown masking: {masking}")

    # Radial profile
    r_edges = np.arange(1, 80, 1)  # pixels
    r_prof, I_prof = [], []
    for j in range(len(r_edges) - 1):
        aper = CircularAnnulus((x_c, y_c), r_in=r_edges[j], r_out=r_edges[j + 1])
        if ap_mask is not None:
            stats = ApertureStats(data, aper, mask=ap_mask)
        else:
            stats = ApertureStats(data, aper)
        if np.isfinite(stats.mean) and stats.mean > 0:
            r_prof.append((r_edges[j] + r_edges[j + 1]) / 2 * pscale)
            I_prof.append(stats.mean)

    r_prof = np.array(r_prof)
    I_prof = np.array(I_prof)
    Re, cog, total = measure_Re_from_profile(r_prof, I_prof)

    return {
        'Re': Re, 'r_profile': r_prof, 'I_profile': I_prof,
        'cog': cog, 'total': total,
        'method': f'HST {band}', 'label': label,
    }


# ─────────────────────────────────────────────
# PSF-convolved contamination map for IFU masking
# ─────────────────────────────────────────────

def build_psf_contamination_map(hdr_ifu, seeing_fwhm=1.27):
    """
    Build a per-spaxel contamination fraction by reprojecting the HST
    arc mask onto the IFU grid after convolving with the KCWI seeing PSF.

    Returns (contamination_map, mask_convolved_hst) or (None, None).
    """
    from scipy.ndimage import map_coordinates

    # Use F140W mask (best arc coverage)
    info = HST_FILES['F140W']
    try:
        mask_hst = fits.getdata(info['mask']).astype(float)
        with fits.open(info['image']) as h:
            wcs_hst = WCS(h[0].header)
    except FileNotFoundError:
        return None, None

    # Convolve with KCWI PSF
    psf_sigma_hst = (seeing_fwhm / info['pixscale']) / 2.355
    mask_conv = gaussian_filter(mask_hst, sigma=psf_sigma_hst)

    # Reproject to IFU grid
    wcs_ifu = WCS(hdr_ifu, naxis=2)
    yy, xx = np.mgrid[0:100, 0:100]
    ra_ifu, dec_ifu = wcs_ifu.pixel_to_world_values(xx.ravel(), yy.ravel())
    x_hst, y_hst = wcs_hst.world_to_pixel_values(ra_ifu, dec_ifu)

    contam = map_coordinates(mask_conv,
                             [y_hst.reshape(100, 100).ravel(),
                              x_hst.reshape(100, 100).ravel()],
                             order=1, mode='constant', cval=0.0)
    contam = contam.reshape(100, 100)

    return contam, mask_conv


# ─────────────────────────────────────────────
# Run all measurements
# ─────────────────────────────────────────────

def run_all(verbose=True):
    """Run all R_e measurements. Returns list of result dicts."""
    results = []

    # Load IFU
    if verbose:
        print("Loading IFU cube...")
    cube, hdr, lam = load_ifu_cube()

    # IFU unmasked
    r = ifu_white_light_Re(cube, hdr)
    results.append(r)
    if verbose:
        print(f"  {r['label']}: R_e = {r['Re']:.2f}\" = {r['Re']*KPC_PER_ARCSEC:.1f} kpc")

    # IFU with PSF-convolved mask at several thresholds
    contam, _ = build_psf_contamination_map(hdr, seeing_fwhm=1.27)
    if contam is not None:
        for thresh in [0.05, 0.15, 0.30]:
            r = ifu_white_light_Re(cube, hdr, mask_contamination=contam,
                                   contamination_threshold=thresh)
            results.append(r)
            if verbose:
                print(f"  {r['label']}: R_e = {r['Re']:.2f}\" = {r['Re']*KPC_PER_ARCSEC:.1f} kpc")

    # HST bands
    for band in ['F140W', 'F200LP']:
        for masking in ['none', 'zeroed', 'proper']:
            try:
                r = hst_Re(band, masking=masking)
                results.append(r)
                if verbose:
                    print(f"  {r['label']}: R_e = {r['Re']:.2f}\" = {r['Re']*KPC_PER_ARCSEC:.1f} kpc")
            except FileNotFoundError:
                if verbose:
                    print(f"  HST {band} ({masking}): file not found, skipping")

    return results


def save_results(results, outdir=None):
    """Save R_e results to .npz file."""
    if outdir is None:
        outdir = os.path.join(REPO, 'results')
    os.makedirs(outdir, exist_ok=True)

    labels = [r['label'] for r in results]
    Re_values = [r['Re'] for r in results]
    methods = [r['method'] for r in results]

    np.savez(os.path.join(outdir, 'Re_measurements.npz'),
             labels=labels, Re_values=Re_values, methods=methods,
             kpc_per_arcsec=KPC_PER_ARCSEC, z_lens=Z_LENS,
             results=results)
    print(f"\nSaved results/Re_measurements.npz")


def plot_comparison(results, outdir=None):
    """Comparison figure: radial profiles + curve of growth + R_e bar chart."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if outdir is None:
        outdir = os.path.join(REPO, 'figures')
    os.makedirs(outdir, exist_ok=True)

    plt.rcParams['figure.facecolor'] = 'white'
    plt.rc('font', family='serif', size=12)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Color scheme by data source
    colors = {
        'IFU': 'C0',
        'HST F140W': 'C1',
        'HST F200LP': 'C2',
    }
    linestyles = {
        'unmasked': '-',
        'none': '-',
        'zeroed': '--',
        'proper': '-',
        'masked': ':',
    }

    # Panel 1: Radial profiles (normalized)
    ax = axes[0]
    for r in results:
        c = colors.get(r['method'], 'gray')
        # Determine linestyle from label
        if 'unmasked' in r['label'] or 'none' in r['label']:
            ls = '-'
            alpha = 0.4
        elif 'zeroed' in r['label']:
            ls = '--'
            alpha = 0.6
        elif 'proper' in r['label']:
            ls = '-'
            alpha = 1.0
        else:
            ls = ':'
            alpha = 0.7

        rp = r['r_profile']
        Ip = r['I_profile'] / r['I_profile'][0] if len(r['I_profile']) > 0 else []
        short_label = r['label'].replace('HST ', '').replace('IFU white-light ', 'IFU ')
        ax.plot(rp, Ip, color=c, ls=ls, alpha=alpha, lw=1.5, label=short_label)
        ax.axvline(r['Re'], color=c, ls=ls, alpha=alpha * 0.5, lw=0.8)

    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Normalized intensity')
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 1.5)
    ax.set_xlim(0, 6)
    ax.set_title('Radial light profiles')
    ax.legend(fontsize=7, loc='upper right', ncol=1)
    ax.grid(alpha=0.3)

    # Panel 2: Curves of growth
    ax = axes[1]
    for r in results:
        c = colors.get(r['method'], 'gray')
        if 'unmasked' in r['label'] or 'none' in r['label']:
            ls, alpha = '-', 0.4
        elif 'zeroed' in r['label']:
            ls, alpha = '--', 0.6
        elif 'proper' in r['label']:
            ls, alpha = '-', 1.0
        else:
            ls, alpha = ':', 0.7

        rp = r['r_profile']
        cog = r['cog']
        total = r['total']
        if total > 0:
            r_cog = np.concatenate([[0], rp])
            ax.plot(r_cog, cog / total, color=c, ls=ls, alpha=alpha, lw=1.5)
            ax.axvline(r['Re'], color=c, ls=ls, alpha=alpha * 0.5, lw=0.8)

    ax.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Cumulative flux fraction')
    ax.set_title('Curve of growth')
    ax.set_xlim(0, 6)
    ax.grid(alpha=0.3)

    # Panel 3: Bar chart of R_e values
    ax = axes[2]
    labels_short = []
    Re_vals = []
    bar_colors = []
    hatches = []
    for r in results:
        lbl = r['label'].replace('HST ', '').replace('IFU white-light ', 'IFU ')
        labels_short.append(lbl)
        Re_vals.append(r['Re'])
        bar_colors.append(colors.get(r['method'], 'gray'))
        if 'zeroed' in r['label']:
            hatches.append('//')
        elif 'unmasked' in r['label'] or 'none' in r['label']:
            hatches.append('..')
        else:
            hatches.append('')

    y_pos = np.arange(len(labels_short))
    bars = ax.barh(y_pos, Re_vals, color=bar_colors, alpha=0.7, edgecolor='black', lw=0.5)
    for bar, h in zip(bars, hatches):
        bar.set_hatch(h)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels_short, fontsize=8)
    ax.set_xlabel(r'$R_e$ (arcsec)')
    ax.set_title(r'$R_e$ comparison')
    ax.grid(alpha=0.3, axis='x')

    # Add kpc axis on top
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim()[0] * KPC_PER_ARCSEC, ax.get_xlim()[1] * KPC_PER_ARCSEC)
    ax2.set_xlabel('$R_e$ (kpc)')

    plt.tight_layout()
    outpath = os.path.join(outdir, 'Re_comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Saved {outpath}")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Measure R_e of AGEL0206 deflector')
    parser.add_argument('--plot', action='store_true', help='Save comparison figure')
    args = parser.parse_args()

    print("=" * 60)
    print("R_e measurements for AGEL J020613-011417 deflector")
    print("=" * 60)

    results = run_all(verbose=True)

    print(f"\n{'Label':<45} {'Re (arcsec)':>12} {'Re (kpc)':>10}")
    print("-" * 70)
    for r in results:
        print(f"{r['label']:<45} {r['Re']:>10.2f}\"  {r['Re']*KPC_PER_ARCSEC:>8.1f}")

    save_results(results)

    if args.plot:
        plot_comparison(results)
