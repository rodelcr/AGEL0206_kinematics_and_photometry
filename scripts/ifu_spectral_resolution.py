#!/usr/bin/env python
"""
Quick diagnostic: spectral resolution of the KCWI IFU data cube.

Reports:
  - DISPSCAL header value (instrumental sigma in pixels per resolution element)
  - Per-pixel wavelength step (A/pix)
  - Instrumental FWHM (A) as a function of wavelength
  - Resolving power R = lambda / FWHM
  - Instrumental velocity dispersion sigma_v_inst = c * (FWHM / lambda) / 2.355  (km/s)
  - Native cube velocity step c * d ln(lam)  (km/s/pix, used as ppxf velscale)
  - Same quantities reported at the ppxf fitting band centroid (~7000 A)

Usage:
  python scripts/ifu_spectral_resolution.py [cube.fits]
  # if no path given, uses the default AGEL0206 KCWI cube in repo root
"""
import os
import sys
import numpy as np
from astropy.io import fits

C_KMS = 299792.458


def main():
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    default_cube = os.path.join(repo, 'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits')
    # Pick first positional arg that's not a flag; default to the standard cube
    positional = [a for a in sys.argv[1:] if not a.startswith('-')]
    cube_path = positional[0] if positional else default_cube

    print('=' * 70)
    print(f'IFU spectral resolution diagnostic')
    print(f'  cube: {cube_path}')
    print('=' * 70)

    with fits.open(cube_path) as h:
        hdr = h[0].header
        cube_shape = h[0].data.shape

    crval = hdr['CRVAL3']
    cdelt = hdr['CD3_3']
    npix = hdr['NAXIS3']
    crpix = hdr.get('CRPIX3', 1.0)
    lam = crval + cdelt * (np.arange(npix) + 1 - crpix)

    print(f'\nCube shape: {cube_shape}  (NAXIS3 x NAXIS2 x NAXIS1)')
    print(f'Wavelength axis: CRVAL3={crval:.3f} A, CD3_3={cdelt:.3f} A/pix, '
          f'CRPIX3={crpix:.1f}, NAXIS3={npix}')
    print(f'  lambda range: {lam[0]:.1f} - {lam[-1]:.1f} A')
    print(f'  d_lambda    : {cdelt:.3f} A/pix  (constant in linear lambda)')

    # DISPSCAL: native KCWI keyword for the instrumental sigma in pixel units
    # (per the existing bootstrap_ppxf.py setup, FWHM_inst = 2.355 * DISPSCAL * d_lambda)
    if 'DISPSCAL' not in hdr:
        print('\nWARNING: DISPSCAL not in header — cannot derive instrumental FWHM.')
        return
    wdisp = float(hdr['DISPSCAL'])
    print(f'\nDISPSCAL header value: {wdisp:.4f}  (instrumental sigma in pixels per resolution element)')

    # Wavelength-dependent dlam_gal (the same expression bootstrap_ppxf uses,
    # via np.gradient in case the wavelength grid isn't perfectly linear).
    dlam = np.gradient(lam)
    fwhm_inst = 2.355 * wdisp * dlam       # Angstroms (FWHM)
    sigma_inst_AA = wdisp * dlam            # Angstroms (sigma)
    R = lam / fwhm_inst                     # resolving power
    sigma_inst_v = C_KMS * sigma_inst_AA / lam  # km/s (1-sigma)

    # Native log-lambda velocity step (= ppxf velscale at this wavelength)
    log_lam = np.log(lam)
    velscale_native = C_KMS * np.gradient(log_lam)  # km/s per native pixel

    # Per-band table
    print(f'\n{"lambda (A)":>10} {"FWHM (A)":>10} {"sigma (A)":>10} {"R":>7} '
          f'{"sigma_v (km/s)":>16} {"velscale (km/s/pix)":>22}')
    print('-' * 80)
    for j in [0, npix // 4, npix // 2, 3 * npix // 4, npix - 1]:
        print(f'{lam[j]:10.1f} {fwhm_inst[j]:10.3f} {sigma_inst_AA[j]:10.3f} '
              f'{R[j]:7.0f} {sigma_inst_v[j]:16.1f} {velscale_native[j]:22.2f}')

    # ppxf fitting band stats (6500-7500 A used by 03/03b/03c/06)
    band_mask = (lam >= 6500) & (lam <= 7500)
    if band_mask.any():
        lam_b = lam[band_mask]
        fwhm_b = fwhm_inst[band_mask]
        sig_b = sigma_inst_v[band_mask]
        vsc_b = velscale_native[band_mask]
        print(f'\nIn the ppxf fitting band 6500-7500 A ({band_mask.sum()} pixels):')
        print(f'  FWHM (A)         : min {fwhm_b.min():.3f}, max {fwhm_b.max():.3f}, '
              f'median {np.median(fwhm_b):.3f}')
        print(f'  R = lam/FWHM     : min {(lam_b/fwhm_b).min():.0f}, max {(lam_b/fwhm_b).max():.0f}, '
              f'median {np.median(lam_b/fwhm_b):.0f}')
        print(f'  sigma_v (km/s)   : min {sig_b.min():.1f}, max {sig_b.max():.1f}, '
              f'median {np.median(sig_b):.1f}')
        print(f'  velscale (km/s)  : median {np.median(vsc_b):.2f}  '
              f'(this is ppxf velscale at the band centroid)')

    # Comparison with z=0.676 rest-frame: the deflector has typical galaxy
    # sigma ~ 200 km/s — well above the instrumental sigma — so we are
    # resolved. Print the deconvolved-equivalent in rest-frame.
    print(f'\nAt z=0.67564 (rest-frame):')
    print(f'  Rest-frame lambda corresponding to obs 6500-7500 A : '
          f'{6500/(1+0.67564):.0f} - {7500/(1+0.67564):.0f} A')
    print(f'  Rest-frame instrumental sigma (band median, deconvolved by 1+z): '
          f'{np.median(sig_b)/(1+0.67564):.1f} km/s')

    # Optional plot
    if '--plot' in sys.argv:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
        axes[0].plot(lam, fwhm_inst, lw=1.2)
        axes[0].axvspan(6500, 7500, color='C2', alpha=0.18, label='ppxf band')
        axes[0].set_xlabel('Wavelength (A)'); axes[0].set_ylabel('FWHM (A)')
        axes[0].set_title('Instrumental FWHM vs lambda'); axes[0].legend(); axes[0].grid(alpha=0.3)
        axes[1].plot(lam, R, lw=1.2)
        axes[1].axvspan(6500, 7500, color='C2', alpha=0.18)
        axes[1].set_xlabel('Wavelength (A)'); axes[1].set_ylabel('R = lam / FWHM')
        axes[1].set_title('Resolving power'); axes[1].grid(alpha=0.3)
        axes[2].plot(lam, sigma_inst_v, lw=1.2, label='sigma_inst (km/s)')
        axes[2].plot(lam, velscale_native, lw=1.2, color='C1', label='velscale (km/s/pix)')
        axes[2].axvspan(6500, 7500, color='C2', alpha=0.18)
        axes[2].set_xlabel('Wavelength (A)'); axes[2].set_ylabel('km/s')
        axes[2].set_title('Velocity resolution'); axes[2].legend(); axes[2].grid(alpha=0.3)
        plt.tight_layout()
        figpath = os.path.join(repo, 'figures', 'ifu_spectral_resolution.png')
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        plt.savefig(figpath, dpi=160, bbox_inches='tight')
        print(f'\nSaved plot: {figpath}')


if __name__ == '__main__':
    main()
