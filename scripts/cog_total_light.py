"""Empirical curve-of-growth (CoG) total light vs the single-Sérsic total (2026-06-16).

Motivated by the cored/cD investigation (notebook 18, Test A): the single-Sérsic total is box-limited
and non-convergent for this cD, so the headline M⋆ (built on a small-box Sérsic) is a LOWER BOUND.
This measures the EMPIRICAL enclosed light directly — cumulative flux in growing elliptical apertures
(masked, sky-subtracted, validated-Sérsic-filled inside the mask) — and compares to the single-Sérsic
total. The CoG−Sérsic gap at a fixed metric radius is a "total-light method" systematic on M⋆.

Sky is estimated from the far outskirts (sigma-clipped); for a cD the corners still hold envelope/ICL
light, so the sky is if anything OVER-estimated → the CoG excess is a conservative LOWER bound.

Run:  python -m scripts.cog_total_light
Outputs: results/cog_total_light.npz, results/figures/cog_total_light.png
"""
import numpy as np
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from scipy.special import gammaincinv, gamma
from scripts.mask_method_comparison import load_band

ORDER = ['F200LP', 'F140W', 'F150W2', 'F322W2']
EXP = {'F200LP':'../velocity_dispersion_from_IFU/AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits',
       'F140W':'../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3_mask.fits',
       'F150W2':'../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f150w2_i2d_mask.fits',
       'F322W2':'../velocity_dispersion_from_IFU/jw05594-o101_t103_nircam_clear-f322w2_i2d_mask.fits'}
# JWST are huge mosaics -> restrict the analysis window so field sources don't swamp the CoG.
RMAX = {'F200LP':11.0, 'F140W':11.0, 'F150W2':8.0, 'F322W2':8.0}
_TAB = np.load('results/sersic_parameter_table.npz', allow_pickle=True)
def tab(b, k): return float({k2: v for k2, v in _TAB[f'{b}_cen']}[k])


def cog_band(band):
    b = load_band(band); ps, zp = b['pix_scale'], b['ab_zp']
    cx, cy = b['cx'], b['cy']; img = np.nan_to_num(b['img'].astype(float))
    em = fits.getdata(EXP[band]).astype(bool)
    ny, nx = img.shape; yy, xx = np.mgrid[:ny, :nx]
    n = tab(band, 'n'); reff_px = tab(band, 'reff_as')/ps; ell = tab(band, '_ell'); pa = tab(band, 'pa')
    ct, st = np.cos(-np.radians(pa)), np.sin(-np.radians(pa))
    xr = (xx-cx)*ct - (yy-cy)*st; yr = (xx-cx)*st + (yy-cy)*ct
    rell = np.hypot(xr, yr/max(1-ell, 0.2)) * ps
    rmax = min(RMAX[band], min(ny, nx)/2*ps - 0.5)
    # sky from the far ring inside the window, sigma-clipped
    skreg = (rell > 0.85*rmax) & ~em & np.isfinite(b['img'])
    s = img[skreg]
    for _ in range(5):
        m, sd = np.median(s), np.std(s); s = s[np.abs(s-m) < 3*sd]
    sky = float(np.median(s)); data = img - sky
    # validated-shape model, amplitude+sky fit on unmasked r<6"
    shape = np.clip(np.asarray(Sersic2D(amplitude=1., r_eff=reff_px, n=n, x_0=cx, y_0=cy,
                    ellip=ell, theta=np.radians(pa))(xx, yy), float), 0, None)
    sel = (rell < 6) & ~em & np.isfinite(b['img'])
    A = np.vstack([shape[sel], np.ones(sel.sum())]).T
    (amp, _), *_ = np.linalg.lstsq(A, img[sel], rcond=None)
    filled = np.where(em, amp*shape, data)
    aa = np.arange(0.5, rmax+0.01, 0.25)
    F = np.array([filled[rell < a].sum() for a in aa])
    mcum = -2.5*np.log10(np.clip(F, 1e-9, None)) + zp
    bn = gammaincinv(2*n, 0.5)
    Fser = 2*np.pi*amp*reff_px**2*n*np.exp(bn)/bn**(2*n)*gamma(2*n)*(1-ell)
    mser = -2.5*np.log10(Fser) + zp
    # sky-sensitivity: redo with sky +/- its scatter
    sksd = float(np.std(s))
    def cog_at(a, dsky):
        f = np.where(em, amp*shape, (img - (sky+dsky)))
        return -2.5*np.log10(max(f[rell < a].sum(), 1e-9)) + zp
    a8 = 8.0 if rmax >= 8 else rmax
    m8, m8p, m8m = cog_at(a8, 0), cog_at(a8, +sksd), cog_at(a8, -sksd)
    return dict(band=band, aa=aa, mcum=mcum, mser=mser, n=n, sky=sky, sksd=sksd,
                a8=a8, m8=m8, m8_skyhi=m8p, m8_skylo=m8m, amp=amp)


def main():
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    res = {n: cog_band(n) for n in ORDER}
    out = {}
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    print('%-7s %8s %8s %9s %8s' % ('band', 'Sersic', 'CoG@a', 'a(")', 'ΔCoG-Ser'))
    for n in ORDER:
        r = res[n]; ax.plot(r['aa'], r['mcum'], '-', label='%s CoG' % n)
        ax.axhline(r['mser'], ls=':', color=ax.lines[-1].get_color(), alpha=0.6)
        print('%-7s %8.2f %8.2f %9.1f %+8.2f  (sky±: %+.2f/%+.2f)' % (
            n, r['mser'], r['m8'], r['a8'], r['m8']-r['mser'], r['m8_skyhi']-r['m8'], r['m8_skylo']-r['m8']))
        for k in ['aa', 'mcum', 'mser', 'm8', 'a8', 'm8_skyhi', 'm8_skylo', 'n', 'sky']:
            out['%s_%s' % (n, k)] = r[k]
    ax.set_xlabel('elliptical semi-major a (arcsec)'); ax.set_ylabel('cumulative mag (AB)')
    ax.invert_yaxis(); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title('Empirical curve-of-growth vs single-Sérsic total (dotted)')
    fig.savefig('results/figures/cog_total_light.png', dpi=130, bbox_inches='tight')

    # Bagpipes M* on the CoG@8" mags (quiescent prior) -> the empirical total-light M*, so the
    # CoG-vs-single-Sérsic gap can be carried as a one-sided 'total-light/envelope' systematic.
    # (DECISION 2026-06-16: do NOT pick estimators by relation-consistency; this is a measurement-
    # quality systematic spanning a model-dependent (Sérsic apcorr) vs empirical-but-bandpass-dependent
    # (CoG) total. Central stays single-Sérsic; CoG carried as the one-sided up term.)
    import os
    import bagpipes as pipes
    from scripts.bagpipes_sersic_refit import FILT_LIST, mag_to_flam
    QUI = {'redshift':(0.674,0.676),
           'exponential':{'age':(4.,15.),'tau':(0.1,1.5),'massformed':(1.,15.),'metallicity':(0.,2.5)},
           'dust':{'type':'Calzetti','Av':(0.,0.6)}}
    qhead = np.load('results/mstar_headline_quiescent.npz', allow_pickle=True)
    piv = qhead['pivot']; mags8 = np.array([res[n]['m8'] for n in ORDER])
    flam = mag_to_flam(mags8, piv); flam_err = 0.10*flam
    g = pipes.galaxy('AGEL0206_COG8', lambda ID: np.array([flam, flam_err]).T, spectrum_exists=False,
                     filt_list=[os.path.abspath(p) for p in FILT_LIST], phot_units='ergscma')
    f = pipes.fit(g, QUI, run='AGEL0206_COG8')
    try: f.fit(verbose=False, sampler='multinest', n_live=400)
    except (AttributeError, OSError): f.fit(verbose=False, sampler='nautilus', n_live=400)
    f.posterior.get_advanced_quantities()
    p = np.percentile(f.posterior.samples['stellar_mass'], [16,50,84])   # H0=70; registry adds +0.0282
    out['logM_cog8_H070'] = p
    print('\nlog M*(CoG@8", quiescent, H0=70) = %.3f -> Planck %.3f' % (p[1], p[1]+0.0282))
    np.savez('results/cog_total_light.npz', **out)
    print('\nΔlogM* if CoG@a adopted (≈ -0.4·ΔCoG-Ser, same M/L):')
    for n in ORDER:
        r = res[n]; print('  %-7s %+.3f dex' % (n, -0.4*(r['m8']-r['mser'])))
    print('Saved -> results/cog_total_light.npz, results/figures/cog_total_light.png')


if __name__ == '__main__':
    main()
