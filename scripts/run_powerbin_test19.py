#!/usr/bin/env python
"""
Standalone runner for Test 19 (PowerBin) from 05x_sigma_discrepancy_diagnostic.
Sets up minimal state, runs PowerBin + ppxf bootstrap, saves results.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, uniform_filter1d
from importlib import resources
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from powerbin import PowerBin
import os, json

plt.rcParams['figure.facecolor'] = 'white'
plt.rc('font', family='serif', size=14)

# ── Load IFU cube ──
ifu_file = os.path.join(os.path.dirname(__file__), '..',
                        'Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits')
ifu_file = os.path.abspath(ifu_file)
print(f"Loading {ifu_file}")
with fits.open(ifu_file) as hdul:
    hdr = hdul[0].header
    cube = np.asarray(hdul[0].data, dtype=float)

crval = hdr['CRVAL3']
cdelt = hdr['CD3_3']
npix_w = hdr['NAXIS3']
crpix = hdr.get('CRPIX3', 1.0)
pix = np.arange(npix_w)
lam = crval + cdelt * (pix + 1 - crpix)
c_kms = 299792.458

print(f'Cube: {cube.shape}, wavelength: {lam[0]:.0f}-{lam[-1]:.0f} Å')

# ── Wide sub-image ──
y_wide = slice(35, 75)
x_wide = slice(30, 70)
cube_wide = cube[:, y_wide, x_wide]
wl_wide = np.sum(cube_wide, axis=0)
ny_w, nx_w = cube_wide.shape[1], cube_wide.shape[2]

# Peak (center of deflector)
wl_full = np.sum(cube, axis=0)
wl_s = gaussian_filter(wl_full, sigma=2)
cy_pk, cx_pk = np.unravel_index(np.argmax(wl_s), wl_s.shape)
dx_pk = cx_pk - x_wide.start
dy_pk = cy_pk - y_wide.start
print(f"Center: pixel ({cx_pk},{cy_pk}), sub-image ({dx_pk},{dy_pk})")

# Radial distances (arcsec)
wcs_ifu_2d = WCS(hdr, naxis=2)
yy_w, xx_w = np.mgrid[y_wide.start:y_wide.stop, x_wide.start:x_wide.stop]
ra_w, dec_w = wcs_ifu_2d.pixel_to_world_values(xx_w.ravel(), yy_w.ravel())
ra_c, dec_c = wcs_ifu_2d.pixel_to_world_values(cx_pk, cy_pk)
dra_w = (ra_w.reshape(ny_w, nx_w) - ra_c) * np.cos(np.radians(dec_c)) * 3600
ddec_w = (dec_w.reshape(ny_w, nx_w) - dec_c) * 3600
r_w = np.sqrt(dra_w**2 + ddec_w**2)

# Noise
noise_sky = np.std(cube[:, 28:40, 45:70], axis=(1, 2))

# ── ppxf helpers ──
sps_name_boot = "emiles"
N_BOOT = 100
degrees_boot = np.arange(15, 26)
ppxf_dir = resources.files('ppxf')
filename_fsps = ppxf_dir / 'sps_models' / 'spectra_fsps_9.0.npz'
filename_emiles = ppxf_dir / 'sps_models' / f'spectra_{sps_name_boot}_9.0.npz'


def quick_ppxf(flux_):
    """Run ppxf at degree=20 on a spectrum, return (sigma, V, chi2)."""
    noise = np.std(cube[:, 28:40, 45:70], axis=(1, 2))
    lam_t = crval + cdelt * (pix + 1 - crpix)
    m = (lam_t >= 6500) & (lam_t <= 7500)
    lam_t2, f_, n_ = lam_t[m], flux_[m], noise[m]
    ll = np.log(lam_t2)
    dll = (ll[-1] - ll[0]) / (len(lam_t2) - 1)
    lln = np.arange(ll[0], ll[-1] + dll, dll)
    f_ = np.interp(lln, ll, f_)
    n_ = np.interp(lln, ll, n_)
    lam_t2 = np.exp(lln)
    m2 = (lam_t2 >= 6000) & (lam_t2 <= 7500)
    lam_t2, f_, n_ = lam_t2[m2], f_[m2], n_[m2]
    gal = f_ / np.median(f_)
    lg = np.copy(lam_t2)
    lg *= np.median(util.vac_to_air(lg) / lg)
    vs = c_kms * (np.log(lg)[-1] - np.log(lg)[0]) / (len(lg) - 1)
    lr = lg / (1 + 0.67511)
    fr = 2.355 * hdr["DISPSCAL"] * np.gradient(lg) / (1 + 0.67511)
    sp = lib.sps_lib(filename_fsps, vs, {"lam": lr, "fwhm": fr}, lam_range=[3500, 5000])
    gp = util.determine_goodpixels(np.log(lr), [3500, 5000])
    pp = ppxf(sp.templates, gal, np.sqrt(n_**2), vs, [0, 300.],
              goodpixels=gp, plot=False, moments=2, trig=False,
              degree=20, lam=lr, lam_temp=sp.lam_temp, mdegree=0)
    return pp.sol[1], pp.sol[0], pp.chi2


def bootstrap_vrms_bin(flux_in, noise_in, n_bootstrap=100, seed=42):
    """Bootstrap V, sigma, V_rms for one radial bin spectrum."""
    z = 0.67511
    mask_wl = (lam >= 6500) & (lam <= 7500)
    lam_t = lam[mask_wl]
    f_t = flux_in[mask_wl]
    n_t = noise_in[mask_wl]
    ll = np.log(lam_t)
    dll = (ll[-1] - ll[0]) / (len(lam_t) - 1)
    lln = np.arange(ll[0], ll[-1] + dll, dll)
    f_t = np.interp(lln, ll, f_t)
    n_t = np.interp(lln, ll, n_t)
    lam_r = np.exp(lln)
    m2 = (lam_r >= 6000) & (lam_r <= 7500)
    lam_r, f_t, n_t = lam_r[m2], f_t[m2], n_t[m2]

    if np.median(f_t) <= 0:
        return None

    galaxy = f_t / np.median(f_t)
    noise_norm = np.sqrt(n_t**2) / np.median(f_t)
    lg = np.copy(lam_r)
    lg *= np.median(util.vac_to_air(lg) / lg)
    vs = c_kms * (np.log(lg)[-1] - np.log(lg)[0]) / (len(lg) - 1)
    lr = lg / (1 + z)
    fr = 2.355 * hdr["DISPSCAL"] * np.gradient(lg) / (1 + z)

    fn = ppxf_dir / "sps_models" / f"spectra_{sps_name_boot}_9.0.npz"
    sp = lib.sps_lib(fn, vs, {"lam": lr, "fwhm": fr}, lam_range=[3500, 5000])
    gp = util.determine_goodpixels(np.log(lr), [3500, 5000])

    # Initial fit at median degree
    deg_mid = int(np.median(degrees_boot))
    pp0 = ppxf(sp.templates, galaxy, noise_norm, vs, [0, 300.],
               goodpixels=gp, plot=False, moments=2, trig=False,
               degree=deg_mid, lam=lr, lam_temp=sp.lam_temp, mdegree=0)
    bestfit = pp0.bestfit
    resid = galaxy - bestfit

    # Local residual scaling
    g_std = np.std(resid)
    if g_std > 0:
        loc_mean = uniform_filter1d(resid, size=75, mode="reflect")
        loc_sq = uniform_filter1d(resid**2, size=75, mode="reflect")
        loc_var = np.maximum(loc_sq - loc_mean**2, 0)
        scale = np.clip(np.sqrt(loc_var) / g_std, 0.2, 5.0)
    else:
        scale = np.ones_like(resid)
    scaled_resid = resid * scale

    # Bootstrap
    rng = np.random.default_rng(seed)
    V_samp = np.full(n_bootstrap, np.nan)
    sig_samp = np.full(n_bootstrap, np.nan)
    vrms_samp = np.full(n_bootstrap, np.nan)

    for i in range(n_bootstrap):
        signs = rng.choice([-1, 1], size=len(galaxy))
        gal_boot = bestfit + scaled_resid * signs
        try:
            pp = ppxf(sp.templates, gal_boot, noise_norm, vs, [0, 300.],
                      goodpixels=gp, plot=False, moments=2, trig=False,
                      degree=deg_mid, lam=lr, lam_temp=sp.lam_temp, mdegree=0)
            V_samp[i] = pp.sol[0]
            sig_samp[i] = pp.sol[1]
            vrms_samp[i] = np.sqrt(pp.sol[0]**2 + pp.sol[1]**2)
        except Exception:
            pass

    valid_s = sig_samp[np.isfinite(sig_samp)]
    valid_v = V_samp[np.isfinite(V_samp)]
    valid_vr = vrms_samp[np.isfinite(vrms_samp)]

    if len(valid_s) < 5:
        return None

    return {
        "sig_p16": np.percentile(valid_s, 16),
        "sig_p50": np.percentile(valid_s, 50),
        "sig_p84": np.percentile(valid_s, 84),
        "V_p16": np.percentile(valid_v, 16),
        "V_p50": np.percentile(valid_v, 50),
        "V_p84": np.percentile(valid_v, 84),
        "vrms_p16": np.percentile(valid_vr, 16),
        "vrms_p50": np.percentile(valid_vr, 50),
        "vrms_p84": np.percentile(valid_vr, 84),
        "n_valid": len(valid_s),
    }


# ══════════════════════════════════════════════════
# TEST 19: PowerBin spatial binning
# ══════════════════════════════════════════════════
print("\n" + "="*60)
print("TEST 19: PowerBin spatial binning")
print("="*60)

wl_pb = wl_wide.copy()
ny_pb, nx_pb = ny_w, nx_w

# Select spaxels: positive signal, within R<4 arcsec of center
r_max_pb = 4.0
pb_mask = (wl_pb > 0) & (r_w <= r_max_pb)
n_pb_input = np.sum(pb_mask)
print(f'PowerBin input: {n_pb_input} spaxels (R<{r_max_pb} arcsec, positive signal)')

# Pixel coordinates
yy_pb, xx_pb = np.mgrid[0:ny_pb, 0:nx_pb]
x_pb = xx_pb[pb_mask].ravel()
y_pb = yy_pb[pb_mask].ravel()
xy_pb = np.column_stack([x_pb, y_pb])

# Per-channel S/N in the ppxf fitting range (6500-7500 Ang)
# This is the relevant S/N for spectral fitting quality
mask_fit_range = (lam >= 6500) & (lam <= 7500)
cube_fit_range = cube_wide[mask_fit_range, :, :]
noise_fit_median = np.median(noise_sky[mask_fit_range])
signal_per_chan = np.median(cube_fit_range, axis=0)  # per-spaxel median flux/channel
sn_per_chan = signal_per_chan / noise_fit_median

sn_input = sn_per_chan[pb_mask].ravel()
print(f'Per-channel S/N (6500-7500 A): median={np.median(sn_input):.2f}, '
      f'min={np.min(sn_input):.2f}, max={np.max(sn_input):.2f}')

# PowerBin capacity = S/N^2 (additive: co-adding N spaxels gives S/N = sqrt(sum(sn_i^2)))
capacity_pb = sn_input**2
target_sn_pb = 12  # gives ~20 bins from total capacity ~3000
target_cap = target_sn_pb**2

print(f'\nTarget S/N = {target_sn_pb} (capacity = {target_cap})')
print(f'Expected bins ~ {np.sum(capacity_pb) / target_cap:.0f}')

pb = PowerBin(xy_pb, capacity_pb, target_capacity=target_cap,
              pixelsize=1.0, verbose=1)

n_pbins = len(pb.npix)
print(f'\nPowerBin result: {n_pbins} bins from {n_pb_input} spaxels')
print(f'  Pixels per bin: median={np.median(pb.npix):.0f}, '
      f'min={np.min(pb.npix)}, max={np.max(pb.npix)}')
print(f'  Fractional S/N scatter: {pb.rms_frac:.1f}%')
print(f'  Single-pixel bins: {np.sum(pb.single)}')

# ── Build bin map ──
pb_map = np.full((ny_pb, nx_pb), -1, dtype=int)
for k in range(len(x_pb)):
    pb_map[y_pb[k], x_pb[k]] = pb.bin_num[k]

pb_r = np.zeros(n_pbins)
pb_n = np.zeros(n_pbins, dtype=int)
for b in range(n_pbins):
    mask_b = pb_map == b
    pb_r[b] = np.mean(r_w[mask_b])
    pb_n[b] = np.sum(mask_b)

# ── Visualization ──
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
delta = 14
y1, y2 = max(0, dy_pk - delta), min(ny_pb, dy_pk + delta)
x1, x2 = max(0, dx_pk - delta), min(nx_pb, dx_pk + delta)

axes[0].imshow(wl_pb[y1:y2, x1:x2], origin='lower', cmap='viridis')
axes[0].set_title('White-light')

axes[1].imshow(np.where(pb_map[y1:y2, x1:x2] >= 0,
                        pb_map[y1:y2, x1:x2], np.nan),
               origin='lower', cmap='tab20')
axes[1].set_title(f'PowerBin ({n_pbins} bins, S/N={target_sn_pb})')

im2 = axes[2].imshow(np.where(pb_map[y1:y2, x1:x2] >= 0,
                               r_w[y1:y2, x1:x2], np.nan),
                     origin='lower', cmap='RdYlBu_r')
axes[2].set_title('Radial distance (arcsec)')
plt.colorbar(im2, ax=axes[2], label='R (arcsec)')

for ax in axes:
    ax.plot(dx_pk - x1, dy_pk - y1, 'w+', ms=12, mew=2)

plt.suptitle(f'PowerBin: {n_pbins} bins from {n_pb_input} spaxels '
             f'(S/N={target_sn_pb}, R<{r_max_pb} arcsec)', fontsize=14)
plt.tight_layout()
figdir = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(figdir, exist_ok=True)
plt.savefig(os.path.join(figdir, 'test19_powerbin_map.png'), dpi=150, bbox_inches='tight')
print(f"\nSaved figures/test19_powerbin_map.png")
plt.close()

sn_achieved = np.sqrt(pb.bin_capacity)
print(f'S/N achieved per bin: median={np.median(sn_achieved):.1f}, '
      f'min={np.min(sn_achieved):.1f}, max={np.max(sn_achieved):.1f}')

# ── PowerBin diagnostic plot ──
fig_pb = plt.figure(figsize=(8, 6))
pb.plot(capacity_scale='sqrt', ylabel='S/N')
plt.savefig(os.path.join(figdir, 'test19_powerbin_diagnostic.png'), dpi=150, bbox_inches='tight')
plt.close()

# ── ppxf + bootstrap on each bin ──
print(f'\nRunning ppxf + bootstrap on {n_pbins} PowerBin bins...')
print(f'  EMILES templates, degrees 15-25, N_bootstrap={N_BOOT}')

pbin_results = []
for b in range(n_pbins):
    mask_b = pb_map == b
    n_spax = np.sum(mask_b)
    if n_spax < 1:
        continue

    flux_ = np.average(cube_wide[:, mask_b], axis=1)
    r_mean = pb_r[b]

    print(f'  Bin {b:3d}: R={r_mean:.2f}", {n_spax:3d} spax...', end=' ')

    sig, V, chi2 = quick_ppxf(flux_)
    boot = bootstrap_vrms_bin(flux_, noise_sky, n_bootstrap=N_BOOT, seed=700 + b)

    if boot is not None:
        result = {
            'bin_id': b, 'r_mean': r_mean, 'n_spax': n_spax,
            'sigma': sig, 'V': V, 'chi2': chi2,
            'vrms': np.sqrt(V**2 + sig**2),
        }
        result.update(boot)
        pbin_results.append(result)
        s = boot['sig_p50']
        sl = boot['sig_p50'] - boot['sig_p16']
        sh = boot['sig_p84'] - boot['sig_p50']
        print(f'sig={s:.0f} -{sl:.0f}/+{sh:.0f} km/s')
    else:
        print('bootstrap failed')

print(f'\nCompleted {len(pbin_results)} / {n_pbins} bins')

# ── Results plots ──
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

sigma_pbmap = np.full((ny_pb, nx_pb), np.nan)
V_pbmap = np.full((ny_pb, nx_pb), np.nan)
vrms_pbmap = np.full((ny_pb, nx_pb), np.nan)
for res in pbin_results:
    mask_b = pb_map == res['bin_id']
    sigma_pbmap[mask_b] = res['sig_p50']
    V_pbmap[mask_b] = res['V_p50']
    vrms_pbmap[mask_b] = res['vrms_p50']

im0 = axes[0, 0].imshow(V_pbmap[y1:y2, x1:x2], origin='lower',
                         cmap='RdBu_r', vmin=-100, vmax=100)
axes[0, 0].set_title('V (km/s) — PowerBin')
plt.colorbar(im0, ax=axes[0, 0])

im1 = axes[0, 1].imshow(sigma_pbmap[y1:y2, x1:x2], origin='lower',
                         cmap='plasma', vmin=100, vmax=400)
axes[0, 1].set_title(r'$\sigma$ (km/s) — PowerBin')
plt.colorbar(im1, ax=axes[0, 1])

im2 = axes[0, 2].imshow(vrms_pbmap[y1:y2, x1:x2], origin='lower',
                         cmap='inferno', vmin=100, vmax=400)
axes[0, 2].set_title(r'$V_{\rm rms}$ (km/s) — PowerBin')
plt.colorbar(im2, ax=axes[0, 2])

for ax in axes[0]:
    ax.contour(wl_pb[y1:y2, x1:x2], levels=5, colors='white', alpha=0.4, linewidths=0.5)
    ax.plot(dx_pk - x1, dy_pk - y1, 'w+', ms=10, mew=2)

r_plot = [res['r_mean'] for res in pbin_results]
sig_plot = [res['sig_p50'] for res in pbin_results]
sig_lo = [res['sig_p50'] - res['sig_p16'] for res in pbin_results]
sig_hi = [res['sig_p84'] - res['sig_p50'] for res in pbin_results]
axes[1, 0].errorbar(r_plot, sig_plot, yerr=[sig_lo, sig_hi], fmt='s',
                    color='C1', capsize=3, markersize=6, label='PowerBin')
axes[1, 0].axhline(210, color='gray', ls='--', lw=0.8, label=r'$\sigma=210$')
axes[1, 0].set_xlabel('R (arcsec)')
axes[1, 0].set_ylabel(r'$\sigma$ (km/s)')
axes[1, 0].set_title(r'$\sigma(R)$ — PowerBin')
axes[1, 0].legend(fontsize=10)
axes[1, 0].grid(alpha=0.3)

V_plot = [res['V_p50'] for res in pbin_results]
V_lo = [res['V_p50'] - res['V_p16'] for res in pbin_results]
V_hi = [res['V_p84'] - res['V_p50'] for res in pbin_results]
axes[1, 1].errorbar(r_plot, V_plot, yerr=[V_lo, V_hi], fmt='s',
                    color='C1', capsize=3, markersize=6)
axes[1, 1].axhline(0, color='gray', ls='--', lw=0.8)
axes[1, 1].set_xlabel('R (arcsec)')
axes[1, 1].set_ylabel('V (km/s)')
axes[1, 1].set_title('V(R) — PowerBin')
axes[1, 1].grid(alpha=0.3)

vrms_plot = [res['vrms_p50'] for res in pbin_results]
vrms_lo = [res['vrms_p50'] - res['vrms_p16'] for res in pbin_results]
vrms_hi = [res['vrms_p84'] - res['vrms_p50'] for res in pbin_results]
axes[1, 2].errorbar(r_plot, vrms_plot, yerr=[vrms_lo, vrms_hi], fmt='s',
                    color='C1', capsize=3, markersize=6)
axes[1, 2].axhline(210, color='gray', ls='--', lw=0.8)
axes[1, 2].set_xlabel('R (arcsec)')
axes[1, 2].set_ylabel(r'$V_{\rm rms}$ (km/s)')
axes[1, 2].set_title(r'$V_{\rm rms}(R)$ — PowerBin')
axes[1, 2].grid(alpha=0.3)

plt.suptitle(f'PowerBin kinematics (S/N={target_sn_pb}, EMILES, N_boot={N_BOOT})', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(figdir, 'test19_powerbin_results.png'), dpi=150, bbox_inches='tight')
print(f"\nSaved figures/test19_powerbin_results.png")
plt.close()

# ── Summary table ──
print(f'\n{"Bin":>4} {"R":>6} {"N":>4} {"sigma":>18} {"V":>18} {"V_rms":>18}')
print('-' * 72)
for res in sorted(pbin_results, key=lambda x: x['r_mean']):
    sl = res['sig_p50'] - res['sig_p16']
    sh = res['sig_p84'] - res['sig_p50']
    vl = res['V_p50'] - res['V_p16']
    vh = res['V_p84'] - res['V_p50']
    vrl = res['vrms_p50'] - res['vrms_p16']
    vrh = res['vrms_p84'] - res['vrms_p50']
    print(f'{res["bin_id"]:4d} {res["r_mean"]:6.2f} {res["n_spax"]:4d} '
          f'{res["sig_p50"]:6.0f} -{sl:4.0f}/+{sh:4.0f} '
          f'{res["V_p50"]:6.0f} -{vl:4.0f}/+{vh:4.0f} '
          f'{res["vrms_p50"]:6.0f} -{vrl:4.0f}/+{vrh:4.0f}')

# ── Save results ──
resdir = os.path.join(os.path.dirname(__file__), '..', 'results')
os.makedirs(resdir, exist_ok=True)
np.savez(os.path.join(resdir, 'test19_powerbin_results.npz'),
         pbin_results=pbin_results,
         n_pbins=n_pbins, target_sn=target_sn_pb, r_max=r_max_pb,
         N_BOOT=N_BOOT, sps_name=sps_name_boot,
         pb_npix=pb.npix, pb_bin_capacity=pb.bin_capacity,
         pb_r=pb_r, pb_n=pb_n, pb_map=pb_map)
print(f"\nSaved results/test19_powerbin_results.npz")
print("\nDone!")
