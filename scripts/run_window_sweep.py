"""Standalone driver for the nb09a wavelength-window σ_e sweep.

Reproduces the cell `03f91c7a` of `notebooks/09a_wavelength_window_sweep.ipynb`
(without nbconvert overhead). Writes the same per-(window, sps) cache files
under `results/nb09a_wavelength_sweep/{wlabel}_{sps}_T{T0}_{T1}_N{N}.npz`.

Run:
    conda run -n ISMgas python scripts/run_window_sweep.py 2>&1 | tee /tmp/sweep.log

Once the cache directory has all caches at the desired N value, re-execute
the notebook (which loads the caches) to refresh the figures and tables.
"""
from __future__ import annotations

import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from scripts.final_sigma_e import (
    load_setup, extract_aperture_spectrum, Z_SYSTEMIC, BOOT_SEED, WINDOW,
)
from scripts.bootstrap_ppxf import (
    setup_ppxf_inputs_from_spectrum, run_bootstrap_single_degree, C_KMS,
)
from ppxf.ppxf import ppxf

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Match the notebook's WINDOWS list exactly.
_R3800_OBS = 3800.0 * (1 + Z_SYSTEMIC)
_R4000_OBS = 4000.0 * (1 + Z_SYSTEMIC)
_R5400_OBS = 5400.0 * (1 + Z_SYSTEMIC)

_R3500_OBS = 3500.0 * (1 + Z_SYSTEMIC)
_R3700_OBS = 3700.0 * (1 + Z_SYSTEMIC)

# Default windows: NO arc mask. Schema: (label, (lo, hi), arc_mask_flag)
WINDOWS = [
    ('w6500_7500', (6500.0, 7500.0), False),
    ('w6000_7500', (6000.0, 7500.0), False),
    ('w5500_7500', (5500.0, 7500.0), False),
    ('w5000_7500', (5000.0, 7500.0), False),
    ('w6500_8000', (6500.0, 8000.0), False),
    ('w6500_8500', (6500.0, 8500.0), False),
    ('w6000_8000', (6000.0, 8000.0), False),
    ('w5500_8000', (5500.0, 8000.0), False),
    ('w6000_8500', (6000.0, 8500.0), False),
    ('w5500_8500', (5500.0, 8500.0), False),
    ('wR4000_5400', (_R4000_OBS, _R5400_OBS), False),
    ('wR3800_5400', (_R3800_OBS, _R5400_OBS), False),
    # Arc-masked variants — same lam range, but with z=1.302 source-emission
    # bands explicitly excluded from goodpixels. Cache filenames carry an
    # `arcmask_` infix so they don't collide with the unmasked caches.
    ('wR3800_5400_arcmask', (_R3800_OBS, _R5400_OBS), True),
    ('wR3700_5400_arcmask', (_R3700_OBS, _R5400_OBS), True),
    ('wR3500_5400_arcmask', (_R3500_OBS, _R5400_OBS), True),
    # Full arc-masked counterparts of every base obs-frame and rest-frame
    # window. Some of the bluer/narrower windows have no arc-emission band
    # falling in their fit range — those will give bit-identical σ to the
    # unmasked variant (validates the masking machinery).
    ('w6500_7500_arcmask', (6500.0, 7500.0), True),
    ('w6000_7500_arcmask', (6000.0, 7500.0), True),
    ('w5500_7500_arcmask', (5500.0, 7500.0), True),
    ('w5000_7500_arcmask', (5000.0, 7500.0), True),
    ('w6500_8000_arcmask', (6500.0, 8000.0), True),
    ('w6500_8500_arcmask', (6500.0, 8500.0), True),
    ('w6000_8000_arcmask', (6000.0, 8000.0), True),
    ('w5500_8000_arcmask', (5500.0, 8000.0), True),
    ('w6000_8500_arcmask', (6000.0, 8500.0), True),
    ('w5500_8500_arcmask', (5500.0, 8500.0), True),
    ('wR4000_5400_arcmask', (_R4000_OBS, _R5400_OBS), True),
]
SPS_LIBS = ('fsps', 'emiles', 'xsl')
DEGREES = np.arange(15, 30)
N_BOOT = 200

SPS_REST_MIN = {'fsps': 100.0, 'emiles': 1680.0, 'xsl': 3500.0}
PADDING = 1.05

def sps_safe_obs_min(sps):
    return SPS_REST_MIN[sps] * PADDING * (1 + Z_SYSTEMIC)

SPS_LAM_RANGE_TEMP = {
    'fsps':   (3000.0, 6000.0),
    'emiles': (3000.0, 6000.0),
    'xsl':    (3500.0, 6000.0),
}

CACHE_DIR = os.path.join(REPO_ROOT, 'results', 'nb09a_wavelength_sweep')
os.makedirs(CACHE_DIR, exist_ok=True)


# Spectral masks for the deflector ppxf fit (rest-frame Å at z_lens=0.67564).
# Three are z=1.302 source-emission contamination; one is a telluric residual.
# Used when arc_mask=True.
ARC_MASKS_REST = [
    (3835.0, 3855.0),  # Mg II 2796/2803 from source (z=1.302)
    (4525.0, 4545.0),  # O2 A-band leading-edge telluric residual at obs 7593-7626 Å
                       # (NOT source emission — see NOTES_4534A_spike_investigation_2026-05-18.md)
    (5115.0, 5135.0),  # [O II] 3727/3729 from source — sits ON Mg b
    (5237.0, 5253.0),  # He I 3819.60 from source (z=1.302) → def-rest 5247.4 Å
                       # Identified 2026-05-27: 3-pixel +ve residual cluster at
                       # def-rest 5244-5248 Å (source-rest 3818-3820 Å) consistent
                       # across NEW _mtwdo_ AND OLD reductions → astrophysical, not CR.
                       # Mask spans ±8 Å around the line center to cover any
                       # source kinematic broadening + adjacent CR spikes at
                       # 5240 Å (NEW) and 5238.6 Å (OLD).
    (5260.0, 5340.0),  # Mg b / [Ne III] 3869 cluster
]

# Bad-pixel mask: local-MAD-flagged cosmic-ray / telluric residual pixels in
# the wide arc-masked window (rest 3800-5400 Å). Identified 2026-05-26 by
# rolling 75-pixel local-MAD on the canonical ppxf-fit residuals at
# |residual|/local_MAD > 3σ. 52 outlier pixels in 26 contiguous ranges
# (total ≈ 46 Å padded). Verified intrinsic to the data (not reduction-
# specific): same pixels flagged in both _mtwdo_ and original reductions
# (M6, NOTES_nb09e_audit_2026-05-20.md addendum). Biggest outlier is a
# 6-pixel cluster at def-rest 5232-5236 Å (obs 8768-8774 Å) at 26σ in
# local-MAD units. Each range padded ±0.5 Å to robustly catch the same
# pixels after log-rebinning. Used when bad_pix_mask=True.
BAD_PIXELS_REST = [
    # Original (2026-05-26) 26 cosmic-ray / sky residual entries identified by
    # local-MAD on ppxf-fit residuals at |residual|/local_MAD > 3σ on the
    # canonical wide-window fit. Biggest is the 6-pix cluster at 5232-5236.
    (4010.5, 4011.5),  (4167.4, 4168.4),
    # 2026-05-27 (M10 sky-line audit): MERGED original (4380.5, 4381.5)
    # with adjacent unmasked sky band 4382-4384 → (4380.0, 4384.0). Sky
    # noise 3.7× median at obs 7343-7345.
    (4380.0, 4384.0),
    (4400.2, 4401.7),  (4489.3, 4490.3),
    # M10 additions (2026-05-27): unmasked sky-line audit across full fit
    # window 3800-5400. Added 9 OH airglow / sky residual bands at >2.5×
    # median sky std, all VERIFIED non-source via the arc spectrum
    # (AGEL0206_arc_source_spectrum.pdf — source-rest at each band is flat).
    # Tag: # M10-sky
    (4553.0, 4554.6),  # M10-sky: obs 7629-7632 (2.6× sky)
    (4602.0, 4610.0),  # M10-sky: obs 7712-7724 (2.9× sky) — user-flagged 4600 Å structure
    (4625.1, 4629.8),  (4632.5, 4633.5),  (4650.2, 4651.2),
    (4652.7, 4654.3),  (4665.6, 4667.2),  (4669.3, 4670.3),
    (4687.3, 4688.3),  # M10-sky: obs 7854-7856 (2.5× sky)
    (4689.6, 4690.6),
    (4691.5, 4693.6),  # M10-sky: MERGED (4691.5, 4692.5) + adjacent obs 7862-7865 (3.5× sky)
    (4721.9, 4724.2),  (4725.6, 4726.6),  (4728.8, 4730.4),
    (4754.4, 4756.0),  (4768.9, 4769.9),
    (4770.8, 4771.8),  # M10-sky: obs 7994-7996 (2.7× sky)
    (4783.3, 4784.3),  (4789.0, 4791.9),
    (4941.8, 4942.8),
    (4951.6, 4955.0),  # M10-sky: obs 8297-8303 (3.8× sky)
    (4979.7, 4984.0),  (4986.3, 4988.6),
    (5011.9, 5015.3),  # M10-sky: obs 8398-8404 (4.0× sky)
    (5029.8, 5033.8),  # M10-sky: obs 8428-8435 (3.1× sky)
    (5052.5, 5054.2),
    (5228.8, 5231.2),  (5232.2, 5236.7),
]

def _apply_arc_mask(goodpixels, lam_rest, masks=ARC_MASKS_REST):
    bad = np.ones(len(lam_rest), dtype=bool); bad[goodpixels] = False
    for lo, hi in masks:
        bad |= ((lam_rest >= lo) & (lam_rest <= hi))
    return np.where(~bad)[0]

def _apply_bad_pixels_mask(goodpixels, lam_rest, masks=BAD_PIXELS_REST):
    """Mask the local-MAD-flagged bad pixels (cosmic-ray residuals).
    Same logic as _apply_arc_mask but for the BAD_PIXELS_REST list. Designed
    to be applied on top of _apply_arc_mask: call _apply_arc_mask first to
    drop the source-emission/telluric bands, then this to drop the
    individual cosmic-ray spikes."""
    bad = np.ones(len(lam_rest), dtype=bool); bad[goodpixels] = False
    for lo, hi in masks:
        bad |= ((lam_rest >= lo) & (lam_rest <= hi))
    return np.where(~bad)[0]


def _ppxf_one_window(flux, noise, hdr, sps, lam_obs_range, lam_fit_range,
                     n_boot, seed_offset, arc_mask=False):
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, hdr,
        sps_name=sps, z=Z_SYSTEMIC,
        lam_obs_range=lam_obs_range,
        lam_fit_range=lam_fit_range,
        lam_range_temp=SPS_LAM_RANGE_TEMP[sps],
        verbose=False, frame_galaxy='auto',
    )
    if arc_mask:
        gp = _apply_arc_mask(inputs['goodpixels'], inputs['lam_gal_rest'])
        inputs['goodpixels'] = gp
    n_pix = len(inputs['galaxy'])
    n_deg = len(DEGREES)
    V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
    V_boot = np.full((n_deg, n_boot), np.nan)
    sig_boot = np.full((n_deg, n_boot), np.nan)
    bf = np.zeros((n_deg, n_pix))

    for i, deg in enumerate(DEGREES):
        pp = ppxf(
            inputs['sps'].templates, inputs['galaxy'], inputs['noise'],
            inputs['velscale'], inputs['start'],
            goodpixels=inputs['goodpixels'], plot=False, moments=2, trig=False,
            degree=int(deg), mdegree=0,
            lam=inputs['lam_gal_rest'], lam_temp=inputs['lam_temp'], quiet=True,
        )
        V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
        bf[i] = pp.bestfit
        rb = run_bootstrap_single_degree(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_boot, window=WINDOW, seed=BOOT_SEED + seed_offset + i,
        )
        V_boot[i] = rb['V_samples']
        sig_boot[i] = rb['sigma_samples']
    return dict(
        sps=sps, lam_obs_range=lam_obs_range, lam_fit_range=lam_fit_range,
        n_pix=n_pix, lam_rest_min=float(inputs['lam_gal_rest'][0]),
        lam_rest_max=float(inputs['lam_gal_rest'][-1]),
        frame_galaxy=inputs['frame_galaxy'],
        V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
        V_boot=V_boot, sig_boot=sig_boot, bf=bf,
    )


def main(n_boot=None, only=None):
    """Run the wavelength-window sweep.

    n_boot: override the default N=200 bootstrap count (e.g. 500 for production)
    only:   iterable of window labels to fit; if None, run the full WINDOWS list
    """
    nb = N_BOOT if n_boot is None else int(n_boot)
    filt = set(only) if only is not None else None
    windows_to_fit = [w for w in WINDOWS if filt is None or w[0] in filt]
    if filt is not None:
        missing = filt - {w[0] for w in WINDOWS}
        if missing:
            raise SystemExit(f"--only includes labels not in WINDOWS: {sorted(missing)}")

    t_global = time.time()
    print(f"Cache directory: {CACHE_DIR}")
    print(f"N_BOOT = {nb}, fitting {len(windows_to_fit)} windows × {len(SPS_LIBS)} SPS", flush=True)
    print(f"Loading setup …", flush=True)
    state = load_setup()
    R_E = state['R_E']
    flux, noise, n_kept, sn_band = extract_aperture_spectrum(
        state, r_max=R_E, mask_weight=0.0
    )
    n_med = float(np.nanmedian(noise[noise > 0]))
    n_floor = n_med * 0.1
    n_zero = int(np.sum(~(noise > 0)))
    noise = np.where(noise > 0, noise, n_floor)
    print(f"Aperture R<R_e ({R_E:.3f}\"): {n_kept} active spaxels, "
          f"S/N = {sn_band:.1f}, noise floor → {n_zero} pix", flush=True)

    n_done = 0
    n_skip = 0
    n_total = len(windows_to_fit) * len(SPS_LIBS)
    for w_idx, (wlabel, lr, arc_mask) in enumerate(windows_to_fit):
        for sps_idx, sps in enumerate(SPS_LIBS):
            lr_eff = (max(lr[0], sps_safe_obs_min(sps)), lr[1])
            T0, T1 = int(SPS_LAM_RANGE_TEMP[sps][0]), int(SPS_LAM_RANGE_TEMP[sps][1])
            cache = os.path.join(CACHE_DIR, f'{wlabel}_{sps}_T{T0}_{T1}_N{nb}.npz')
            if os.path.exists(cache):
                d = np.load(cache, allow_pickle=True)
                print(f"  [skip] {wlabel}/{sps}  σ_orig={d['sig_orig'].min():.0f}-{d['sig_orig'].max():.0f}", flush=True)
                n_skip += 1
                continue
            t0 = time.time()
            d = _ppxf_one_window(flux, noise, state['hdr'], sps,
                                 lam_obs_range=lr_eff, lam_fit_range=lr_eff,
                                 n_boot=nb,
                                 seed_offset=(50_000 + 10_000 * w_idx + 100 * sps_idx),
                                 arc_mask=arc_mask)
            np.savez(cache, **d)
            n_done += 1
            elapsed = time.time() - t0
            clip = '' if lr_eff[0] == lr[0] else f' (CLIP→{lr_eff[0]:.0f})'
            print(f"  [done] {wlabel}/{sps}  σ={d['sig_orig'].min():.0f}-{d['sig_orig'].max():.0f}  "
                  f"V={d['V_orig'].mean():+.0f}  npix={d['n_pix']}  "
                  f"rest=[{d['lam_rest_min']:.0f},{d['lam_rest_max']:.0f}]  "
                  f"({elapsed:.0f}s){clip}  [{n_skip+n_done}/{n_total}]", flush=True)

    elapsed_global = time.time() - t_global
    print(f"\nTotal: {n_done} computed, {n_skip} from cache, {n_total} total — {elapsed_global:.0f}s ({elapsed_global/60:.1f} min)")


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description="ppxf σ_e sweep over fit-wavelength windows")
    p.add_argument('--n_boot', type=int, default=None,
                   help=f"Bootstrap count (default {N_BOOT}).")
    p.add_argument('--only', nargs='+', default=None,
                   help="Restrict to specific window labels (e.g. wR3800_5400_arcmask).")
    args = p.parse_args()
    main(n_boot=args.n_boot, only=args.only)
