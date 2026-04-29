"""Soft-mask (arc-spaxel I-weight × 0.5) σ_e track at N=500.

Tests the interpolation between hard-mask headline (Track A, 267.8 km/s) and
no-mask sensitivity (Track B, 252.8 km/s) by KEEPING arc-flagged spaxels in
the aperture but DOWN-WEIGHTING their contribution to the I-weighted spectrum
extraction by a factor of 0.5.

Caches to results/final_sigma_e_paper/{label}_{sps}_N500_softmask_w0p5.npz
matching the schema of `scripts.final_sigma_e.run_aperture_sps`.

Usage:
    python scripts/soft_mask_track.py
"""
import os
import sys
import warnings
from pathlib import Path
from time import perf_counter as clock

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)
warnings.filterwarnings("ignore", category=RuntimeWarning)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.final_sigma_e import (   # noqa: E402
    load_setup, CACHE_DIR, APERTURE_LABELS, APERTURE_FRACS, SPS_LIBS,
    Z_SYSTEMIC, DEGREES, BOOT_SEED, WINDOW,
)
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402
from scripts.bootstrap_ppxf_parallel import run_bootstrap_single_degree_parallel  # noqa: E402

N_BOOTSTRAP = 500
N_JOBS = 8
MASK_WEIGHT = 0.5
SUFFIX = "_softmask_w0p5"


def extract_aperture_spectrum_soft(state, r_max, mask_weight):
    """Same as final_sigma_e.extract_aperture_spectrum but with arc-spaxel
    I-weight × `mask_weight` instead of hard-drop."""
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


def main():
    state = load_setup()
    print(f"\nSoft mask track: arc-spaxel I-weight × {MASK_WEIGHT}")
    print(f"  N_bootstrap = {N_BOOTSTRAP}, parallel n_jobs = {N_JOBS}")
    print(f"  Suffix     = {SUFFIX}")

    for label_idx, (frac, label) in enumerate(zip(APERTURE_FRACS, APERTURE_LABELS)):
        r_max = frac * state["R_E"]
        flux, noise, n_active, sn_band = extract_aperture_spectrum_soft(
            state, r_max, MASK_WEIGHT,
        )
        print(f"\n  Aperture {label}  R<{r_max:.3f}\"  "
              f"({n_active} active spaxels, S/N={sn_band:.1f})")

        for sps_name in SPS_LIBS:
            cache = CACHE_DIR / f"{label}_{sps_name}_N{N_BOOTSTRAP}{SUFFIX}.npz"
            if cache.exists():
                print(f"    Skip {sps_name}: {cache.name} exists")
                continue
            print(f"    === {label}/{sps_name}/softmask_w{MASK_WEIGHT} ===")
            inputs = setup_ppxf_inputs_from_spectrum(
                flux, noise, state["hdr"], sps_name=sps_name, z=Z_SYSTEMIC,
                verbose=False, frame_galaxy="auto",
            )
            n_pix = len(inputs["galaxy"])
            n_deg = len(DEGREES)
            V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
            V_boot = np.full((n_deg, N_BOOTSTRAP), np.nan)
            sig_boot = np.full((n_deg, N_BOOTSTRAP), np.nan)
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
                rb = run_bootstrap_single_degree_parallel(
                    inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
                    n_bootstrap=N_BOOTSTRAP, window=WINDOW,
                    seed=BOOT_SEED + 2000 + 1000 * label_idx + i, n_jobs=N_JOBS,
                )
                V_boot[i] = rb["V_samples"]
                sig_boot[i] = rb["sigma_samples"]
            elapsed = clock() - t0
            np.savez(
                cache,
                label=label, sps_name=sps_name, r_max=float(r_max),
                mask_on=False, mask_weight=float(MASK_WEIGHT),
                n_spax=int(n_active), sn_band=float(sn_band),
                degrees=np.asarray(DEGREES),
                V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
                V_boot=V_boot, sig_boot=sig_boot,
                frame_galaxy=inputs["frame_galaxy"],
                velscale=float(inputs["velscale"]),
                n_pix=int(n_pix), n_bootstrap=int(N_BOOTSTRAP),
                elapsed_s=float(elapsed),
                galaxy=inputs["galaxy"], noise=inputs["noise"],
                lam_gal_rest=inputs["lam_gal_rest"], best_fit=bf_arr,
            )
            print(f"      done in {elapsed:.1f}s ({elapsed/60:.1f} min); "
                  f"σ={sig_orig.min():.0f}-{sig_orig.max():.0f} km/s")

    print("\n  All soft-mask N=500 fits committed to cache.")


if __name__ == "__main__":
    main()
