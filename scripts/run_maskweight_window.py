"""F200 spatial mask-weight sweep at a configurable fit window + spectral arc mask.

The headline narrow-window budget includes a ±16 km/s "mask systematic"
measured by comparing mask_weight=0.0 (hard-drop arc spaxels) vs
mask_weight=1.0 (keep all spaxels) at the headline 6500-7500 obs window.

At the wide arc-masked window, the SPECTRAL arc mask already removes the
worst contaminated pixels — so the F200 SPATIAL mask sensitivity is likely
smaller. This script measures it: same wlabel as run_window_sweep.WINDOWS
but with mask_weight varied across {0.0, 0.5, 1.0}.

Run:
    python scripts/run_maskweight_window.py --wlabel wR3800_5400_arcmask --n-bootstrap 100
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import numpy as np

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402
from scripts.bootstrap_ppxf_parallel import (  # noqa: E402
    run_bootstrap_single_degree_parallel,
)
from scripts.final_sigma_e import load_setup, extract_aperture_spectrum  # noqa: E402
from scripts.run_window_sweep import (  # noqa: E402
    WINDOWS as WIN_LIST, _apply_arc_mask, _apply_bad_pixels_mask,
    SPS_LAM_RANGE_TEMP, sps_safe_obs_min, Z_SYSTEMIC,
)

SPS_LIBS = ["fsps", "emiles", "xsl"]
DEGREES = np.arange(15, 30)
WINDOW = 75
BOOT_SEED_BASE = 42 + 110000

MASK_WEIGHTS = [0.0, 0.5, 1.0]   # 0=hard mask (headline), 0.5=soft, 1.0=no mask


def _lookup_window(wlabel):
    for w, lr, _arc in WIN_LIST:
        if w == wlabel:
            return lr, wlabel.endswith("_arcmask")
    raise SystemExit(f"--wlabel {wlabel!r} not in run_window_sweep.WINDOWS")


def run_one(mask_weight, sps, state, lam_obs_range, arc_mask,
            cache_dir, n_bootstrap, n_jobs, seed_offset, prefix=""):
    R_E = state["R_E"]
    mw_tag = f"w{int(round(mask_weight*100)):02d}"  # w00, w50, w100
    cache_fn = os.path.join(cache_dir, f"{mw_tag}_{sps}_N{n_bootstrap}.npz")
    if os.path.exists(cache_fn):
        d = np.load(cache_fn, allow_pickle=True)
        if "sig_boot" in d.files and d["sig_boot"].shape == (len(DEGREES), n_bootstrap):
            print(f"{prefix}[skip] cached {cache_fn}", flush=True)
            return cache_fn

    flux, noise, n_kept, sn = extract_aperture_spectrum(
        state, r_max=R_E, mask_weight=mask_weight,
    )
    n_floor = float(np.nanmedian(noise[noise > 0])) * 0.1
    noise = np.where(noise > 0, noise, n_floor)

    lr_eff = (max(lam_obs_range[0], sps_safe_obs_min(sps)), lam_obs_range[1])
    print(f"{prefix}[mask_weight={mask_weight:.2f} | {sps:6s}] N_kept={n_kept} S/N={sn:.1f} "
          f"obs=[{lr_eff[0]:.0f},{lr_eff[1]:.0f}] arc={'Y' if arc_mask else 'N'}", flush=True)

    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, state['hdr'], sps_name=sps, z=Z_SYSTEMIC,
        lam_obs_range=lr_eff, lam_fit_range=lr_eff,
        lam_range_temp=SPS_LAM_RANGE_TEMP[sps],
        verbose=False, frame_galaxy='auto',
        mask_balmer=not arc_mask,
    )
    if arc_mask:
        inputs['goodpixels'] = _apply_arc_mask(inputs['goodpixels'], inputs['lam_gal_rest'])
        inputs['goodpixels'] = _apply_bad_pixels_mask(inputs['goodpixels'], inputs['lam_gal_rest'])

    n_pix = len(inputs["galaxy"])
    n_deg = len(DEGREES)
    bf = np.zeros((n_deg, n_pix))
    V_o = np.zeros(n_deg); sig_o = np.zeros(n_deg)
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
        bf[i] = pp.bestfit; V_o[i] = pp.sol[0]; sig_o[i] = pp.sol[1]
        rb = run_bootstrap_single_degree_parallel(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=n_bootstrap, window=WINDOW,
            seed=BOOT_SEED_BASE + seed_offset + i, n_jobs=n_jobs,
        )
        V_b[i] = rb["V_samples"]; sig_b[i] = rb["sigma_samples"]

    elapsed = time.time() - t0
    out = dict(
        V_orig=V_o, sig_orig=sig_o, V_boot=V_b, sig_boot=sig_b,
        best_fit=bf, galaxy=inputs["galaxy"], noise=inputs["noise"],
        lam_gal_rest=inputs["lam_gal_rest"], goodpixels=inputs["goodpixels"],
        degrees=np.asarray(DEGREES), r_max=R_E, sps=sps,
        mask_weight=float(mask_weight), n_kept=int(n_kept), sn_band=sn,
        n_bootstrap=n_bootstrap, lam_obs_range=np.asarray(lr_eff),
        arc_mask=arc_mask,
    )
    os.makedirs(cache_dir, exist_ok=True)
    np.savez(cache_fn, **out)
    print(f"{prefix}  → σ_orig={sig_o.min():.0f}-{sig_o.max():.0f} | "
          f"sig_boot p50={np.median(sig_b[np.isfinite(sig_b)]):.1f} | {elapsed:.1f}s | "
          f"saved {cache_fn}", flush=True)
    return cache_fn


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wlabel", required=True)
    parser.add_argument("--n-bootstrap", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()

    lam_obs_range, arc_mask = _lookup_window(args.wlabel)
    cache_dir = os.path.join("results", f"maskweight_sweep_{args.wlabel}")
    print(f"\nF200 spatial mask-weight sweep at wlabel={args.wlabel}")
    print(f"  obs range = {lam_obs_range}; spectral arc_mask = {arc_mask}")
    print(f"  mask weights = {MASK_WEIGHTS}")
    print(f"  N_BOOT = {args.n_bootstrap}; n_jobs = {args.n_jobs}")
    print(f"  cache dir = {cache_dir}")

    state = load_setup()
    print(f"\nLoaded state.\n")

    total = len(MASK_WEIGHTS) * len(SPS_LIBS)
    counter = 0
    t_total = time.time()
    for mw_idx, mw in enumerate(MASK_WEIGHTS):
        for sps_idx, sps in enumerate(SPS_LIBS):
            counter += 1
            seed_offset = 100 * mw_idx + 10 * sps_idx
            prefix = f"[{counter}/{total}] "
            try:
                run_one(mw, sps, state, lam_obs_range, arc_mask, cache_dir,
                        n_bootstrap=args.n_bootstrap, n_jobs=args.n_jobs,
                        seed_offset=seed_offset, prefix=prefix)
            except Exception as e:
                print(f"{prefix}FAILED on mw={mw}/{sps}: {type(e).__name__}: {e}", flush=True)
    print(f"\n=== DONE in {(time.time()-t_total)/60:.1f} min ===")


if __name__ == "__main__":
    main()
