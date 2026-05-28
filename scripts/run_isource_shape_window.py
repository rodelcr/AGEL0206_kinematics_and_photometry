"""I(r)-shape sweep at a configurable fit window, with optional arc mask.

Variant of `scripts/run_isource_shape_sweep.py` that:
  - takes --wlabel (the window key from run_window_sweep.WINDOWS, used as
    suffix for the cache dir and to look up the obs-range)
  - applies the arc-emission masks if the wlabel ends with `_arcmask`
  - reuses the 10-shape construction from run_isource_shape_sweep
  - caches separately so it doesn't collide with the headline I-shape run

Run:
    python scripts/run_isource_shape_window.py --wlabel wR3800_5400_arcmask --n-bootstrap 100
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
from scripts.run_isource_sweep import load_setup  # noqa: E402
from scripts.run_isource_shape_sweep import build_all_shapes  # noqa: E402
from scripts.run_window_sweep import (  # noqa: E402
    WINDOWS as WIN_LIST, ARC_MASKS_REST, _apply_arc_mask,
    BAD_PIXELS_REST, _apply_bad_pixels_mask,
    SPS_LAM_RANGE_TEMP, sps_safe_obs_min, Z_SYSTEMIC,
)

SPS_LIBS = ["fsps", "emiles", "xsl"]
DEGREES = np.arange(15, 30)
WINDOW = 75
BOOT_SEED_BASE = 42 + 90000


def _lookup_window(wlabel):
    for w, lr, _arc in WIN_LIST:
        if w == wlabel:
            return lr, wlabel.endswith("_arcmask")
    raise SystemExit(f"--wlabel {wlabel!r} not in run_window_sweep.WINDOWS")


def run_one(shape_name, I_map, sps, state, lam_obs_range, arc_mask,
            cache_dir, n_bootstrap, n_jobs, seed_offset, prefix=""):
    cache_fn = os.path.join(cache_dir, f"{shape_name}_{sps}_N{n_bootstrap}.npz")
    if os.path.exists(cache_fn):
        d = np.load(cache_fn, allow_pickle=True)
        if "sig_boot" in d.files and d["sig_boot"].shape == (len(DEGREES), n_bootstrap):
            print(f"{prefix}[skip] cached {cache_fn}", flush=True)
            return cache_fn

    cube = state["cube"]
    r_spax = state["r_spax"]
    arc_spax_mask = state["arc_spax_mask"]
    R_E = state["R_E"]
    band_mask = state["band_mask"]
    noise_sky = state["noise_sky"]

    sel = (r_spax < R_E) & ~arc_spax_mask
    w = I_map[sel]
    if w.sum() <= 0:
        print(f"{prefix}FAIL: I-weight sum <= 0 for {shape_name}", flush=True)
        return None
    w_norm = w / w.sum()
    flux = np.sum(cube[:, sel] * w_norm[None, :], axis=1)
    n_kept = int(sel.sum())
    sn = float(np.median(flux[band_mask]) / np.median(noise_sky[band_mask]))

    # Per-SPS lower-bound clip (XSL floor at obs ~6158 Å)
    lr_eff = (max(lam_obs_range[0], sps_safe_obs_min(sps)), lam_obs_range[1])
    print(f"{prefix}[{shape_name:<18s} | {sps:6s}] N_kept={n_kept:3d} S/N={sn:5.1f}  "
          f"obs=[{lr_eff[0]:.0f},{lr_eff[1]:.0f}]  arc={'Y' if arc_mask else 'N'}", flush=True)

    # noise floor (cube has zero-noise pix at blue edge)
    n_floor = float(np.nanmedian(noise_sky[noise_sky > 0])) * 0.1
    noise_floor = np.where(noise_sky > 0, noise_sky, n_floor)

    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise_floor.copy(), state["hdr"], sps_name=sps, z=Z_SYSTEMIC,
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
        ishape=shape_name, n_kept=n_kept, sn_band=sn,
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
    parser.add_argument("--wlabel", required=True,
                        help="Window label from run_window_sweep.WINDOWS")
    parser.add_argument("--n-bootstrap", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--only-shape", default=None)
    parser.add_argument("--only-sps", default=None)
    args = parser.parse_args()

    lam_obs_range, arc_mask = _lookup_window(args.wlabel)
    cache_dir = os.path.join("results", f"ishape_sweep_{args.wlabel}")
    print(f"\nI-shape sweep at wlabel={args.wlabel}")
    print(f"  obs range = {lam_obs_range}; arc_mask = {arc_mask}")
    print(f"  N_BOOT = {args.n_bootstrap}; n_jobs = {args.n_jobs}")
    print(f"  cache dir = {cache_dir}")

    state = load_setup()
    I_shapes = build_all_shapes(state)
    print(f"\nBuilt {len(I_shapes)} I-shape maps.\n")

    shapes = list(I_shapes.keys())
    if args.only_shape is not None:
        shapes = [s for s in shapes if s == args.only_shape]
    sps_libs = SPS_LIBS if args.only_sps is None else [args.only_sps]

    total = len(shapes) * len(sps_libs)
    counter = 0
    t_total = time.time()
    print(f"=== Running {len(shapes)} I-shapes × {len(sps_libs)} SPS × N={args.n_bootstrap} ===\n")
    for s_idx, shape_name in enumerate(shapes):
        for sps_idx, sps in enumerate(sps_libs):
            counter += 1
            seed_offset = 100 * s_idx + 10 * sps_idx
            prefix = f"[{counter}/{total}] "
            try:
                run_one(shape_name, I_shapes[shape_name], sps, state,
                        lam_obs_range, arc_mask, cache_dir,
                        n_bootstrap=args.n_bootstrap, n_jobs=args.n_jobs,
                        seed_offset=seed_offset, prefix=prefix)
            except Exception as e:
                print(f"{prefix}FAILED on {shape_name}/{sps}: {type(e).__name__}: {e}", flush=True)
    print(f"\n=== DONE in {(time.time()-t_total)/60:.1f} min ===")


if __name__ == "__main__":
    main()
