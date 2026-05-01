#!/usr/bin/env python3
"""
Refresh nb07c §7 (discrete Gültekin) and §7b (flat-σ extrapolation) at
N=500 with frame-aware ppxf (SPS_NATIVE_FRAME = {fsps: vacuum, emiles:
air, xsl: air}).

The original nb07c caches in results/annular_bootstrap_07c/ are from
2026-04-24 — pre-frame-fix (commit 9983720, 2026-04-28). The σ shift
from the frame fix is small (≤5 km/s per SPS, ≤2 km/s pooled) but
should be propagated for the published cross-checks H1 / H2 in
TESTS_AND_DIAGNOSTICS.

Equal-N inner binning (5 bins inside R_safe = 3R_e/4) + 1 outer flagged
bin (R_safe < R < R_e) — same as nb07c.

Outputs:
    results/refresh_07c_gultekin_postframe.npz
"""
from __future__ import annotations
import sys
from pathlib import Path
from time import perf_counter as clock

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402
from scripts.bootstrap_ppxf_parallel import (  # noqa: E402
    run_bootstrap_single_degree_parallel as run_bootstrap_single_degree,
)
from scripts.final_sigma_e import (  # noqa: E402
    SPS_LIBS, DEGREES, BOOT_SEED, WINDOW, Z_SYSTEMIC,
    KPC_PER_ARCSEC, load_setup,
)


N_BOOT = 500
N_INNER = 5
R_SAFE_FRAC = 0.75
SIGMA_MAX_CAP = 400.0  # arc-contamination flag
N_MC = 500             # MC draws for the Gültekin integral
CACHE_DIR = REPO / "results" / "refresh_07c_gultekin"
SUMMARY = REPO / "results" / "refresh_07c_gultekin_postframe.npz"


def annular_aperture_spectrum(state, r_lo, r_hi):
    """I-weighted spectrum of spaxels in [r_lo, r_hi] with F200 mask hard-applied."""
    sel = (state["r_spax"] >= r_lo) & (state["r_spax"] < r_hi) & ~state["arc_spax_mask"]
    n = int(sel.sum())
    I = np.clip(state["ifu_band"][sel], 0, None)
    w = I / max(I.sum(), 1e-30)
    flux = np.sum(state["cube"][:, sel] * w[None, :], axis=1)
    sn = float(
        np.median(flux[state["band_mask"]])
        / np.median(state["noise_sky"][state["band_mask"]])
    )
    return flux, state["noise_sky"].copy(), n, sn


def run_annulus_sps(state, j, r_lo, r_hi, sps, force=False):
    cache = CACHE_DIR / f"ann{j}_{sps}_N{N_BOOT}.npz"
    if cache.exists() and not force:
        d = dict(np.load(cache, allow_pickle=True))
        for k in list(d.keys()):
            if d[k].shape == ():
                d[k] = d[k].item()
        return d
    flux, noise, n_kept, sn = annular_aperture_spectrum(state, r_lo, r_hi)
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, state["hdr"], sps_name=sps, z=Z_SYSTEMIC,
        verbose=False, frame_galaxy="auto",
    )
    n_pix = len(inputs["galaxy"])
    n_deg = len(DEGREES)
    V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
    V_b = np.full((n_deg, N_BOOT), np.nan)
    sig_b = np.full((n_deg, N_BOOT), np.nan)
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
        rb = run_bootstrap_single_degree(
            inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
            n_bootstrap=N_BOOT, window=WINDOW,
            seed=BOOT_SEED + 100 * j + i, n_jobs=8,
        )
        V_b[i] = rb["V_samples"]
        sig_b[i] = rb["sigma_samples"]
    elapsed = clock() - t0
    out = dict(
        ann_idx=j, r_lo=float(r_lo), r_hi=float(r_hi),
        sps_name=sps, n_spax=int(n_kept), sn_band=float(sn),
        degrees=np.asarray(DEGREES),
        V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
        V_boot=V_b, sig_boot=sig_b,
        frame_galaxy=inputs["frame_galaxy"],
        n_bootstrap=int(N_BOOT), elapsed_s=float(elapsed),
    )
    np.savez(cache, **out)
    print(f"    [ann{j}/{sps:6s}] σ={sig_orig.min():.0f}-{sig_orig.max():.0f}  "
          f"V={V_orig.mean():+.1f}  frame={inputs['frame_galaxy']:7s}  "
          f"{elapsed:.0f}s  ({n_kept} spax, S/N={sn:.1f})")
    return out


def main():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print("═" * 70)
    print("REFRESH nb07c §7 / §7b at N=500 with frame-aware ppxf")
    print("═" * 70)
    state = load_setup()

    # Equal-N inner bins inside R_safe + 1 outer bin
    R_SAFE = R_SAFE_FRAC * state["R_E"]
    r_inside = np.sort(state["r_spax"][state["r_spax"] < R_SAFE].ravel())
    inner_edges = np.quantile(r_inside, np.linspace(0, 1, N_INNER + 1))
    r_edges = np.concatenate([inner_edges, [state["R_E"]]])
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    n_ann = len(r_edges) - 1
    print(f"\n  R_safe = {R_SAFE:.3f}\"  ({N_INNER} equal-N inner + 1 outer = {n_ann})")
    print(f"  r_edges = [{', '.join(f'{e:.3f}' for e in r_edges)}]")

    # ── Per-annulus per-SPS bootstrap ─────────────────────────────────────
    print(f"\n§B  Per-annulus ppxf+bootstrap, frame-aware, N={N_BOOT}")
    print("─" * 70)
    annuli = []
    for j in range(n_ann):
        print(f"  ann{j}:  r ∈ [{r_edges[j]:.3f}, {r_edges[j+1]:.3f}]\"")
        per_sps = {}
        for sps in SPS_LIBS:
            per_sps[sps] = run_annulus_sps(state, j, r_edges[j], r_edges[j+1], sps)
        annuli.append(dict(
            j=j, r_lo=float(r_edges[j]), r_hi=float(r_edges[j+1]),
            r_mid=float(r_mid[j]), per_sps=per_sps,
        ))

    # ── Per-SPS V_sys subtraction (median of innermost annulus) ───────────
    print(f"\n§C  Per-SPS V_sys offsets (median of ann0)")
    V_sys_per_sps = {}
    for sps in SPS_LIBS:
        Vs = annuli[0]["per_sps"][sps]["V_boot"].ravel()
        Vs = Vs[np.isfinite(Vs)]
        V_sys_per_sps[sps] = float(np.median(Vs))
        print(f"    {sps:6s}: V_sys = {V_sys_per_sps[sps]:+.2f} km/s  "
              f"(frame={annuli[0]['per_sps'][sps]['frame_galaxy']})")

    # ── F_j per annulus from IFU band weights ─────────────────────────────
    F_j = []
    for j in range(n_ann):
        sel = ((state["r_spax"] >= r_edges[j])
               & (state["r_spax"] < r_edges[j+1])
               & ~state["arc_spax_mask"])
        F_j.append(float(np.clip(state["ifu_band"][sel], 0, None).sum()))
    F_j = np.asarray(F_j)
    print(f"\n  F_j (IFU-band integrated flux per annulus): "
          f"[{', '.join(f'{f:.2f}' for f in F_j)}]")

    # ── §7 discrete Gültekin sum (arc-filtered) ───────────────────────────
    # σ_e²(<R) = Σ_j F_j (V_j² + σ_j²) / Σ_j F_j
    # arc filter: drop annuli where σ > SIGMA_MAX_CAP at any SPS (median σ)
    print(f"\n§7  Discrete Gültekin annular sum (arc-filter at σ > {SIGMA_MAX_CAP} km/s)")
    rng = np.random.default_rng(BOOT_SEED + 5000)
    arc_flag = np.zeros(n_ann, bool)
    for j in range(n_ann):
        median_sigs = []
        for sps in SPS_LIBS:
            d = annuli[j]["per_sps"][sps]
            median_sigs.append(float(np.median(d["sig_boot"][np.isfinite(d["sig_boot"])])))
        if max(median_sigs) > SIGMA_MAX_CAP:
            arc_flag[j] = True
            print(f"    ann{j} flagged (max σ = {max(median_sigs):.0f} km/s)")
    keep_mask = ~arc_flag
    F_keep = F_j[keep_mask]
    F_total = F_keep.sum()

    # MC pool: per-iteration draw σ_j and V_j from each annulus's bootstrap
    # samples, subtract per-SPS V_sys, compute Gültekin sum
    sigma_e_samples = []
    for _ in range(N_MC):
        # Pool: pick one bootstrap iter per (annulus, SPS) and average
        num = 0.0
        for j in range(n_ann):
            if not keep_mask[j]:
                continue
            term = 0.0
            n_sps_used = 0
            for sps in SPS_LIBS:
                d = annuli[j]["per_sps"][sps]
                vs = d["V_boot"].ravel()
                ss = d["sig_boot"].ravel()
                ok = np.isfinite(vs) & np.isfinite(ss)
                if ok.sum() == 0:
                    continue
                idx = rng.integers(ok.sum())
                v = vs[ok][idx] - V_sys_per_sps[sps]
                s = ss[ok][idx]
                term += v**2 + s**2
                n_sps_used += 1
            if n_sps_used > 0:
                num += F_j[j] * (term / n_sps_used)
        sigma_e_samples.append(np.sqrt(num / F_total))
    sigma_e_samples = np.asarray(sigma_e_samples)
    s7_p16, s7_p50, s7_p84 = np.percentile(sigma_e_samples, [16, 50, 84])
    print(f"\n  §7  σ_e (filtered) = {s7_p50:.2f} -{s7_p50-s7_p16:.2f} +{s7_p84-s7_p50:.2f} km/s "
          f"(N_MC={N_MC})")

    # ── §7b flat-σ extrapolation: keep all bins, override outer (V, σ) ─────
    last_clean_j = max(j for j in range(n_ann) if keep_mask[j])
    print(f"\n§7b  Flat-σ extrapolation — outer bin(s) inherit (V, σ) from ann{last_clean_j}")
    sigma_e_b_samples = []
    for _ in range(N_MC):
        num = 0.0
        for j in range(n_ann):
            term = 0.0
            n_sps_used = 0
            j_src = j if keep_mask[j] else last_clean_j
            for sps in SPS_LIBS:
                d = annuli[j_src]["per_sps"][sps]
                vs = d["V_boot"].ravel()
                ss = d["sig_boot"].ravel()
                ok = np.isfinite(vs) & np.isfinite(ss)
                if ok.sum() == 0:
                    continue
                idx = rng.integers(ok.sum())
                v = vs[ok][idx] - V_sys_per_sps[sps]
                s = ss[ok][idx]
                term += v**2 + s**2
                n_sps_used += 1
            if n_sps_used > 0:
                num += F_j[j] * (term / n_sps_used)
        sigma_e_b_samples.append(np.sqrt(num / F_j.sum()))
    sigma_e_b_samples = np.asarray(sigma_e_b_samples)
    s7b_p16, s7b_p50, s7b_p84 = np.percentile(sigma_e_b_samples, [16, 50, 84])
    print(f"  §7b σ_e            = {s7b_p50:.2f} -{s7b_p50-s7b_p16:.2f} +{s7b_p84-s7b_p50:.2f} km/s "
          f"(N_MC={N_MC})")

    # ── Pre-fix anchors from TESTS_AND_DIAGNOSTICS catalogue ──────────────
    s7_pre = (254.99, 24.2, 28.4)
    s7b_pre = (271.0, 33.0, 35.0)
    print(f"\n§D  Comparison with pre-frame-fix nb07c (Apr 24)")
    print("─" * 70)
    print(f"  §7  σ_e:  pre = {s7_pre[0]:7.2f} -{s7_pre[1]:.2f}/+{s7_pre[2]:.2f}   "
          f"post = {s7_p50:7.2f} -{s7_p50-s7_p16:.2f}/+{s7_p84-s7_p50:.2f}   "
          f"Δ = {s7_p50 - s7_pre[0]:+.2f}")
    print(f"  §7b σ_e:  pre = {s7b_pre[0]:7.2f} -{s7b_pre[1]:.2f}/+{s7b_pre[2]:.2f}   "
          f"post = {s7b_p50:7.2f} -{s7b_p50-s7b_p16:.2f}/+{s7b_p84-s7b_p50:.2f}   "
          f"Δ = {s7b_p50 - s7b_pre[0]:+.2f}")

    np.savez(
        SUMMARY,
        z_systemic=Z_SYSTEMIC, R_E=state["R_E"], R_safe=R_SAFE,
        r_edges=r_edges, r_mid=r_mid, F_j=F_j, arc_flag=arc_flag,
        last_clean_j=last_clean_j, n_bootstrap=N_BOOT, n_mc=N_MC,
        sps_libs=np.array(SPS_LIBS),
        V_sys_per_sps=np.array([V_sys_per_sps[s] for s in SPS_LIBS]),
        s7_p16=s7_p16, s7_p50=s7_p50, s7_p84=s7_p84,
        s7b_p16=s7b_p16, s7b_p50=s7b_p50, s7b_p84=s7b_p84,
        s7_samples=sigma_e_samples, s7b_samples=sigma_e_b_samples,
        s7_pre=s7_pre, s7b_pre=s7b_pre,
        delta_s7=float(s7_p50 - s7_pre[0]),
        delta_s7b=float(s7b_p50 - s7b_pre[0]),
    )
    print(f"\nSaved → {SUMMARY.relative_to(REPO)}")


if __name__ == "__main__":
    main()
