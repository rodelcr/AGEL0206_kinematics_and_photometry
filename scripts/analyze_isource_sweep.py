"""Analyze the §6cum I(r) sweep results.

Loads all per-config caches in results/annular_bootstrap_07c_isource/, computes
combined-SPS posteriors (concatenation across SPS × degrees), prints a summary
table, and produces a 4-panel diagnostic figure.

Usage
-----
    conda activate ISMgas
    python scripts/analyze_isource_sweep.py
"""

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.run_isource_sweep import (  # noqa: E402
    NB06_VIZ_COMBOS, SEL_STRATEGIES, SPS_LIBS, CACHE_DIR,
    load_setup, build_I_weight_map,
)

FIG_DIR = "results/figures"
HEADLINE_KEY = ("IFU_band", "unmasked", "arcmask")


def collect_posterior(isrc, maskstrat, selstrat):
    """Concatenate sig_boot across SPS_LIBS × degrees → 1D combined posterior."""
    chunks = []
    for sps in SPS_LIBS:
        fn = os.path.join(CACHE_DIR, f"{isrc}_{maskstrat}_{selstrat}_{sps}.npz")
        if not os.path.exists(fn):
            return None, None
        d = np.load(fn, allow_pickle=True)
        sb = d["sig_boot"]
        chunks.append(sb.ravel())
    pool = np.concatenate(chunks)
    pool = pool[np.isfinite(pool)]
    p16, p50, p84 = np.percentile(pool, [16, 50, 84])
    n_kept = int(d["n_kept"]) if "n_kept" in d.files else None
    return pool, dict(p16=p16, p50=p50, p84=p84, n_kept=n_kept)


def per_sps_p50(isrc, maskstrat, selstrat):
    out = {}
    for sps in SPS_LIBS:
        fn = os.path.join(CACHE_DIR, f"{isrc}_{maskstrat}_{selstrat}_{sps}.npz")
        if not os.path.exists(fn):
            out[sps] = np.nan
            continue
        d = np.load(fn, allow_pickle=True)
        sb = d["sig_boot"].ravel()
        sb = sb[np.isfinite(sb)]
        out[sps] = float(np.median(sb)) if sb.size else np.nan
    return out


def main():
    # --------------------------------------------------------------
    # Build the summary table
    # --------------------------------------------------------------
    rows = []
    for isrc, maskstrat in NB06_VIZ_COMBOS:
        for selstrat in SEL_STRATEGIES:
            pool, stats = collect_posterior(isrc, maskstrat, selstrat)
            if pool is None:
                rows.append((isrc, maskstrat, selstrat, None))
                continue
            psps = per_sps_p50(isrc, maskstrat, selstrat)
            rows.append((isrc, maskstrat, selstrat, dict(stats=stats, psps=psps)))

    # Headline reference (from the sweep itself if present)
    headline = None
    for r in rows:
        if (r[0], r[1], r[2]) == HEADLINE_KEY and r[3] is not None:
            headline = r[3]["stats"]
            break

    if headline is None:
        # Fallback to the masked production cache if sweep doesn't yet have it
        try:
            chunks = []
            for sps in SPS_LIBS:
                fn = f"results/annular_bootstrap_07c/cumR_2p305_{sps}.npz"
                d = np.load(fn, allow_pickle=True)
                chunks.append(d["sig_boot"].ravel())
            pool = np.concatenate(chunks)
            pool = pool[np.isfinite(pool)]
            p16, p50, p84 = np.percentile(pool, [16, 50, 84])
            headline = dict(p16=p16, p50=p50, p84=p84, n_kept=None)
            print(f"[headline reference] from production cache: {p50:.2f} -{p50-p16:.2f}/+{p84-p50:.2f}")
        except Exception:
            headline = dict(p16=243.0, p50=267.32, p84=291.0, n_kept=147)
            print("[headline reference] hard-coded")

    head_p50 = headline["p50"]

    # --------------------------------------------------------------
    # Print table
    # --------------------------------------------------------------
    sep = "=" * 110
    print(f"\n{sep}")
    print(f"§6cum I(r) sweep × {{arcmask, nomask}} at R<R_e   N=500 boot × 3 SPS × 15 deg")
    print(f"Headline σ_e ({HEADLINE_KEY}): {head_p50:.2f} km/s")
    print(sep)
    print(f"{'I-source':<10} {'mask_strategy':<14} {'selstrat':<8} "
          f"{'p50':>7} {'-1σ':>6} {'+1σ':>6} "
          f"{'fsps':>7} {'emiles':>7} {'xsl':>7} "
          f"{'Δheadline':>10} {'N_kept':>7}")
    print("-" * 110)
    for isrc, mask, sel, info in rows:
        if info is None:
            print(f"{isrc:<10} {mask:<14} {sel:<8} {'(not yet cached)':>50}")
            continue
        s = info["stats"]
        psps = info["psps"]
        nk = s["n_kept"] if s["n_kept"] is not None else 0
        print(
            f"{isrc:<10} {mask:<14} {sel:<8} "
            f"{s['p50']:>7.1f} {s['p50']-s['p16']:>6.1f} {s['p84']-s['p50']:>6.1f} "
            f"{psps['fsps']:>7.1f} {psps['emiles']:>7.1f} {psps['xsl']:>7.1f} "
            f"{s['p50']-head_p50:>+10.1f} {nk:>7d}"
        )
    print(sep)

    # Spreads
    have = [r for r in rows if r[3] is not None]
    arc_rows = [r for r in have if r[2] == "arcmask"]
    no_rows = [r for r in have if r[2] == "nomask"]
    arc_p50 = np.array([r[3]["stats"]["p50"] for r in arc_rows])
    no_p50 = np.array([r[3]["stats"]["p50"] for r in no_rows])
    print(f"\nI-source spread within arcmask: min={arc_p50.min():.1f} max={arc_p50.max():.1f} "
          f"range={arc_p50.max()-arc_p50.min():.1f} km/s")
    print(f"I-source spread within nomask:  min={no_p50.min():.1f} max={no_p50.max():.1f} "
          f"range={no_p50.max()-no_p50.min():.1f} km/s")

    arcmask_avg = np.mean(arc_p50)
    nomask_avg = np.mean(no_p50)
    print(f"\nMean σ_e (arcmask): {arcmask_avg:.1f} km/s")
    print(f"Mean σ_e (nomask):  {nomask_avg:.1f} km/s")
    print(f"Mask effect (mean):  {nomask_avg-arcmask_avg:+.1f} km/s")

    # --------------------------------------------------------------
    # 4-panel figure
    # --------------------------------------------------------------
    state = load_setup()

    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.0, 1.0], hspace=0.30, wspace=0.30)

    # Panel A: bar chart of σ_e with errorbars (mask state side-by-side)
    axA = fig.add_subplot(gs[0, :])
    labels = []
    p50s, p16s, p84s, colors = [], [], [], []
    color_map = {"arcmask": "C0", "nomask": "C3"}
    for isrc, mask in NB06_VIZ_COMBOS:
        for sel in SEL_STRATEGIES:
            info = next((r for r in rows if r[0] == isrc and r[1] == mask and r[2] == sel), None)
            if info is None or info[3] is None:
                continue
            s = info[3]["stats"]
            labels.append(f"{isrc}\n{mask}\n{sel}")
            p50s.append(s["p50"])
            p16s.append(s["p50"] - s["p16"])
            p84s.append(s["p84"] - s["p50"])
            colors.append(color_map[sel])
    x = np.arange(len(p50s))
    axA.bar(x, p50s, color=colors, alpha=0.75, edgecolor="black")
    axA.errorbar(x, p50s, yerr=[p16s, p84s], fmt="none", color="black", capsize=3, lw=1)
    axA.axhline(head_p50, color="black", ls="--", lw=1.0, alpha=0.6,
                label=f"Headline σ_e = {head_p50:.1f}")
    axA.set_xticks(x)
    axA.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axA.set_ylabel("σ_e (<R_e) [km/s]")
    axA.set_title("§6cum I(r) sweep × {arcmask, nomask}  —  combined SPS, N=500 bootstrap")
    axA.legend(loc="upper left")
    axA.grid(alpha=0.3, axis="y")

    # Panel B: posterior overlays for representative cases
    axB = fig.add_subplot(gs[1, 0:2])
    show_cases = [
        ("IFU_band", "unmasked", "arcmask",  "C0", "IFU_band·unmasked·arcmask (HEADLINE)"),
        ("IFU_band", "unmasked", "nomask",   "C3", "IFU_band·unmasked·nomask"),
        ("IFU_band", "15pct_psf", "nomask",  "C2", "IFU_band·15pct_psf·nomask"),
        ("F140W",    "unmasked", "arcmask",  "C4", "F140W·unmasked·arcmask"),
        ("F140W",    "arc_only_ifu", "nomask", "C5", "F140W·arc_only_ifu·nomask"),
        ("F200LP",   "arc_only_ifu", "nomask", "C6", "F200LP·arc_only_ifu·nomask"),
    ]
    for isrc, mask, sel, c, lab in show_cases:
        pool, _ = collect_posterior(isrc, mask, sel)
        if pool is None:
            continue
        axB.hist(pool, bins=80, density=True, histtype="step", color=c, lw=1.3, alpha=0.85, label=lab)
    axB.axvline(head_p50, color="black", ls="--", lw=1, alpha=0.6)
    axB.set_xlabel("σ_e (<R_e) [km/s]")
    axB.set_ylabel("density")
    axB.set_title("Combined-SPS σ_e posteriors (representative cases)")
    axB.legend(loc="upper right", fontsize=8)
    axB.grid(alpha=0.3)

    # Panel C: I-maps small multiples
    cx, cy = state["cx_ifu"], state["cy_ifu"]
    R_E = state["R_E"]
    pix_scale = float(np.abs(np.diag(state["wcs_ifu"].pixel_scale_matrix))[0]) * 3600 if False else 0.30
    delta = 18
    for k, (isrc, mask) in enumerate(NB06_VIZ_COMBOS):
        if k >= 4:
            break  # show first 4 in first row of small multiples
        ax = fig.add_subplot(gs[1, 2 + k] if k < 2 else gs[2, k-2])
        try:
            I = build_I_weight_map(isrc, mask, state)
        except Exception as e:
            ax.set_title(f"{isrc}/{mask}\nERROR: {type(e).__name__}", fontsize=8)
            ax.axis("off")
            continue
        ax.imshow(I, origin="lower", cmap="viridis",
                  vmin=np.percentile(I[I > 0], 5) if (I > 0).any() else 0,
                  vmax=np.percentile(I[I > 0], 99) if (I > 0).any() else 1)
        ax.add_patch(Circle((cx, cy), R_E / pix_scale, fill=False, edgecolor="white", lw=1.4))
        ax.plot(cx, cy, "w+", ms=10, mew=1.5)
        ax.set_xlim(int(cx-delta), int(cx+delta))
        ax.set_ylim(int(cy-delta), int(cy+delta))
        ax.set_title(f"{isrc}\n{mask}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])

    # Bottom row: aperture spectra of three representative cases
    axS = fig.add_subplot(gs[2, 2:])
    show_specs = [
        ("IFU_band", "unmasked", "arcmask", "C0"),
        ("F140W",    "unmasked", "nomask",  "C4"),
        ("F200LP",   "arc_only_ifu", "nomask", "C6"),
    ]
    for isrc, mask, sel, c in show_specs:
        fn = os.path.join(CACHE_DIR, f"{isrc}_{mask}_{sel}_fsps.npz")
        if not os.path.exists(fn):
            continue
        d = np.load(fn, allow_pickle=True)
        gal = d["galaxy"]
        lam = d["lam_gal_rest"]
        axS.plot(lam, gal / np.median(gal), lw=0.8, color=c, alpha=0.9,
                 label=f"{isrc}/{mask}/{sel}")
    axS.set_xlabel("rest-frame λ [Å]")
    axS.set_ylabel("flux (normalized)")
    axS.set_title("Aperture spectra (FSPS bin) — representative")
    axS.legend(fontsize=8)
    axS.grid(alpha=0.3)

    fig.suptitle(
        f"§6cum I(r) systematic sweep   |   headline σ_e = {head_p50:.1f} km/s "
        f"(IFU_band, unmasked, arcmask)",
        fontsize=14,
    )
    out_png = os.path.join(FIG_DIR, "nb07c_isource_sweep.png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"\nSaved → {out_png}")

    # Also save the table as npz for downstream loading
    out_npz = "results/sigma_e_isource_sweep.npz"
    arr_isrc = []; arr_mask = []; arr_sel = []
    arr_p16 = []; arr_p50 = []; arr_p84 = []
    arr_n_kept = []; arr_dh = []
    for isrc, mask, sel, info in rows:
        if info is None:
            continue
        s = info["stats"]
        arr_isrc.append(isrc); arr_mask.append(mask); arr_sel.append(sel)
        arr_p16.append(s["p16"]); arr_p50.append(s["p50"]); arr_p84.append(s["p84"])
        arr_n_kept.append(s["n_kept"] if s["n_kept"] is not None else -1)
        arr_dh.append(s["p50"] - head_p50)
    np.savez(
        out_npz,
        isource=np.array(arr_isrc), maskstrat=np.array(arr_mask), selstrat=np.array(arr_sel),
        p16=np.array(arr_p16), p50=np.array(arr_p50), p84=np.array(arr_p84),
        n_kept=np.array(arr_n_kept), delta_headline=np.array(arr_dh),
        headline_p50=head_p50, headline_p16=headline["p16"], headline_p84=headline["p84"],
    )
    print(f"Saved → {out_npz}")


if __name__ == "__main__":
    main()
