"""Analyze the §6cum I(r)-shape sweep results.

Loads results/annular_bootstrap_07c_ishape/, builds combined-SPS posteriors
and a summary table, and produces a 4-panel figure focused on the I-shape
systematic (mask is held fixed at F200-raw across all 10 shapes).
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.run_isource_shape_sweep import (  # noqa: E402
    build_all_shapes, SPS_LIBS, CACHE_DIR,
)
from scripts.run_isource_sweep import load_setup  # noqa: E402

FIG_DIR = "results/figures"
HEADLINE_SHAPE = "IFU_band"

SHAPE_ORDER = [
    "IFU_band", "IFU_wl",
    "F140W_raw", "F200LP_raw",
    "F140W_arcmasked", "F200LP_arcmasked",
    "F140W_1Dcog", "F200LP_1Dcog",
    "F140W_Sersic2D", "F200LP_Sersic2D",
]

GROUP_OF = {
    "IFU_band": "IFU",
    "IFU_wl": "IFU",
    "F140W_raw": "raw HST",
    "F200LP_raw": "raw HST",
    "F140W_arcmasked": "arc-masked HST",
    "F200LP_arcmasked": "arc-masked HST",
    "F140W_1Dcog": "1D CoG smooth",
    "F200LP_1Dcog": "1D CoG smooth",
    "F140W_Sersic2D": "2D Sersic smooth",
    "F200LP_Sersic2D": "2D Sersic smooth",
}

GROUP_COLORS = {
    "IFU": "#1f77b4",
    "raw HST": "#2ca02c",
    "arc-masked HST": "#9467bd",
    "1D CoG smooth": "#ff7f0e",
    "2D Sersic smooth": "#d62728",
}


def collect(shape, sps_libs=SPS_LIBS):
    chunks = []
    n_kept = None
    for sps in sps_libs:
        fn = os.path.join(CACHE_DIR, f"{shape}_{sps}.npz")
        if not os.path.exists(fn):
            return None
        d = np.load(fn, allow_pickle=True)
        chunks.append(d["sig_boot"].ravel())
        n_kept = int(d["n_kept"]) if "n_kept" in d.files else None
    pool = np.concatenate(chunks)
    pool = pool[np.isfinite(pool)]
    p16, p50, p84 = np.percentile(pool, [16, 50, 84])
    psps = {}
    for sps in sps_libs:
        d = np.load(os.path.join(CACHE_DIR, f"{shape}_{sps}.npz"), allow_pickle=True)
        sb = d["sig_boot"].ravel()
        sb = sb[np.isfinite(sb)]
        psps[sps] = float(np.median(sb))
    return dict(pool=pool, p16=p16, p50=p50, p84=p84, psps=psps, n_kept=n_kept)


def main():
    rows = []
    for shape in SHAPE_ORDER:
        info = collect(shape)
        rows.append((shape, info))

    head_p50 = next(r[1]["p50"] for r in rows if r[0] == HEADLINE_SHAPE)

    # Read actual N_bootstrap from a representative cache for honest reporting
    _probe = np.load(os.path.join(CACHE_DIR, f"{HEADLINE_SHAPE}_{SPS_LIBS[0]}.npz"),
                      allow_pickle=True)
    n_boot_actual = (int(_probe["n_bootstrap"]) if "n_bootstrap" in _probe.files
                     else int(_probe["sig_boot"].shape[1]))

    sep = "=" * 110
    print(f"\n{sep}")
    print(f"§6cum I(r)-SHAPE sweep at R<R_e   (mask FIXED at F200-raw)   "
          f"N={n_boot_actual} boot × 3 SPS × 15 deg")
    print(f"Headline σ_e ({HEADLINE_SHAPE}): {head_p50:.2f} km/s")
    print(sep)
    print(f"{'I-shape':<20} {'group':<19} "
          f"{'p50':>7} {'-1σ':>5} {'+1σ':>5} "
          f"{'fsps':>7} {'emiles':>7} {'xsl':>7} "
          f"{'Δheadline':>10} {'N_kept':>7}")
    print("-" * 110)
    for shape, info in rows:
        if info is None:
            print(f"{shape:<20} (not yet cached)")
            continue
        psps = info["psps"]
        nk = info["n_kept"] or 0
        print(
            f"{shape:<20} {GROUP_OF[shape]:<19} "
            f"{info['p50']:>7.1f} {info['p50']-info['p16']:>5.1f} {info['p84']-info['p50']:>5.1f} "
            f"{psps['fsps']:>7.1f} {psps['emiles']:>7.1f} {psps['xsl']:>7.1f} "
            f"{info['p50']-head_p50:>+10.1f} {nk:>7d}"
        )
    print(sep)

    have = [(s, i) for s, i in rows if i is not None]
    p50s = np.array([i["p50"] for _, i in have])
    print(f"\nFull spread (10 shapes): min={p50s.min():.1f}  max={p50s.max():.1f}  "
          f"range={p50s.max()-p50s.min():.1f} km/s")

    # Excluding the F200LP_Sersic2D outlier (poor fit, n=0.30 boundary)
    p50s_ok = np.array([i["p50"] for s, i in have if s != "F200LP_Sersic2D"])
    print(f"Excluding F200LP_Sersic2D (poor fit): min={p50s_ok.min():.1f}  max={p50s_ok.max():.1f}  "
          f"range={p50s_ok.max()-p50s_ok.min():.1f} km/s")

    # Sub-group spreads
    for grp in GROUP_COLORS:
        vals = [i["p50"] for s, i in have if GROUP_OF[s] == grp]
        if len(vals) >= 2:
            print(f"  {grp:<18s} spread = {max(vals)-min(vals):.1f} km/s "
                  f"({len(vals)} shapes)")

    # ─── 4-panel figure ───
    state = load_setup()
    print("\nBuilding I-shape maps for figure (re-runs Sersic fits)...")
    I_shapes = build_all_shapes(state)

    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 5, height_ratios=[1.0, 1.0, 1.0], hspace=0.35, wspace=0.30)

    # Panel A: bar chart with grouped colors
    axA = fig.add_subplot(gs[0, :])
    p50_arr, p16_arr, p84_arr, lbls, cols = [], [], [], [], []
    for shape, info in rows:
        if info is None:
            continue
        p50_arr.append(info["p50"])
        p16_arr.append(info["p50"] - info["p16"])
        p84_arr.append(info["p84"] - info["p50"])
        lbls.append(shape)
        cols.append(GROUP_COLORS[GROUP_OF[shape]])
    x = np.arange(len(p50_arr))
    axA.bar(x, p50_arr, color=cols, alpha=0.8, edgecolor="black")
    axA.errorbar(x, p50_arr, yerr=[p16_arr, p84_arr], fmt="none", color="black", capsize=4, lw=1.2)
    axA.axhline(head_p50, color="black", ls="--", lw=1.2, alpha=0.6,
                label=f"Headline {HEADLINE_SHAPE} = {head_p50:.1f}")
    axA.set_xticks(x)
    axA.set_xticklabels(lbls, rotation=30, ha="right", fontsize=10)
    axA.set_ylabel(r"$\sigma_e (<R_e)$ [km/s]")
    axA.set_title(f"§6cum I(r)-shape sweep — mask FIXED at F200-raw (N={n_boot_actual} × 3 SPS × 15 deg)\n"
                  f"All shapes use sel = (r_spax < R_e) ∧ ¬arc_spax_mask")
    # Legend handles for groups
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=c, label=g) for g, c in GROUP_COLORS.items()]
    legend_handles.insert(0, plt.Line2D([0], [0], color="black", ls="--", label=f"Headline σ_e"))
    axA.legend(handles=legend_handles, loc="upper left", fontsize=9, ncol=2)
    axA.grid(alpha=0.3, axis="y")

    # Panel B: posterior overlay (one per group)
    axB = fig.add_subplot(gs[1, 0:3])
    rep = ["IFU_band", "F140W_raw", "F140W_arcmasked", "F140W_1Dcog", "F140W_Sersic2D",
           "F200LP_raw", "F200LP_Sersic2D"]
    color_per_shape = {
        "IFU_band": "C0", "F140W_raw": "C2", "F140W_arcmasked": "C5",
        "F140W_1Dcog": "C1", "F140W_Sersic2D": "C3",
        "F200LP_raw": "C4", "F200LP_Sersic2D": "C6",
    }
    for shape in rep:
        info = next((i for s, i in rows if s == shape), None)
        if info is None:
            continue
        axB.hist(info["pool"], bins=80, density=True, histtype="step",
                 color=color_per_shape[shape], lw=1.4, alpha=0.85,
                 label=f"{shape}: {info['p50']:.0f}")
    axB.axvline(head_p50, color="black", ls="--", lw=1, alpha=0.6)
    axB.set_xlabel(r"$\sigma_e (<R_e)$ [km/s]")
    axB.set_ylabel("density")
    axB.set_title("Combined-SPS σ_e posteriors")
    axB.legend(loc="upper right", fontsize=9)
    axB.grid(alpha=0.3)

    # Panel C: I-maps small multiples (5 representative shapes)
    cx, cy = state["cx_ifu"], state["cy_ifu"]
    R_E = state["R_E"]
    pix = 0.30
    delta = 18

    # Show 5 of the 10 shapes, picked to span the groups
    show_shapes = ["IFU_band", "F140W_raw", "F140W_arcmasked", "F140W_1Dcog", "F140W_Sersic2D"]
    for k, sh in enumerate(show_shapes):
        ax = fig.add_subplot(gs[1, 3 + k] if k < 2 else gs[2, k-2])
        I = I_shapes[sh]
        pos = I[I > 0]
        ax.imshow(I, origin="lower", cmap="viridis",
                  vmin=np.percentile(pos, 5) if pos.size else 0,
                  vmax=np.percentile(pos, 99) if pos.size else 1)
        ax.add_patch(Circle((cx, cy), R_E / pix, fill=False, edgecolor="white", lw=1.5))
        ax.plot(cx, cy, "w+", ms=10, mew=1.5)
        ax.set_xlim(int(cx - delta), int(cx + delta))
        ax.set_ylim(int(cy - delta), int(cy + delta))
        info = next((i for s, i in rows if s == sh), None)
        sig = info["p50"] if info else float("nan")
        ax.set_title(f"{sh}\nσ_e={sig:.0f}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])

    # Bottom-right: F200LP shapes
    show_f200 = ["F200LP_raw", "F200LP_arcmasked", "F200LP_1Dcog", "F200LP_Sersic2D"]
    for k, sh in enumerate(show_f200):
        ax = fig.add_subplot(gs[2, 1 + k])
        I = I_shapes[sh]
        pos = I[I > 0]
        ax.imshow(I, origin="lower", cmap="viridis",
                  vmin=np.percentile(pos, 5) if pos.size else 0,
                  vmax=np.percentile(pos, 99) if pos.size else 1)
        ax.add_patch(Circle((cx, cy), R_E / pix, fill=False, edgecolor="white", lw=1.5))
        ax.plot(cx, cy, "w+", ms=10, mew=1.5)
        ax.set_xlim(int(cx - delta), int(cx + delta))
        ax.set_ylim(int(cy - delta), int(cy + delta))
        info = next((i for s, i in rows if s == sh), None)
        sig = info["p50"] if info else float("nan")
        ax.set_title(f"{sh}\nσ_e={sig:.0f}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(
        f"§6cum I(r)-shape systematic — mask fixed at F200-raw   "
        f"|   headline σ_e = {head_p50:.1f} km/s\n"
        f"Full 10-shape spread = {p50s.max()-p50s.min():.1f} km/s   |   "
        f"Excluding F200LP_Sersic2D (poor fit, n→0.30): {p50s_ok.max()-p50s_ok.min():.1f} km/s",
        fontsize=13,
    )
    out = os.path.join(FIG_DIR, "nb07c_ishape_sweep.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved → {out}")

    # Save table
    out_npz = "results/sigma_e_ishape_sweep.npz"
    np.savez(
        out_npz,
        ishape=np.array([s for s, _ in have]),
        group=np.array([GROUP_OF[s] for s, _ in have]),
        p16=np.array([i["p16"] for _, i in have]),
        p50=np.array([i["p50"] for _, i in have]),
        p84=np.array([i["p84"] for _, i in have]),
        n_kept=np.array([i["n_kept"] or -1 for _, i in have]),
        delta_headline=np.array([i["p50"] - head_p50 for _, i in have]),
        headline_p50=head_p50,
    )
    print(f"Saved → {out_npz}")


if __name__ == "__main__":
    main()
