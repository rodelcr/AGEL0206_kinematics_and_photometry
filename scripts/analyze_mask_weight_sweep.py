"""5-point mask-weight sweep analysis: σ_e(w) for w ∈ {0.0, 0.25, 0.5, 0.75, 1.0}.

Loads the per-(aperture, SPS, mask_weight) caches written by
`scripts/final_sigma_e.py` and `scripts/soft_mask_track.py`, pools σ
posteriors across SPS libraries (FSPS/EMILES/XSL) per weight, and
characterizes the σ_e(w) curve:

  * Linear fit:   σ(w) = σ_0 + slope × w
  * Quadratic fit: σ(w) = a + b·w + c·w²
  * Curvature ratio: σ(0.5) − [(σ(0) + σ(1)) / 2]
        positive  →  sub-linear (concave-down) → masking is gentle near w=0
        zero      →  exactly linear
        negative  →  super-linear (concave-up)  → bias dominated by a few
                                                  highly-contaminating spaxels

Outputs:
  results/mask_weight_sweep.npz                 — summary + samples per weight
  results/figures/nb09_mask_weight_sweep.png    — σ_e(w) curve

Usage:
    python scripts/analyze_mask_weight_sweep.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)

CACHE_DIR = REPO / "results" / "final_sigma_e_paper"
SPS_LIBS = ("fsps", "emiles", "xsl")
APERTURE_LABELS = ("Re_2", "Re")
N_BOOTSTRAP = 500
WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)


def _suffix(w: float) -> str:
    """Cache filename suffix for a given mask weight (matches final_sigma_e.py
    and soft_mask_track.py conventions)."""
    if w == 0.0:
        return ""
    if w == 1.0:
        return "_nomask"
    s = f"{w:g}".replace(".", "p")
    return f"_softmask_w{s}"


def load_per_sps(label: str, w: float, n: int) -> dict | None:
    """Load all SPS caches for one (label, weight). Returns None if any missing."""
    out = {}
    for sps in SPS_LIBS:
        cache = CACHE_DIR / f"{label}_{sps}_N{n}{_suffix(w)}.npz"
        if not cache.exists():
            return None
        d = dict(np.load(cache, allow_pickle=True))
        for k in list(d.keys()):
            if d[k].shape == ():
                d[k] = d[k].item()
        out[sps] = d
    return out


def pool(per_sps: dict) -> dict:
    """Pool σ across SPS by concatenating bootstrap samples (same convention as
    final_sigma_e.pool_sps)."""
    sigs = []
    Vs = []
    per = {}
    for sps, d in per_sps.items():
        s = d["sig_boot"].ravel()
        V = d["V_boot"]
        V_sys = float(np.nanmedian(V))
        sigs.append(s)
        Vs.append((V - V_sys).ravel())
        per[sps] = dict(
            sigma_p50=float(np.nanpercentile(s, 50)),
            sigma_p16=float(np.nanpercentile(s, 16)),
            sigma_p84=float(np.nanpercentile(s, 84)),
            V_sys=V_sys,
            n_bootstrap=int(d["n_bootstrap"]),
        )
    sigma_pool = np.concatenate(sigs)
    sigma_pool = sigma_pool[np.isfinite(sigma_pool)]
    return dict(
        sigma_samples=sigma_pool,
        sigma_p16=float(np.percentile(sigma_pool, 16)),
        sigma_p50=float(np.percentile(sigma_pool, 50)),
        sigma_p84=float(np.percentile(sigma_pool, 84)),
        per_sps=per,
    )


def fit_curves(weights: np.ndarray, sig: np.ndarray, sig_err: np.ndarray) -> dict:
    """Fit linear and quadratic models to σ(w). Returns coefficients + residuals."""
    w = np.asarray(weights, dtype=float)
    y = np.asarray(sig, dtype=float)
    e = np.asarray(sig_err, dtype=float)
    # Linear: weighted least squares
    W = 1.0 / e**2
    A1 = np.vstack([np.ones_like(w), w]).T
    coef_lin, *_ = np.linalg.lstsq(A1 * np.sqrt(W)[:, None], y * np.sqrt(W), rcond=None)
    yhat_lin = A1 @ coef_lin
    resid_lin = y - yhat_lin
    chi2_lin = float(np.sum((resid_lin / e) ** 2))
    # Quadratic
    A2 = np.vstack([np.ones_like(w), w, w**2]).T
    coef_quad, *_ = np.linalg.lstsq(A2 * np.sqrt(W)[:, None], y * np.sqrt(W), rcond=None)
    yhat_quad = A2 @ coef_quad
    resid_quad = y - yhat_quad
    chi2_quad = float(np.sum((resid_quad / e) ** 2))
    # Curvature: σ(0.5) − midpoint of σ(0) and σ(1)
    if len(w) >= 3:
        i0 = int(np.argmin(np.abs(w - 0.0)))
        i1 = int(np.argmin(np.abs(w - 1.0)))
        i_mid = int(np.argmin(np.abs(w - 0.5)))
        mid_chord = 0.5 * (y[i0] + y[i1])
        curvature = float(y[i_mid] - mid_chord)
        curvature_err = float(np.sqrt(0.25 * e[i0] ** 2 + 0.25 * e[i1] ** 2 + e[i_mid] ** 2))
    else:
        curvature = float("nan")
        curvature_err = float("nan")
    return dict(
        lin_intercept=float(coef_lin[0]),
        lin_slope=float(coef_lin[1]),
        lin_chi2=chi2_lin,
        quad_a=float(coef_quad[0]),
        quad_b=float(coef_quad[1]),
        quad_c=float(coef_quad[2]),
        quad_chi2=chi2_quad,
        curvature=curvature,
        curvature_err=curvature_err,
    )


def main():
    print("=" * 72)
    print(f"Mask-weight sweep — N={N_BOOTSTRAP}, weights={WEIGHTS}")
    print(f"Cache dir: {CACHE_DIR.relative_to(REPO)}")
    print("=" * 72)

    # Load all (label, weight) combinations; fall back to N=50 if N=500 missing
    sweep = {lab: {} for lab in APERTURE_LABELS}
    for lab in APERTURE_LABELS:
        print(f"\n[{lab}]")
        for w in WEIGHTS:
            per_sps = load_per_sps(lab, w, N_BOOTSTRAP)
            n_used = N_BOOTSTRAP
            if per_sps is None:
                per_sps = load_per_sps(lab, w, 50)
                n_used = 50
            if per_sps is None:
                print(f"  w={w:>5.2f}  MISSING (skip)")
                continue
            pool_d = pool(per_sps)
            sweep[lab][w] = dict(pool=pool_d, n_bootstrap=n_used, per_sps=per_sps)
            sigma_p50 = pool_d["sigma_p50"]
            sigma_p16 = pool_d["sigma_p16"]
            sigma_p84 = pool_d["sigma_p84"]
            half_err = 0.5 * ((sigma_p84 - sigma_p50) + (sigma_p50 - sigma_p16))
            sps_str = "  ".join(f"{s}:{pool_d['per_sps'][s]['sigma_p50']:.1f}"
                                for s in SPS_LIBS)
            print(f"  w={w:>5.2f}  σ={sigma_p50:6.2f} "
                  f"(-{sigma_p50-sigma_p16:.2f}/+{sigma_p84-sigma_p50:.2f})  "
                  f"⟨err⟩={half_err:.2f}  N={n_used}  per-SPS [{sps_str}]")

    # Build sweep arrays per aperture
    summary = {}
    for lab in APERTURE_LABELS:
        weights_loaded = sorted(sweep[lab].keys())
        if len(weights_loaded) < 3:
            print(f"\n[{lab}] only {len(weights_loaded)} weights loaded; need ≥3 for fits")
            continue
        ws = np.asarray(weights_loaded)
        p50 = np.array([sweep[lab][w]["pool"]["sigma_p50"] for w in ws])
        p16 = np.array([sweep[lab][w]["pool"]["sigma_p16"] for w in ws])
        p84 = np.array([sweep[lab][w]["pool"]["sigma_p84"] for w in ws])
        err = 0.5 * ((p84 - p50) + (p50 - p16))
        fits = fit_curves(ws, p50, err)
        summary[lab] = dict(
            weights=ws, sigma_p50=p50, sigma_p16=p16, sigma_p84=p84, sigma_err=err,
            **fits,
        )
        print(f"\n[{lab}] linearity fits")
        print(f"  Linear:    σ(w) = {fits['lin_intercept']:.2f} + ({fits['lin_slope']:+.2f}) × w     "
              f"χ²={fits['lin_chi2']:.2f}")
        print(f"  Quadratic: σ(w) = {fits['quad_a']:.2f} + ({fits['quad_b']:+.2f})×w + "
              f"({fits['quad_c']:+.2f})×w²   χ²={fits['quad_chi2']:.2f}")
        print(f"  Curvature: σ(0.5) − mid-chord = {fits['curvature']:+.2f} "
              f"± {fits['curvature_err']:.2f} km/s")
        if abs(fits['curvature']) < fits['curvature_err']:
            print(f"             → consistent with linear within 1σ")
        elif fits['curvature'] > 0:
            print(f"             → sub-linear (concave-down): masking gentle, soft is enough")
        else:
            print(f"             → super-linear (concave-up): few bright arc spaxels dominate")

    # Plot σ_e(w) curve, two panels (one per aperture)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    colors = {"Re_2": "C0", "Re": "C3"}
    for ax, lab, lab_disp in zip(axes, APERTURE_LABELS,
                                  [r"$R<R_e/2$", r"$R<R_e$ (headline)"]):
        if lab not in summary:
            ax.set_title(f"{lab_disp} — incomplete sweep", fontsize=12)
            continue
        s = summary[lab]
        ws = s["weights"]; p50 = s["sigma_p50"]; err = s["sigma_err"]
        # Bootstrap envelopes
        lo = p50 - s["sigma_p16"]; hi = s["sigma_p84"] - p50
        ax.errorbar(ws, p50, yerr=[lo, hi], fmt="o", ms=10, lw=1.8, capsize=5,
                    color=colors[lab], label=f"pooled σ_e per weight")
        # Fit curves on a fine grid
        w_grid = np.linspace(0, 1, 101)
        y_lin = s["lin_intercept"] + s["lin_slope"] * w_grid
        y_qua = s["quad_a"] + s["quad_b"] * w_grid + s["quad_c"] * w_grid**2
        ax.plot(w_grid, y_lin, "--", color="gray", lw=1.5, alpha=0.8,
                label=f"linear: σ(w) = {s['lin_intercept']:.0f} + ({s['lin_slope']:+.0f})·w")
        ax.plot(w_grid, y_qua, "-", color="black", lw=1.8, alpha=0.85,
                label=f"quadratic (c={s['quad_c']:+.1f})")
        # Linear chord between σ(0) and σ(1) for visual reference
        if 0.0 in ws and 1.0 in ws:
            i0 = int(np.argmin(np.abs(ws - 0.0)))
            i1 = int(np.argmin(np.abs(ws - 1.0)))
            ax.plot([0, 1], [p50[i0], p50[i1]], ":", color="gray", lw=1, alpha=0.5)

        ax.set_xlabel(r"Arc-spaxel I-weight  $w$", fontsize=12)
        ax.set_ylabel(r"$\sigma_e$ pooled posterior median  [km/s]", fontsize=12)
        ax.set_title(f"{lab_disp}: curvature = {s['curvature']:+.1f} ± "
                     f"{s['curvature_err']:.1f} km/s", fontsize=12)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc="best")
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.text(0.02, 0.02, "w=0: hard mask (drop)\nw=1: no mask (full I-weight)",
                transform=ax.transAxes, fontsize=9, va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))

    fig.suptitle(f"σ_e mask-weight sweep — pooled SPS posteriors (N={N_BOOTSTRAP})",
                 fontsize=13)
    plt.tight_layout()
    out = REPO / "results" / "figures" / "nb09_mask_weight_sweep.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved figure → {out.relative_to(REPO)}")

    # Save summary npz (samples are large; save per-weight pooled samples too)
    save = dict(weights=np.asarray(WEIGHTS), n_bootstrap=N_BOOTSTRAP)
    for lab in APERTURE_LABELS:
        if lab not in summary:
            continue
        s = summary[lab]
        save[f"{lab}_weights_loaded"] = s["weights"]
        save[f"{lab}_sigma_p16"] = s["sigma_p16"]
        save[f"{lab}_sigma_p50"] = s["sigma_p50"]
        save[f"{lab}_sigma_p84"] = s["sigma_p84"]
        save[f"{lab}_sigma_err_mean"] = s["sigma_err"]
        save[f"{lab}_lin_intercept"] = s["lin_intercept"]
        save[f"{lab}_lin_slope"] = s["lin_slope"]
        save[f"{lab}_lin_chi2"] = s["lin_chi2"]
        save[f"{lab}_quad_a"] = s["quad_a"]
        save[f"{lab}_quad_b"] = s["quad_b"]
        save[f"{lab}_quad_c"] = s["quad_c"]
        save[f"{lab}_quad_chi2"] = s["quad_chi2"]
        save[f"{lab}_curvature"] = s["curvature"]
        save[f"{lab}_curvature_err"] = s["curvature_err"]
        # also save per-weight pooled samples (for histogram in nb09)
        for w in s["weights"]:
            samp = sweep[lab][float(w)]["pool"]["sigma_samples"]
            wstr = _suffix(float(w)).replace("_", "") or "w0"
            save[f"{lab}_samples_{wstr}"] = samp
    out_npz = REPO / "results" / "mask_weight_sweep.npz"
    np.savez(out_npz, **save)
    print(f"Saved summary → {out_npz.relative_to(REPO)}")


if __name__ == "__main__":
    main()
