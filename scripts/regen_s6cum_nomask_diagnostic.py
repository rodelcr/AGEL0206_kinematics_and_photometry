"""Regenerate nb07c_s6cum_nomask_diagnostic.png with the corrected HST_mean center.

The original cell in nb07c used stale `cx_f, cy_f` variables that drifted
~0.09" from the HST-mean center used by the actual analysis. This standalone
script rebuilds the diagnostic from the production caches using the HST_mean
center (cx_ifu, cy_ifu) — exactly the center the σ_e numbers were computed at.

Output: results/figures/nb07c_s6cum_nomask_diagnostic.png
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.run_isource_sweep import load_setup  # noqa: E402

ANNULAR_DIR = "results/annular_bootstrap_07c"
ANNULAR_DIR_NOMASK = "results/annular_bootstrap_07c_nomask"
RE_TAG = "2p305"
SPS_LIBS = ["fsps", "emiles", "xsl"]
FIG_DIR = "results/figures"


def main():
    # Setup: cube, HST_mean center, ifu_band, arc_spax_mask, R_E
    state = load_setup()
    cube = state["cube"]
    cx, cy = state["cx_ifu"], state["cy_ifu"]
    ifu_band = state["ifu_band"]
    r_spax = state["r_spax"]
    arc_spax_mask = state["arc_spax_mask"]
    R_E = state["R_E"]
    nx, ny = state["nx"], state["ny"]
    pix = 0.30

    # I-weight as used in §6cum
    I_weight = np.clip(ifu_band, 0, None)
    sel_Re = r_spax < R_E
    n_kept_masked = int((sel_Re & ~arc_spax_mask).sum())
    n_kept_nomask = int(sel_Re.sum())

    I_total = float(I_weight[sel_Re].sum())
    I_lost = float(I_weight[sel_Re & arc_spax_mask].sum())
    I_kept_pct = 100 * (1 - I_lost / I_total)

    print(f"At R<R_e: F200 mask removes {(sel_Re & arc_spax_mask).sum()}/{sel_Re.sum()} spaxels "
          f"({100*I_lost/I_total:.1f}% of I-weight)")
    print(f"  → headline keeps {I_kept_pct:.1f}% of IFU-band light at R<R_e")

    # Load galaxy spectra from cache (one SPS — they all share the rest-grid)
    m = np.load(os.path.join(ANNULAR_DIR, f"cumR_{RE_TAG}_fsps.npz"))
    n = np.load(os.path.join(ANNULAR_DIR_NOMASK, f"cumR_{RE_TAG}_fsps.npz"))
    gal_masked = m["galaxy"]
    gal_nomask = n["galaxy"]
    lam_rest = m["lam_gal_rest"]

    # Combined-SPS posteriors at R_e
    pool_m, pool_n = [], []
    for sps in SPS_LIBS:
        d = np.load(os.path.join(ANNULAR_DIR, f"cumR_{RE_TAG}_{sps}.npz"))
        pool_m.append(d["sig_boot"].ravel())
        d = np.load(os.path.join(ANNULAR_DIR_NOMASK, f"cumR_{RE_TAG}_{sps}.npz"))
        pool_n.append(d["sig_boot"].ravel())
    pool_m = np.concatenate(pool_m); pool_m = pool_m[np.isfinite(pool_m)]
    pool_n = np.concatenate(pool_n); pool_n = pool_n[np.isfinite(pool_n)]
    m16, m50, m84 = np.percentile(pool_m, [16, 50, 84])
    n16, n50, n84 = np.percentile(pool_n, [16, 50, 84])

    # Centre in arcsec on the imshow extent grid (extent uses [0, nx*pix])
    cx_as = (cx + 0.5) * pix
    cy_as = (cy + 0.5) * pix

    # ─── Figure ───
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.2], width_ratios=[1, 1, 1.2],
                          hspace=0.35, wspace=0.30)

    # Panel A: F200LP-masked I-weight at R<R_e
    axA = fig.add_subplot(gs[0, 0])
    img_A = np.where(sel_Re & ~arc_spax_mask, I_weight, np.nan)
    vmin = max(np.nanmin(img_A[img_A > 0]), 1e-3)
    vmax = float(np.nanmax(I_weight))
    axA.imshow(
        img_A, origin="lower", cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=[0, nx * pix, 0, ny * pix],
    )
    arc_overlay = np.where(sel_Re & arc_spax_mask, 1.0, np.nan)
    axA.imshow(arc_overlay, origin="lower", cmap="Reds", alpha=0.55,
               extent=[0, nx * pix, 0, ny * pix])
    axA.add_patch(Circle((cx_as, cy_as), R_E, fc="none", ec="cyan", lw=1.5, ls="--",
                         label=f"R_e = {R_E:.2f}\""))
    axA.plot(cx_as, cy_as, "c+", ms=14, mew=2)
    axA.set_xlim(cx_as - 1.5 * R_E, cx_as + 1.5 * R_E)
    axA.set_ylim(cy_as - 1.5 * R_E, cy_as + 1.5 * R_E)
    axA.set_xlabel(r"X (arcsec)")
    axA.set_ylabel(r"Y (arcsec)")
    axA.set_title(
        f"F200LP-masked  (N_kept = {n_kept_masked})\n"
        f"red overlay = arc spaxels excluded\n"
        f"I-weight kept = {I_kept_pct:.1f}%  •  spaxels masked = "
        f"{int((sel_Re & arc_spax_mask).sum())}/{int(sel_Re.sum())} "
        f"({100*(sel_Re & arc_spax_mask).sum()/sel_Re.sum():.1f}%)",
        fontsize=10,
    )

    # Panel B: no mask
    axB = fig.add_subplot(gs[0, 1])
    img_B = np.where(sel_Re, I_weight, np.nan)
    axB.imshow(
        img_B, origin="lower", cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=[0, nx * pix, 0, ny * pix],
    )
    axB.contour(
        np.arange(nx) * pix + 0.5 * pix,
        np.arange(ny) * pix + 0.5 * pix,
        (sel_Re & arc_spax_mask).astype(float),
        levels=[0.5], colors=["#3a5a85"], linewidths=1.5, linestyles="--",
    )
    axB.add_patch(Circle((cx_as, cy_as), R_E, fc="none", ec="cyan", lw=1.5, ls="--"))
    axB.plot(cx_as, cy_as, "c+", ms=14, mew=2)
    axB.set_xlim(cx_as - 1.5 * R_E, cx_as + 1.5 * R_E)
    axB.set_ylim(cy_as - 1.5 * R_E, cy_as + 1.5 * R_E)
    axB.set_xlabel(r"X (arcsec)")
    axB.set_ylabel(r"Y (arcsec)")
    axB.set_title(
        f"No arc mask  (N_kept = {n_kept_nomask})\n"
        f"blue dashed outline = previously-masked region (now admitted)\n"
        f"+{n_kept_nomask - n_kept_masked} arc spaxels admitted",
        fontsize=10,
    )

    # Panel C: σ_e posterior overlay
    axC = fig.add_subplot(gs[0, 2])
    axC.hist(pool_m, bins=60, density=True, alpha=0.55, color="#e63946",
             label=f"F200LP masked: {m50:.0f} −{m50-m16:.0f}/+{m84-m50:.0f}")
    axC.hist(pool_n, bins=60, density=True, alpha=0.55, color="#4a6fa5",
             label=f"no arc mask: {n50:.0f} −{n50-n16:.0f}/+{n84-n50:.0f}")
    axC.axvline(m50, color="#a8253a", lw=1.5)
    axC.axvline(n50, color="#3a5a85", lw=1.5)
    axC.axvline(m50 + (m50 - n50) * 0, color="black", ls=":", lw=0.5)
    axC.set_xlabel(r"$\sigma_e (<R_e)$ (km/s)")
    axC.set_ylabel("PDF")
    axC.legend(fontsize=10, loc="upper right")
    axC.set_title(
        f"σ_e joint posterior\n3 SPS × 15 deg × 500 boot, N={len(pool_m)}\n"
        f"Δσ_e = {n50-m50:+.1f} km/s",
        fontsize=10,
    )
    axC.grid(alpha=0.3)

    # Panel D: aperture spectra overlay (rest-frame)
    axD = fig.add_subplot(gs[1, :])
    axD.step(lam_rest, gal_masked, color="#e63946", lw=1.2,
             label="F200LP-masked aperture (R<R_e)", zorder=3)
    axD.step(lam_rest, gal_nomask, color="#4a6fa5", lw=1.0, alpha=0.9,
             label="No-mask aperture (R<R_e)", zorder=2)
    diff = gal_masked - gal_nomask
    axD.plot(lam_rest, diff, color="gray", lw=0.8, alpha=0.7,
             label="masked − nomask = removed arc-light residual")
    for name, wl in [("Ca K", 3933.66), ("Ca H", 3968.47), ("G-band", 4304.40)]:
        axD.axvline(wl, color="k", ls=":", lw=0.6, alpha=0.4)
        axD.text(wl, 0.02, name, ha="center", va="bottom",
                 fontsize=9, color="k", fontstyle="italic",
                 transform=axD.get_xaxis_transform())
    axD.set_xlabel(r"Rest wavelength ($\AA$)")
    axD.set_ylabel("I-weighted flux (normalized)")
    axD.legend(fontsize=10, loc="upper right")
    axD.grid(alpha=0.3)
    axD.set_title(
        f"I-weighted aperture spectrum at R<R_e — masked vs no-mask  "
        f"(F200 mask removes {100-I_kept_pct:.1f}% of I-weight)",
        fontsize=11,
    )
    axD.set_xlim(lam_rest[0], lam_rest[-1])

    fig.suptitle(
        f"§6cum F200LP arc-mask sensitivity — masked vs no-mask diagnostic\n"
        f"HST_mean center IFU pix=({cx:.2f}, {cy:.2f})  ·  R_e={R_E:.3f}\"  ·  "
        f"F200 raw mask = {arc_spax_mask.sum()} of {nx*ny} spaxels ({100*arc_spax_mask.mean():.2f}%)",
        fontsize=13, y=1.00,
    )
    out = os.path.join(FIG_DIR, "nb07c_s6cum_nomask_diagnostic.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    main()
