"""Summary figure for nb08: aperture vs Sersic-total photometry → M★ comparison.

Two panels:
  A. Photometry comparison — F_lambda vs wavelength for both choices
  B. log(M★) posterior overlay (aperture from nb02, Sersic from nb08)
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)


def main():
    s = np.load("results/bagpipes_sersic_refit.npz", allow_pickle=True)
    p = np.load("results/sersic_total_photometry.npz", allow_pickle=True)
    ap = np.load("results/bagpipes_sed_results.npz", allow_pickle=True)

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.30)

    # Panel A: photometry comparison
    axA = fig.add_subplot(gs[0, 0])
    pivot = s["pivot_AA"]
    flam_ap = s["flam_aperture"]
    flam_se = s["flam_sersic"]
    axA.plot(pivot / 1e4, flam_ap, "o-", color="C0", ms=10, lw=1.5,
             label="aperture (nb02 headline)")
    axA.plot(pivot / 1e4, flam_se, "s-", color="C3", ms=10, lw=1.5,
             label="Sersic-total (nb08)")
    for i, name in enumerate(s["filter_names"]):
        d = s["mag_sersic"][i] - s["mag_aperture"][i]
        axA.annotate(f"{str(name)}\nΔ={d:+.2f}",
                     (pivot[i] / 1e4, flam_se[i]),
                     textcoords="offset points", xytext=(8, 5), fontsize=9)
    axA.set_yscale("log")
    axA.set_xlabel("Wavelength [μm]")
    axA.set_ylabel(r"F$_\lambda$ [erg/s/cm²/Å]")
    axA.set_title("AGEL0206 deflector photometry — aperture vs Sersic-total")
    axA.legend(loc="lower right")
    axA.grid(alpha=0.3)

    # Panel B: M★ posterior overlay
    axB = fig.add_subplot(gs[0, 1])
    M_ap = ap["stellar_mass"]
    M_se = s["log_M_sersic_samples"]

    bins = np.linspace(min(M_ap.min(), M_se.min()) - 0.05,
                       max(M_ap.max(), M_se.max()) + 0.05, 60)
    axB.hist(M_ap, bins=bins, density=True, alpha=0.55, color="C0",
             label=f"aperture (nb02): {np.median(M_ap):.3f} "
                   f"−{np.median(M_ap)-np.percentile(M_ap, 16):.3f}/"
                   f"+{np.percentile(M_ap, 84)-np.median(M_ap):.3f}")
    axB.hist(M_se, bins=bins, density=True, alpha=0.55, color="C3",
             label=f"Sersic-total (nb08): {np.median(M_se):.3f} "
                   f"−{np.median(M_se)-np.percentile(M_se, 16):.3f}/"
                   f"+{np.percentile(M_se, 84)-np.median(M_se):.3f}")
    axB.axvline(np.median(M_ap), color="C0", lw=1.5, ls="--")
    axB.axvline(np.median(M_se), color="C3", lw=1.5, ls="--")

    delta = np.median(M_se) - np.median(M_ap)
    axB.set_xlabel(r"log$_{10}$(M$_\star$ / M$_\odot$)")
    axB.set_ylabel("PDF")
    axB.set_title(f"Bagpipes M★ posterior — Δlog M = {delta:+.3f} dex (×{10**delta:.2f})")
    axB.legend(fontsize=10, loc="upper right")
    axB.grid(alpha=0.3)

    fig.suptitle(
        "nb08 — Aperture systematic on stellar mass via 2D Sersic-total photometry",
        fontsize=13, y=1.02,
    )

    out = "results/figures/nb08_sersic_vs_aperture.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
