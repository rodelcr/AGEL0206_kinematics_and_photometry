"""Summary figure for nb08: aperture vs Sersic-total photometry → M★ comparison.

Two panels:
  A. SED comparison — for BOTH photometry choices: the best-fit Bagpipes model
     spectrum (line + 16-84 band), the MEASURED photometry (points + errorbars),
     and the filter-convolved MODEL photometry (open squares at the filter pivots).
  B. log(M★) posterior overlay (aperture vs Sersic-total).
All SED/photometry arrays come from results/bagpipes_sersic_refit.npz.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)


def _plot_sed(ax, s, prefix, color, label, mum_lo=0.35, mum_hi=5.0):
    """Overlay one choice: model spectrum band + line, measured photometry (filled,
    errorbars), filter-convolved model photometry (open squares)."""
    wav = s[f"{prefix}_wav_obs"] / 1e4                       # μm (observed)
    sel = (wav > mum_lo) & (wav < mum_hi)
    ax.fill_between(wav[sel], s[f"{prefix}_spec_p16"][sel], s[f"{prefix}_spec_p84"][sel],
                    color=color, alpha=0.18, lw=0)
    ax.plot(wav[sel], s[f"{prefix}_spec_p50"][sel], "-", color=color, lw=1.2, alpha=0.9,
            label=f"{label}: model SED")
    eff = s["eff_wavs"] / 1e4
    # measured photometry (what was fed to the fit) — filled, with errorbars
    ax.errorbar(eff, s[f"{prefix}_flam"], yerr=s[f"{prefix}_flam_err"], fmt="o",
                color=color, ms=8, capsize=3, lw=1.3, mec="k", mew=0.6, zorder=5,
                label=f"{label}: measured phot")
    # filter-convolved model photometry — open squares
    ax.plot(eff, s[f"{prefix}_model_phot_p50"], "s", mfc="none", mec=color, mew=1.6,
            ms=12, zorder=6, label=f"{label}: model phot (filter-conv.)")


def main():
    s = np.load("results/bagpipes_sersic_refit.npz", allow_pickle=True)

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.30)

    # Panel A: SED + measured + filter-convolved model photometry, both choices
    axA = fig.add_subplot(gs[0, 0])
    _plot_sed(axA, s, "ap", "C0", "aperture (masked, not filled)")
    _plot_sed(axA, s, "ser", "C3", "Sersic-total")
    for i, name in enumerate(s["filter_names"]):
        d = s["mag_sersic"][i] - s["mag_aperture"][i]
        axA.annotate(f"{str(name)}\nΔ={d:+.2f}",
                     (s["eff_wavs"][i] / 1e4, s["ser_flam"][i]),
                     textcoords="offset points", xytext=(6, 8), fontsize=8)
    axA.set_yscale("log")
    axA.set_xlabel("Observed wavelength [μm]")
    axA.set_ylabel(r"F$_\lambda$ [erg/s/cm²/Å]")
    axA.set_title("AGEL0206 deflector SED — aperture (empirical) vs Sersic-total")
    axA.legend(loc="lower center", fontsize=8, ncol=2)
    axA.grid(alpha=0.3)

    # Panel B: M★ posterior overlay
    axB = fig.add_subplot(gs[0, 1])
    M_ap = s["log_M_aperture_samples"]
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
