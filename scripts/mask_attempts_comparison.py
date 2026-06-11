"""Side-by-side comparison of the arc-mask ATTEMPTS, to document why we settled on
the single-Sersic + color-gated mask.

Three attempts (deflector light model + IR-extension growth rule):
  A. 2-component (bulge+disk) Sersic, residual-grown IR extension (NO color gate)
     -> results/photometry_systematics_2comp_backup.npz
  B. single-component Sersic with a BROKEN fit (amplitude->0) + color gate  [DISCARDED
     artifact; not plotted — model was ~zero so masks were color-only]
  C. single-component Sersic (fixed amplitude init) + color-gated IR extension  [FINAL]
     -> results/photometry_systematics.npz

For each band, overlay the A (red) and C (lime) per-band mask contours on the image,
annotate mask sizes. Saves results/figures/mask_attempts_comparison.png.

Usage: conda activate ISMgas; python scripts/mask_attempts_comparison.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.mask_method_comparison import load_band

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]


def main():
    A = np.load("results/photometry_systematics_2comp_backup.npz", allow_pickle=True)
    C = np.load("results/photometry_systematics.npz", allow_pickle=True)
    bands = {n: load_band(n) for n in ORDER}

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.4))
    for ax, n in zip(axes, ORDER):
        b = bands[n]
        ps = b["pix_scale"]
        half = int(np.ceil(6.0 / ps))
        cy, cx = int(b["cy"]), int(b["cx"])
        sl = (slice(max(0, cy - half), cy + half), slice(max(0, cx - half), cx + half))
        img = np.nan_to_num(b["img"])[sl]
        vmax = float(np.nanpercentile(np.abs(img), 99))
        ax.imshow(img, origin="lower", cmap="magma", vmin=0, vmax=vmax)
        mA = A[f"{n}_perband_mask"][sl].astype(float)
        mC = C[f"{n}_perband_mask"][sl].astype(float)
        ax.contour(mA, levels=[0.5], colors="red", linewidths=1.1)
        ax.contour(mC, levels=[0.5], colors="lime", linewidths=1.1)
        nA = int(A[f"{n}_perband_mask"].sum())
        nC = int(C[f"{n}_perband_mask"].sum())
        ax.set_title(f"{n}\nred = 2-comp ({nA}px)   lime = 1-Sersic+color ({nC}px)",
                     fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Arc-mask attempts (per-band): 2-component residual-grown (red) vs "
                 "single-Sersic + color-gated (lime, FINAL)", fontsize=13, y=1.02)
    plt.tight_layout()
    out = "results/figures/mask_attempts_comparison.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Saved → {out}")

    # numeric summary
    print(f"\n{'band':<8}{'2comp px':>10}{'1Ser px':>10}{'raw 2comp':>11}{'raw 1Ser':>11}"
          f"{'fill 2comp':>12}{'fill 1Ser':>12}")
    for i, n in enumerate(ORDER):
        print(f"{n:<8}{int(A[f'{n}_perband_mask'].sum()):>10}"
              f"{int(C[f'{n}_perband_mask'].sum()):>10}"
              f"{A['perband_raw'][i]:>11.3f}{C['perband_raw'][i]:>11.3f}"
              f"{A['perband_filled'][i]-A['perband_raw'][i]:>+12.3f}"
              f"{C['perband_filled'][i]-C['perband_raw'][i]:>+12.3f}")


if __name__ == "__main__":
    main()
