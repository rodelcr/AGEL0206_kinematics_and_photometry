"""Stellar-mass masking budget from the (single-Sersic + color-gated) photometry.

Reads the per-band/global raw/filled AB-mag vectors written by
scripts/photometry_systematics.py and runs Bagpipes (nb02 priors) on each, at 10%
and 20% fractional flux errors, then recomputes the masking systematic on log M*.

This makes the M* budget SCRIPT-reproducible (previously the *_Mstar.npz files were
notebook one-offs). Run names carry a `_1ser_` tag so they do NOT collide with the
old 2-component-Sersic caches in pipes/posterior/.

Outputs (consumed by notebook 12 §9):
    results/photometry_systematics_Mstar.npz   (logM percentiles per vector/err)
    results/Mstar_masking_systematic.npz       (masking systematic decomposition)

Usage: conda activate ISMgas; python scripts/mstar_masking_budget.py
"""
import os
import sys

import numpy as np

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

from scripts.bagpipes_sersic_refit import run_bagpipes_for_mags

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
VECTORS = ("perband_raw", "perband_filled", "global_raw", "global_filled")
ERR_FRACS = {"10pct": 0.10, "20pct": 0.20}


def main():
    P = np.load("results/photometry_systematics.npz", allow_pickle=True)
    pivot = P["pivot"]

    out = {"filter_names": np.array(ORDER), "pivot": pivot}
    pcts = {}  # (vector, err) -> [p16,p50,p84]
    for ef, frac in ERR_FRACS.items():
        for vec in VECTORS:
            mags = P[vec]
            run = f"AGEL0206_1ser_{vec}_{ef}"
            print(f"[{ef}] {vec}: mags={np.round(mags,3).tolist()}  run={run}")
            samp = run_bagpipes_for_mags(mags, pivot, run, err_frac=frac)
            p = np.percentile(samp, [16, 50, 84])
            pcts[(vec, ef)] = p
            out[f"logM_{vec}_{ef}"] = p
            out[f"samples_{vec}_{ef}"] = samp

    # headline = empirical raw-central, one-sided-UP (fill-in) systematic; per-band
    for ef in ERR_FRACS:
        c = pcts[("perband_raw", ef)]
        f = pcts[("perband_filled", ef)]
        stat = (c[2] - c[0]) / 2
        up = float(np.hypot(stat, f[1] - c[1]))  # stat (+) under-arc one-sided up
        out[f"headline_{ef}_central"] = c[1]
        out[f"headline_{ef}_stat"] = stat
        out[f"headline_{ef}_sys_up"] = up
        out[f"headline_{ef}_reach"] = f[1]
    np.savez("results/photometry_systematics_Mstar.npz", **out)

    # masking systematic = peak-to-peak/2 over the 4 vectors' medians, per err frac
    msk = {}
    for ef in ERR_FRACS:
        meds = np.array([pcts[(v, ef)][1] for v in VECTORS])
        span = (meds.min(), meds.max())
        masking_sys = (span[1] - span[0]) / 2
        # decomposition
        under_arc = abs(pcts[("perband_filled", ef)][1] - pcts[("perband_raw", ef)][1]) / 2
        maskdef = abs(pcts[("perband_raw", ef)][1] - pcts[("global_raw", ef)][1]) / 2
        msk[f"masking_sys_{ef}"] = masking_sys
        msk[f"span_{ef}"] = np.array(span)
        msk[f"under_arc_{ef}"] = under_arc
        msk[f"maskdef_{ef}"] = maskdef
        for v in VECTORS:
            msk[f"{v}_{ef}"] = pcts[(v, ef)]
        print(f"\n[{ef}] masking_sys = ±{masking_sys:.3f} dex  "
              f"(span {span[0]:.3f}-{span[1]:.3f}; under-arc ±{under_arc:.3f}, "
              f"mask-def ±{maskdef:.3f})")
    np.savez("results/Mstar_masking_systematic.npz", **msk)
    print("\nSaved → results/photometry_systematics_Mstar.npz")
    print("Saved → results/Mstar_masking_systematic.npz")


if __name__ == "__main__":
    main()
