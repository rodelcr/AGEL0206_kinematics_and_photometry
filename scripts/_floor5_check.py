"""5% vs 10% flux-floor sensitivity on the matched-aperture estimators (2 R_e).
Headline stays 10%; 5% is the robustness check."""
import os, sys, numpy as np
sys.path.insert(0, os.getcwd())
from scripts.bagpipes_sersic_refit import run_bagpipes_for_mags
d = np.load("results/aperture_matched_photometry.npz", allow_pickle=True)
pivot = d["pivot"]
out = {}
for kind in ("raw", "raw_apcorr", "filled", "total"):
    mags = d[f"mag_{kind}_2"]
    s = run_bagpipes_for_mags(mags, pivot, f"AGEL0206_aper_{kind}_2Re_5pct", err_frac=0.05)
    p = np.percentile(s, [16, 50, 84])
    out[kind] = p
    p10 = d[f"logM_{kind}_2"]
    print(f"  {kind:<11} 5%: {p[1]:.3f} (-{p[1]-p[0]:.3f}/+{p[2]-p[1]:.3f})   "
          f"10%: {p10[1]:.3f} (-{p10[1]-p10[0]:.3f}/+{p10[2]-p10[1]:.3f})   Δcentral={p[1]-p10[1]:+.3f}")
np.savez("results/aperture_floor5_check.npz", **{f"logM_{k}_5pct": v for k, v in out.items()})
print("Saved → results/aperture_floor5_check.npz")
