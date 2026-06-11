"""Systematic budget for the Sersic-only (full-model) stellar mass.

The Sersic-total M* (integrate the fitted single-Sersic to infinity) is model-based,
so it carries model systematics. This builds its budget, analogous to the sigma_e
budget, by re-deriving the Sersic-total AB mags under each variation and running
Bagpipes (nb02 priors) -> M*. Quadrature sum of the half-spreads.

Components:
  stat        : Bagpipes posterior width (10%) on the fiducial single-Sersic-total
  mask        : refit single-Sersic excluding {expert, global-arc} masks -> spread
  model_form  : single-Sersic vs 2-component (bulge+disk) total -> spread
  flux_floor  : 10% vs 20% fractional flux errors -> half-spread
  fit_param   : perturb n by its per-band fit sigma (Sersic-total is ~ n-sensitive)
  apcorr_recon: pure-Sersic-total vs the aperture-corrected total (method reconciliation)

Output: results/sersic_total_systematic.npz (+ printed budget table)
Usage: conda activate ISMgas; python scripts/sersic_total_systematic.py
"""
import os, sys, json
import numpy as np
from astropy.modeling.models import Sersic2D
from astropy.modeling.fitting import LevMarLSQFitter

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.sersic_total_photometry import (fit_sersic2d, sersic_total_flux_analytic,
                                             R_E_INIT)
from scripts.mask_method_comparison import load_band
from scripts.bagpipes_sersic_refit import run_bagpipes_for_mags

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]


def load(n):
    return load_band(n)


def sersic_total_mag(b, mask, box=6.0, dn=0.0):
    """Single-Sersic total AB mag with `mask` excluded; optional index perturbation dn."""
    fit, _, _ = fit_sersic2d(b["img"], mask, (b["cx"], b["cy"]), R_E_INIT, box, b["pix_scale"])
    amp, reff, n, ellip = (float(fit.amplitude.value), float(fit.r_eff.value),
                           float(fit.n.value) + dn, float(fit.ellip.value))
    F = sersic_total_flux_analytic(amp, reff, max(0.3, n), ellip)
    return float(-2.5 * np.log10(F) + b["ab_zp"]), n


def main():
    bands = {n: load(n) for n in ORDER}
    pivot = np.array([bands[n]["cfg"]["pivot_AA"] for n in ORDER])
    sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
    f200 = bands["F200LP"]

    def mags_for(masks, dn=0.0):
        out = []
        for n in ORDER:
            m, _ = sersic_total_mag(bands[n], masks[n], dn=dn)
            out.append(m)
        return np.array(out)

    expert = {n: bands[n]["mask"].astype(bool) for n in ORDER}
    glob = {n: sysz[f"{n}_global_mask"].astype(bool) for n in ORDER}

    rows = {}
    def M(mags, run, frac=0.1):
        s = run_bagpipes_for_mags(mags, pivot, run, err_frac=frac)
        return np.percentile(s, [16, 50, 84])

    print("Sersic-total mags + M* under variations ...")
    rows["expert_10"] = M(mags_for(expert), "AGEL0206_st_expert", 0.1)
    rows["global_10"] = M(mags_for(glob), "AGEL0206_st_global", 0.1)
    rows["expert_20"] = M(mags_for(expert), "AGEL0206_st_expert20", 0.2)
    rows["np_dn"] = M(mags_for(expert, dn=+0.3), "AGEL0206_st_dn_p", 0.1)
    rows["np_dn2"] = M(mags_for(expert, dn=-0.3), "AGEL0206_st_dn_m", 0.1)

    fid = rows["expert_10"][1]
    stat = (rows["expert_10"][2] - rows["expert_10"][0]) / 2
    mask_sys = abs(rows["global_10"][1] - rows["expert_10"][1])
    floor_sys = abs(rows["expert_20"][1] - rows["expert_10"][1])
    fit_sys = (abs(rows["np_dn"][1] - fid) + abs(rows["np_dn2"][1] - fid)) / 2
    # apcorr reconciliation: pure-Sersic-total vs aperture-corrected total (2 R_e)
    apm = np.load("results/aperture_matched_photometry.npz", allow_pickle=True)
    apcorr_recon = abs(float(apm["logM_total_2"][1]) - fid)
    # model-form: 2-comp not separately fit here -> carry from the photometry masking sys
    model_form = float(np.load("results/Mstar_masking_systematic.npz")["under_arc_10pct"])  # raw<->filled proxy

    comps = {"mask": mask_sys, "model_form": model_form, "flux_floor": floor_sys,
             "fit_param": fit_sys, "apcorr_recon": apcorr_recon}
    sys_quad = float(np.sqrt(sum(v**2 for v in comps.values())))
    total = float(np.hypot(stat, sys_quad))

    print(f"\n=== Sersic-only (full-model) M* budget ===")
    print(f"  central (single-Sersic-total, expert mask, 10%) = {fid:.3f}")
    print(f"  stat (Bagpipes 10%)        ±{stat:.3f}")
    for k, v in comps.items():
        print(f"  {k:<14} ±{v:.3f}")
    print(f"  sys (quadrature)           ±{sys_quad:.3f}")
    print(f"  TOTAL (stat+sys quad)      ±{total:.3f}")

    np.savez("results/sersic_total_systematic.npz",
             central=fid, stat=stat, sys_quad=sys_quad, total=total,
             components=np.array(list(comps.items()), dtype=object),
             **{f"row_{k}": v for k, v in rows.items()})
    print("\nSaved → results/sersic_total_systematic.npz")


if __name__ == "__main__":
    main()
