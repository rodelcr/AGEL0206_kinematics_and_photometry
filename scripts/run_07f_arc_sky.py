#!/usr/bin/env python3
"""
07f — Arc-sky-template σ_e at R<R_e (no-mask + arc-as-sky).

Tests whether the −16 km/s no-mask vs masked Δ from nb09 is purely
continuum dilution by ppxf-fitting the no-mask R<R_e aperture spectrum
with the bright outer-arc spectrum passed in as a `sky` template (free-
amplitude additive component, NOT convolved with the deflector LOSVD).

Three points of comparison (R<R_e, masked headline pipeline):
    σ_A (masked)             from results/final_sigma_e_paper.npz
    σ_B (no-mask)            from results/final_sigma_e_paper.npz
    σ_D (no-mask + arc-sky)  this script

Recovery fraction = (σ_D − σ_B) / (σ_A − σ_B). 1.0 ⇒ pure dilution; 0.0
⇒ none of the Δ is dilution; intermediate ⇒ mixed mechanism.

Arc template: I-weighted average of IFU spaxels that are
    (a) inside the reprojected F200LP arc mask AND
    (b) beyond R_safe = 3 R_e / 4
so they contain mostly arc light with minimal direct deflector contribution
(seeing-PSF leak from the bright deflector centre is the only residual
contaminant — see TESTS_AND_DIAGNOSTICS.md §E7 caveat list).

Caches: results/arc_sky_07f/{sps}_N{n}.npz
Summary: results/arc_sky_07f.npz
"""
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path
from time import perf_counter as clock

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.bootstrap_ppxf import (  # noqa: E402
    SPS_NATIVE_FRAME,
    setup_ppxf_inputs_from_spectrum,
    compute_local_residual_scaling,
    C_KMS,
)
from scripts.final_sigma_e import (  # noqa: E402
    SPS_LIBS,
    DEGREES,
    BOOT_SEED,
    WINDOW,
    Z_SYSTEMIC,
    KPC_PER_ARCSEC,
    load_setup,
    extract_aperture_spectrum,
    pool_sps,
)


CACHE_DIR = REPO / "results" / "arc_sky_07f"
SUMMARY = REPO / "results" / "arc_sky_07f.npz"


# ─────────────────────────────────────────────────────────────────────────────
def build_arc_template(state, r_safe_frac=0.75):
    """Build the arc template spectrum.

    Spaxels inside the reprojected F200LP arc mask AND beyond
    R_safe = r_safe_frac × R_e are I-weighted averaged using the IFU band
    map as the weight. This matches the deflector aperture extraction
    convention for consistency.
    """
    r_safe = r_safe_frac * state["R_E"]
    arc_sel = state["arc_spax_mask"] & (state["r_spax"] > r_safe)
    n = int(arc_sel.sum())
    if n == 0:
        raise RuntimeError("No outer-arc spaxels — adjust r_safe_frac")
    I = np.clip(state["ifu_band"], 0, None)
    w = I[arc_sel]
    w_norm = w / max(w.sum(), 1e-30)
    arc_flux = np.sum(state["cube"][:, arc_sel] * w_norm[None, :], axis=1)
    sn = float(
        np.median(arc_flux[state["band_mask"]])
        / np.median(state["noise_sky"][state["band_mask"]])
    )
    print(f"  arc template: N_spax={n}  R>{r_safe:.3f}\"  S/N(band)={sn:.1f}")
    return arc_flux, n, sn, r_safe


# ─────────────────────────────────────────────────────────────────────────────
def run_one_sps(state, arc_flux, sps_name, n_bootstrap, force):
    """Run ppxf at no-mask R<R_e with arc-as-sky for one SPS library."""
    cache = CACHE_DIR / f"{sps_name}_N{n_bootstrap}.npz"
    if cache.exists() and not force:
        d = dict(np.load(cache, allow_pickle=True))
        for k in list(d.keys()):
            if d[k].shape == ():
                d[k] = d[k].item()
        return d

    R_E = state["R_E"]
    flux_apt, noise_apt, n_apt, sn_apt = extract_aperture_spectrum(
        state, R_E, mask_weight=1.0,
    )
    inputs = setup_ppxf_inputs_from_spectrum(
        flux_apt, noise_apt, state["hdr"], sps_name=sps_name, z=Z_SYSTEMIC,
        verbose=False, frame_galaxy="auto",
    )
    arc_inputs = setup_ppxf_inputs_from_spectrum(
        arc_flux, noise_apt, state["hdr"], sps_name=sps_name, z=Z_SYSTEMIC,
        verbose=False, frame_galaxy="auto",
    )
    sky = np.atleast_2d(arc_inputs["galaxy"]).T  # (n_pix, 1)

    n_pix = len(inputs["galaxy"])
    n_deg = len(DEGREES)
    V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
    sky_amp = np.zeros(n_deg)
    V_boot = np.full((n_deg, n_bootstrap), np.nan)
    sig_boot = np.full((n_deg, n_bootstrap), np.nan)
    bf_arr = np.zeros((n_deg, n_pix))
    t0 = clock()

    for i, deg in enumerate(DEGREES):
        pp = ppxf(
            inputs["sps"].templates, inputs["galaxy"], inputs["noise"],
            inputs["velscale"], inputs["start"],
            goodpixels=inputs["goodpixels"], plot=False, moments=2, trig=False,
            degree=int(deg), mdegree=0,
            lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"],
            sky=sky, quiet=True,
        )
        V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
        bf_arr[i] = pp.bestfit
        # ppxf flattens the (n_pix, n_age, n_Z) template grid internally
        # then appends sky-template columns. The single arc-sky weight is
        # the last element of pp.weights (we pass sky.shape[1] == 1).
        sky_amp[i] = float(pp.weights[-1])

        # Bootstrap with sky template held fixed
        residuals = inputs["galaxy"] - pp.bestfit
        scale = compute_local_residual_scaling(residuals, window=WINDOW)
        scaled = residuals * scale
        rng = np.random.default_rng(BOOT_SEED + 7000 + 100 * i
                                    + (abs(hash(sps_name)) % 100))
        for b in range(n_bootstrap):
            signs = rng.choice([-1, 1], size=n_pix)
            galaxy_boot = pp.bestfit + scaled * signs
            try:
                pp_b = ppxf(
                    inputs["sps"].templates, galaxy_boot, inputs["noise"],
                    inputs["velscale"], inputs["start"],
                    goodpixels=inputs["goodpixels"], plot=False, moments=2, trig=False,
                    degree=int(deg), mdegree=0,
                    lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"],
                    sky=sky, quiet=True,
                )
                V_boot[i, b] = pp_b.sol[0]
                sig_boot[i, b] = pp_b.sol[1]
            except Exception:
                continue

    elapsed = clock() - t0
    out = dict(
        sps_name=sps_name, r_max=float(R_E),
        n_spax=int(n_apt), sn_band=float(sn_apt),
        degrees=np.asarray(DEGREES),
        V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
        sky_amp=sky_amp,
        V_boot=V_boot, sig_boot=sig_boot,
        frame_galaxy=inputs["frame_galaxy"],
        velscale=float(inputs["velscale"]),
        n_pix=int(n_pix), n_bootstrap=int(n_bootstrap),
        elapsed_s=float(elapsed),
        galaxy=inputs["galaxy"], noise=inputs["noise"],
        arc_template=arc_inputs["galaxy"],
        lam_gal_rest=inputs["lam_gal_rest"], best_fit=bf_arr,
    )
    np.savez(cache, **out)
    print(f"    [arc-sky    ] R_e/{sps_name:6s}  "
          f"σ={sig_orig.min():.0f}-{sig_orig.max():.0f} "
          f"V={V_orig.mean():+.1f}  α_arc={sky_amp.mean():.3f}  "
          f"frame={inputs['frame_galaxy']:7s}  "
          f"{elapsed:.0f}s  ({n_apt} spax, S/N={sn_apt:.1f})")
    return out


# ─────────────────────────────────────────────────────────────────────────────
def main(n_bootstrap, force):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n{'═'*70}\n07f ARC-SKY σ_e — N_BOOTSTRAP={n_bootstrap}, force={force}\n{'═'*70}")
    print(f"  cache → {CACHE_DIR.relative_to(REPO)}")

    state = load_setup()

    print("\n§B  Build arc template (outer-arc spaxels, R > 3R_e/4, in F200 mask)")
    print("─" * 70)
    arc_flux, n_arc, sn_arc, r_safe = build_arc_template(state)

    print(f"\n§C  ppxf at R<R_e with arc-as-sky, per SPS")
    print("─" * 70)
    per_sps = {}
    for sps in SPS_LIBS:
        per_sps[sps] = run_one_sps(state, arc_flux, sps, n_bootstrap, force)

    pool = pool_sps(per_sps)

    # Comparison anchors from nb09 npz (Track A and Track B at R<R_e)
    nb09 = np.load(REPO / "results" / "final_sigma_e_paper.npz", allow_pickle=True)
    sigma_A = float(nb09["sigma_p50"][-1])           # masked headline
    sigma_B = float(nb09["nomask_sigma_p50"][-1])    # no-mask
    sigma_D = float(pool["sigma_p50"])               # no-mask + arc-sky

    delta_AB = sigma_A - sigma_B
    delta_DB = sigma_D - sigma_B
    recovery = (delta_DB / delta_AB) if abs(delta_AB) > 1e-6 else float("nan")

    print(f"\n§D  Comparison at R<R_e = {state['R_E']:.3f}\"")
    print("─" * 70)
    print(f"  σ_A (masked headline)        = {sigma_A:7.2f} km/s  [+0.0 baseline]")
    print(f"  σ_B (no-mask)                = {sigma_B:7.2f} km/s  "
          f"[Δ from A = {sigma_B - sigma_A:+.2f}]")
    print(f"  σ_D (no-mask + arc-sky, 07f) = {sigma_D:7.2f} km/s  "
          f"[Δ from A = {sigma_D - sigma_A:+.2f},  Δ from B = {delta_DB:+.2f}]")
    print(f"\n  Recovery fraction (σ_D − σ_B)/(σ_A − σ_B) = {recovery:.2%}")
    if recovery > 0.9:
        verdict = "≈ 1: arc continuum dilution explains the full no-mask Δ"
    elif recovery > 0.5:
        verdict = "≈ majority: dilution dominant, modest residual mechanism"
    elif recovery > 0.1:
        verdict = "partial: dilution + secondary mechanism comparable"
    else:
        verdict = "≈ 0: arc-sky did NOT recover the masked σ — secondary mechanism dominant"
    print(f"  Verdict: {verdict}")

    np.savez(
        SUMMARY,
        z_systemic=Z_SYSTEMIC,
        R_E=state["R_E"],
        r_safe=r_safe,
        n_arc_spax=n_arc, sn_arc=sn_arc,
        kpc_per_arcsec=KPC_PER_ARCSEC,
        n_bootstrap=int(n_bootstrap),
        sps_libs=np.array(SPS_LIBS),
        sigma_p16=pool["sigma_p16"], sigma_p50=pool["sigma_p50"], sigma_p84=pool["sigma_p84"],
        V_p16=pool["V_p16"], V_p50=pool["V_p50"], V_p84=pool["V_p84"],
        sigma_samples=pool["sigma_samples"],
        per_sps_sigma_p50=np.array([per_sps[s]["sig_orig"].mean() for s in SPS_LIBS]),
        per_sps_sky_amp=np.array([per_sps[s]["sky_amp"].mean() for s in SPS_LIBS]),
        sigma_A_masked=sigma_A,
        sigma_B_nomask=sigma_B,
        sigma_D_arc_sky=sigma_D,
        delta_AB=delta_AB, delta_DB=delta_DB,
        recovery_fraction=recovery,
        verdict=verdict,
        arc_template_galaxy=per_sps[SPS_LIBS[0]]["arc_template"],
        lam_gal_rest=per_sps[SPS_LIBS[0]]["lam_gal_rest"],
    )
    print(f"\nSaved → {SUMMARY.relative_to(REPO)}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n_bootstrap", type=int, default=50)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    main(args.n_bootstrap, args.force)
