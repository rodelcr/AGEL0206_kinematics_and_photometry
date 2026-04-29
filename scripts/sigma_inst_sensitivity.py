"""σ_inst sensitivity check for the nb09 headline aperture.

Re-fits the headline R<R_e masked I-weighted spectrum with FWHM_inst
inflated by 1.5× to verify that ppxf's LSF deconvolution is robust.
Expected (quadrature deconvolution σ_obs² = σ_int² + σ_inst²):
  σ_galaxy ≈ 268 km/s, σ_inst = 13 → 19.5 km/s gives Δσ ≈ −0.5 km/s.

If the empirical Δσ is < 5 km/s, we declare the LSF subtraction sub-dominant
relative to the SPS systematic and confirm the ±32 km/s budget stands.
"""
import copy
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.final_sigma_e import (   # noqa: E402
    load_setup, extract_aperture_spectrum,
    SPS_LIBS, Z_SYSTEMIC, DEGREES,
)
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402


def fit_one(flux, noise, hdr, sps_name, degree):
    inputs = setup_ppxf_inputs_from_spectrum(
        flux, noise, hdr, sps_name=sps_name, z=Z_SYSTEMIC,
        verbose=False, frame_galaxy="auto",
    )
    pp = ppxf(
        inputs["sps"].templates, inputs["galaxy"], inputs["noise"],
        inputs["velscale"], inputs["start"],
        goodpixels=inputs["goodpixels"], plot=False, moments=2, trig=False,
        degree=int(degree), mdegree=0,
        lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"], quiet=True,
    )
    return float(pp.sol[0]), float(pp.sol[1]), float(pp.chi2)


def main():
    state = load_setup()
    r_max = state["R_E"]
    flux, noise, n_kept, sn_band = extract_aperture_spectrum(state, r_max, mask_on=True)
    print(f"\nHeadline aperture: R<{r_max:.3f}\"  ({n_kept} spax, S/N={sn_band:.1f})")

    hdr_baseline = state["hdr"]
    fwhm_baseline = float(hdr_baseline["DISPSCAL"])
    print(f"DISPSCAL (baseline) = {fwhm_baseline} → FWHM_inst = {2.355*fwhm_baseline:.4f} Å, "
          f"σ_v_inst ≈ 13 km/s")

    # Inflated header: DISPSCAL × 1.5 → FWHM × 1.5 → σ_v_inst × 1.5
    hdr_inflated = copy.deepcopy(hdr_baseline)
    hdr_inflated["DISPSCAL"] = fwhm_baseline * 1.5
    print(f"DISPSCAL (inflated) = {hdr_inflated['DISPSCAL']:.4f} "
          f"→ FWHM_inst × 1.5, σ_v_inst ≈ 19.5 km/s\n")

    # Run at three reference degrees
    degrees_to_test = [DEGREES[0], int(np.median(DEGREES)), DEGREES[-1]]
    print(f"Test degrees: {degrees_to_test}")
    print(f"{'SPS':<8}  {'deg':>4}  {'σ baseline':>12}  {'σ inflated':>12}  "
          f"{'Δσ':>8}  {'V baseline':>12}  {'V inflated':>12}")
    print("-" * 88)
    deltas = []
    for sps in SPS_LIBS:
        for deg in degrees_to_test:
            V_b, sig_b, _ = fit_one(flux, noise, hdr_baseline, sps, deg)
            V_i, sig_i, _ = fit_one(flux, noise, hdr_inflated, sps, deg)
            ds = sig_i - sig_b
            deltas.append(ds)
            print(f"{sps:<8}  {deg:>4d}  {sig_b:>12.2f}  {sig_i:>12.2f}  "
                  f"{ds:>+8.2f}  {V_b:>+12.2f}  {V_i:>+12.2f}")

    deltas = np.asarray(deltas)
    print(f"\nΔσ (σ_inst×1.5 − baseline): "
          f"median={np.median(deltas):+.2f}, "
          f"mean={np.mean(deltas):+.2f}, "
          f"|max|={np.max(np.abs(deltas)):.2f} km/s")
    print(f"\nQuadrature prediction: σ_obs² = σ_int² + σ_inst²")
    print(f"  σ_int = √(268² − 13² + 19.5²) − 268 = {np.sqrt(268**2 - 13**2 + 19.5**2) - 268:+.2f} km/s")
    print(f"  → Empirical |Δσ| ≤ 5 km/s = LSF deconvolution is robust.")

    out = REPO / "results" / "sigma_inst_sensitivity.npz"
    np.savez(out,
             r_max=float(r_max), n_spax=int(n_kept), sn_band=float(sn_band),
             dispscal_baseline=float(fwhm_baseline),
             dispscal_inflated=float(hdr_inflated["DISPSCAL"]),
             deltas_kms=deltas,
             sps=np.array(SPS_LIBS),
             degrees=np.array(degrees_to_test))
    print(f"\nSaved → {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
