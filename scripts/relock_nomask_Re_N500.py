"""Run the 3 missing no-mask Re fits at N=500 using the parallel bootstrap runner.

After the original N=500 production hit a 4-hour numerical pathology on the
no-mask Re_2/xsl fit, the no-mask Re track was cached only at N=50. This
driver fills in the gap using `bootstrap_ppxf_parallel.run_bootstrap_single_degree_parallel`
(joblib pool with BLAS=1 per worker), which runs ~5× faster.

Cache writes match the schema of `scripts.final_sigma_e.run_aperture_sps`,
so on the next `final_sigma_e.py` run the fallback chain auto-promotes
to N=500.
"""
import os
import sys
from pathlib import Path
from time import perf_counter as clock

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from ppxf.ppxf import ppxf  # noqa: E402

from scripts.final_sigma_e import (   # noqa: E402
    load_setup, extract_aperture_spectrum,
    CACHE_DIR, APERTURE_LABELS, APERTURE_FRACS, SPS_LIBS,
    Z_SYSTEMIC, DEGREES, BOOT_SEED, WINDOW,
)
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum  # noqa: E402
from scripts.bootstrap_ppxf_parallel import run_bootstrap_single_degree_parallel  # noqa: E402

N_BOOTSTRAP = 500
LABEL = "Re"
MASK_ON = False
N_JOBS = 8


def main():
    state = load_setup()
    r_max = state["R_E"]
    flux, noise, n_kept, sn_band = extract_aperture_spectrum(
        state, r_max, mask_on=MASK_ON,
    )
    print(f"\nRelock target: {LABEL}/nomask  R<{r_max:.3f}\"  ({n_kept} spax, S/N={sn_band:.1f})")

    label_idx = APERTURE_LABELS.index(LABEL)

    for sps_name in SPS_LIBS:
        cache = CACHE_DIR / f"{LABEL}_{sps_name}_N{N_BOOTSTRAP}_nomask.npz"
        if cache.exists():
            print(f"  Skip {sps_name}: {cache.name} exists")
            continue
        print(f"\n  === {LABEL}/{sps_name}/nomask  N={N_BOOTSTRAP}  parallel ===")
        inputs = setup_ppxf_inputs_from_spectrum(
            flux, noise, state["hdr"], sps_name=sps_name, z=Z_SYSTEMIC,
            verbose=False, frame_galaxy="auto",
        )
        n_pix = len(inputs["galaxy"])
        n_deg = len(DEGREES)
        V_orig = np.zeros(n_deg); sig_orig = np.zeros(n_deg); chi2_orig = np.zeros(n_deg)
        V_boot = np.full((n_deg, N_BOOTSTRAP), np.nan)
        sig_boot = np.full((n_deg, N_BOOTSTRAP), np.nan)
        bf_arr = np.zeros((n_deg, n_pix))
        t0 = clock()
        for i, deg in enumerate(DEGREES):
            pp = ppxf(
                inputs["sps"].templates, inputs["galaxy"], inputs["noise"],
                inputs["velscale"], inputs["start"],
                goodpixels=inputs["goodpixels"], plot=False, moments=2, trig=False,
                degree=int(deg), mdegree=0,
                lam=inputs["lam_gal_rest"], lam_temp=inputs["lam_temp"], quiet=True,
            )
            V_orig[i], sig_orig[i], chi2_orig[i] = pp.sol[0], pp.sol[1], pp.chi2
            bf_arr[i] = pp.bestfit
            rb = run_bootstrap_single_degree_parallel(
                inputs, degree=int(deg), best_fit_spectrum=pp.bestfit,
                n_bootstrap=N_BOOTSTRAP, window=WINDOW,
                seed=BOOT_SEED + 1000 * label_idx + i, n_jobs=N_JOBS,
            )
            V_boot[i] = rb["V_samples"]
            sig_boot[i] = rb["sigma_samples"]
        elapsed = clock() - t0
        np.savez(
            cache,
            label=LABEL, sps_name=sps_name, r_max=float(r_max),
            mask_on=MASK_ON, n_spax=int(n_kept), sn_band=float(sn_band),
            degrees=np.asarray(DEGREES),
            V_orig=V_orig, sig_orig=sig_orig, chi2_orig=chi2_orig,
            V_boot=V_boot, sig_boot=sig_boot,
            frame_galaxy=inputs["frame_galaxy"],
            velscale=float(inputs["velscale"]),
            n_pix=int(n_pix), n_bootstrap=int(N_BOOTSTRAP),
            elapsed_s=float(elapsed),
            galaxy=inputs["galaxy"], noise=inputs["noise"],
            lam_gal_rest=inputs["lam_gal_rest"], best_fit=bf_arr,
        )
        print(f"    done in {elapsed:.1f}s ({elapsed/60:.1f} min); "
              f"σ={sig_orig.min():.0f}-{sig_orig.max():.0f} km/s")
    print("\nAll missing N=500 nomask Re fits committed to cache.")


if __name__ == "__main__":
    main()
