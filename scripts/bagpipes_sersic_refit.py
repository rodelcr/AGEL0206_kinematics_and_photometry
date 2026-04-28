"""Bagpipes refit using Sersic-total photometry → ΔM_star vs aperture headline.

Inputs (from results/sersic_total_photometry.npz):
    mag_aperture, mag_sersic_analytic, mag_aperture_err per band
Same Bagpipes priors as nb02:
    z prior (0.674, 0.676), exponential SFH, Calzetti dust, fractional 10% errors

Output: results/bagpipes_sersic_refit.npz with M_star posterior + comparison
        to nb02 aperture posterior from results/bagpipes_sed_results.npz
"""

import os
import sys
import numpy as np

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO)
os.chdir(REPO)

import astropy.units as u
from astropy.constants import c
import bagpipes as pipes


def mag_to_flam(mags_AB, pivot_AA):
    """AB magnitudes → F_lambda in erg/s/cm²/Å."""
    fnu_units = u.erg / u.s / u.Hz / u.cm**2
    flam_units = u.erg / (u.s * u.cm**2 * u.AA)
    fnu = (np.asarray(mags_AB) * u.ABmag).to(fnu_units)
    pivot = pivot_AA * u.AA
    flam = (fnu * (c / pivot**2)).to(flam_units)
    return flam.value


def main():
    data = np.load("results/sersic_total_photometry.npz", allow_pickle=True)
    filter_names = data["filter_names"]
    pivot = data["pivot_wavelengths_AA"]
    mag_aperture = data["mag_aperture"]
    mag_sersic = data["mag_sersic_analytic"]

    # Order must match the Bagpipes filt_list (F200LP, F140W, F150W2, F322W2)
    expected = ["F200LP", "F140W", "F150W2", "F322W2"]
    idx = [int(np.where(filter_names == n)[0][0]) for n in expected]
    filter_names = filter_names[idx]
    pivot = pivot[idx]
    mag_aperture = mag_aperture[idx]
    mag_sersic = mag_sersic[idx]

    flam_aperture = mag_to_flam(mag_aperture, pivot)
    flam_sersic = mag_to_flam(mag_sersic, pivot)
    print(f"Filter        pivot(Å)   mag_ap    mag_ser   Δmag")
    print("-" * 60)
    for i, name in enumerate(filter_names):
        print(f"{name:<12} {pivot[i]:>8.0f}  {mag_aperture[i]:>7.4f}  "
              f"{mag_sersic[i]:>7.4f}  {mag_sersic[i]-mag_aperture[i]:>+7.4f}")

    # 10% fractional errors (nb02 convention)
    flam_err = 0.1 * flam_sersic

    # Bagpipes resolves filt_list relative to its install_dir, so use absolute paths.
    filt_list = [
        os.path.abspath(p) for p in [
            "HST_WFC3_UVIS1.F200LP.dat",
            "HST_WFC3_IR.F140W.dat",
            "JWST_NIRCam.F150W2.dat",
            "JWST_NIRCam.F322W2.dat",
        ]
    ]

    # nb02 fit_instructions
    exp = {
        "age": (0.1, 15.),
        "tau": (0.3, 10.),
        "massformed": (1., 15.),
        "metallicity": (0., 2.5),
    }
    dust = {"type": "Calzetti", "Av": (0., 2.)}
    fit_instructions = {
        "redshift": (0.674, 0.676),
        "exponential": exp,
        "dust": dust,
    }

    def load_data(ID):
        return np.array([flam_sersic, flam_err]).T

    galaxy = pipes.galaxy("0206_sersic", load_data, spectrum_exists=False,
                          filt_list=filt_list, phot_units="ergscma")
    fit = pipes.fit(galaxy, fit_instructions, run="AGEL0206_sersic")

    print("\nRunning Bagpipes fit on Sersic-total photometry...")
    try:
        fit.fit(verbose=True, sampler="multinest", n_live=400)
    except (AttributeError, OSError) as e:
        print(f"MultiNest unavailable ({e}); using nautilus.")
        fit.fit(verbose=True, sampler="nautilus", n_live=400)

    samples = fit.posterior.samples
    log_M_sersic = samples["stellar_mass"]

    # Load nb02 aperture posterior
    ap_path = "results/bagpipes_sed_results.npz"
    if os.path.exists(ap_path):
        ap = np.load(ap_path, allow_pickle=True)
        log_M_aperture = ap["stellar_mass"]
        ap_p16, ap_p50, ap_p84 = np.percentile(log_M_aperture, [16, 50, 84])
    else:
        log_M_aperture = None
        ap_p50 = 11.33

    s_p16, s_p50, s_p84 = np.percentile(log_M_sersic, [16, 50, 84])

    print()
    print("=" * 70)
    print("STELLAR MASS COMPARISON")
    print("=" * 70)
    print(f"{'source':<25} {'p16':>7} {'p50':>7} {'p84':>7} {'-1σ':>6} {'+1σ':>6}")
    print("-" * 70)
    if log_M_aperture is not None:
        a16, a50, a84 = np.percentile(log_M_aperture, [16, 50, 84])
        print(f"{'aperture (nb02)':<25} {a16:>7.3f} {a50:>7.3f} {a84:>7.3f} "
              f"{a50-a16:>6.3f} {a84-a50:>6.3f}")
    print(f"{'Sersic-total (nb08)':<25} {s_p16:>7.3f} {s_p50:>7.3f} {s_p84:>7.3f} "
          f"{s_p50-s_p16:>6.3f} {s_p84-s_p50:>6.3f}")
    if log_M_aperture is not None:
        print(f"{'Δlog M* (Sersic−ap)':<25} {s_p50-a50:>+22.3f}")
        print(f"{'Δ in solar masses':<25} {10**(s_p50-a50):>22.2f}× more massive")
    print("=" * 70)

    # Save
    out = "results/bagpipes_sersic_refit.npz"
    save = dict(
        filter_names=filter_names,
        pivot_AA=pivot,
        mag_aperture=mag_aperture,
        mag_sersic=mag_sersic,
        flam_aperture=flam_aperture,
        flam_sersic=flam_sersic,
        flam_err_used=flam_err,
        log_M_sersic_samples=log_M_sersic,
        log_M_sersic_p16=s_p16,
        log_M_sersic_p50=s_p50,
        log_M_sersic_p84=s_p84,
    )
    if log_M_aperture is not None:
        save.update(dict(
            log_M_aperture_samples=log_M_aperture,
            log_M_aperture_p16=ap_p16, log_M_aperture_p50=ap_p50, log_M_aperture_p84=ap_p84,
            delta_log_M=s_p50 - ap_p50,
        ))
    np.savez(out, **save)
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    main()
