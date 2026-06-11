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


FILT_LIST = ["HST_WFC3_UVIS1.F200LP.dat", "HST_WFC3_IR.F140W.dat",
             "JWST_NIRCam.F150W2.dat", "JWST_NIRCam.F322W2.dat"]
FIT_INSTRUCTIONS = {
    "redshift": (0.674, 0.676),
    "exponential": {"age": (0.1, 15.), "tau": (0.3, 10.),
                    "massformed": (1., 15.), "metallicity": (0., 2.5)},
    "dust": {"type": "Calzetti", "Av": (0., 2.)},
}


def run_bagpipes_for_mags(mags_AB, pivot_AA, run_name, n_live=400, err_frac=0.1):
    """Fit Bagpipes to a 4-band AB-mag vector (order F200LP, F140W, F150W2, F322W2)
    with the nb02 priors and `err_frac` fractional flux errors (0.1 = 10%, default;
    pass 0.2 for the conservative 20%). Returns the stellar_mass posterior samples.
    Factored out of main() so callers (e.g. scripts/arc_mask_verification.py,
    scripts/photometry_systematics.py) can test alternate masks' photometry."""
    return fit_and_extract(mags_AB, pivot_AA, run_name,
                           n_live=n_live, err_frac=err_frac)["stellar_mass"]


def fit_and_extract(mags_AB, pivot_AA, run_name, n_live=400, err_frac=0.1):
    """Fit (or reload cached) Bagpipes and return a dict with the stellar-mass
    posterior AND the best-fit model SED + filter-convolved model photometry, so a
    caller can plot the measured points, the model spectrum, and the model photometry
    together. Keys:
        stellar_mass            : (n,) log10(M*/Msun) samples
        flam, flam_err          : (4,) measured F_lambda + error fed to the fit
        eff_wavs                : (4,) observed-frame filter effective wavelengths [AA]
        model_phot_p16/50/84    : (4,) filter-convolved MODEL photometry percentiles
        wav_obs                 : (nwav,) observed-frame model wavelength grid [AA]
        spec_p16/50/84          : (nwav,) model spectrum F_lambda percentiles
        z_p50                   : median redshift
    """
    flam = mag_to_flam(mags_AB, pivot_AA)
    flam_err = err_frac * flam
    filt_list = [os.path.abspath(p) for p in FILT_LIST]

    def load_data(ID):
        return np.array([flam, flam_err]).T

    galaxy = pipes.galaxy(run_name, load_data, spectrum_exists=False,
                          filt_list=filt_list, phot_units="ergscma")
    fit = pipes.fit(galaxy, FIT_INSTRUCTIONS, run=run_name)
    try:
        fit.fit(verbose=False, sampler="multinest", n_live=n_live)
    except (AttributeError, OSError) as e:
        print(f"MultiNest unavailable ({e}); using nautilus.")
        fit.fit(verbose=False, sampler="nautilus", n_live=n_live)

    fit.posterior.get_advanced_quantities()
    s = fit.posterior.samples
    z_p50 = float(np.median(s["redshift"]))
    wav_rest = fit.posterior.model_galaxy.wavelengths        # rest-frame [AA]
    wav_obs = wav_rest * (1.0 + z_p50)                       # observed-frame [AA]
    spec_p16, spec_p50, spec_p84 = np.percentile(s["spectrum_full"], [16, 50, 84], axis=0)
    mp16, mp50, mp84 = np.percentile(s["photometry"], [16, 50, 84], axis=0)
    return dict(
        stellar_mass=s["stellar_mass"],
        flam=flam, flam_err=flam_err,
        eff_wavs=np.asarray(fit.galaxy.filter_set.eff_wavs, dtype=float),
        model_phot_p16=mp16, model_phot_p50=mp50, model_phot_p84=mp84,
        wav_obs=wav_obs, spec_p16=spec_p16, spec_p50=spec_p50, spec_p84=spec_p84,
        z_p50=z_p50,
    )


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

    # Fit BOTH photometry choices, extracting the model SED + filter-convolved model
    # photometry (so the figure can show measured points, model spectrum, model phot).
    # AGEL0206_sersic is cached (reloads instantly); the aperture fig-run fits fresh.
    print("\nBagpipes fit on Sersic-total photometry...")
    ser = fit_and_extract(mag_sersic, pivot, "AGEL0206_sersic")
    print("Bagpipes fit on empirical aperture photometry (masked, not filled)...")
    apr = fit_and_extract(mag_aperture, pivot, "AGEL0206_aperture_fig")

    log_M_sersic = ser["stellar_mass"]
    log_M_aperture = apr["stellar_mass"]
    ap_p16, ap_p50, ap_p84 = np.percentile(log_M_aperture, [16, 50, 84])
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

    # Save: M* posteriors + measured/model photometry + model SED for BOTH choices
    out = "results/bagpipes_sersic_refit.npz"
    save = dict(
        filter_names=filter_names,
        pivot_AA=pivot,
        mag_aperture=mag_aperture,
        mag_sersic=mag_sersic,
        flam_aperture=flam_aperture,
        flam_sersic=flam_sersic,
        flam_err_used=apr["flam_err"],
        log_M_sersic_samples=log_M_sersic,
        log_M_sersic_p16=s_p16, log_M_sersic_p50=s_p50, log_M_sersic_p84=s_p84,
        log_M_aperture_samples=log_M_aperture,
        log_M_aperture_p16=ap_p16, log_M_aperture_p50=ap_p50, log_M_aperture_p84=ap_p84,
        delta_log_M=s_p50 - ap_p50,
        # model SED + filter-convolved model photometry (for the figure)
        eff_wavs=apr["eff_wavs"],
        ap_flam=apr["flam"], ap_flam_err=apr["flam_err"],
        ap_model_phot_p16=apr["model_phot_p16"], ap_model_phot_p50=apr["model_phot_p50"],
        ap_model_phot_p84=apr["model_phot_p84"],
        ap_wav_obs=apr["wav_obs"], ap_spec_p16=apr["spec_p16"],
        ap_spec_p50=apr["spec_p50"], ap_spec_p84=apr["spec_p84"], ap_z_p50=apr["z_p50"],
        ser_flam=ser["flam"], ser_flam_err=ser["flam_err"],
        ser_model_phot_p16=ser["model_phot_p16"], ser_model_phot_p50=ser["model_phot_p50"],
        ser_model_phot_p84=ser["model_phot_p84"],
        ser_wav_obs=ser["wav_obs"], ser_spec_p16=ser["spec_p16"],
        ser_spec_p50=ser["spec_p50"], ser_spec_p84=ser["spec_p84"], ser_z_p50=ser["z_p50"],
    )
    np.savez(out, **save)
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    main()
