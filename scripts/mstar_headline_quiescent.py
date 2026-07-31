"""Spectrum-consistent (quiescence-constrained) headline M* + stellar properties (2026-06-15).

The 4-band SED is degenerate (age-dust-M/L): the fiducial flat-age prior slides onto a
young+dusty branch (SFR~57) that the KCWI absorption-line spectrum (strong Ca H&K/G/Mg b/Fe,
NO emission) rules out (degeneracy check: ΔlnZ=-0.18, identical χ², scripts/sed_quiescence_check.py).
Per user decision (2026-06-15), the HEADLINE adopts a passive SFH prior (old age, short tau,
low dust) motivated by the spectrum; the fiducial-vs-quiescent spread is carried as the
SFH-prior systematic.

Runs the validated-photometry estimator ladder (raw/raw_apcorr/filled/total) + the Sérsic-total
+ the auto-fit total all under the QUIESCENT prior, plus the fiducial validated-total for the
SFH-prior systematic. Outputs results/mstar_headline_quiescent.npz.

Run:  python -m scripts.mstar_headline_quiescent
"""
import numpy as np, os
import bagpipes as pipes
from scripts.bagpipes_sersic_refit import FILT_LIST, mag_to_flam

FILT = [os.path.abspath(p) for p in FILT_LIST]

QUI = {"redshift": (0.674, 0.676),
       "exponential": {"age": (4., 15.), "tau": (0.1, 1.5), "massformed": (1., 15.), "metallicity": (0., 2.5)},
       "dust": {"type": "Calzetti", "Av": (0., 0.6)}}
FID = {"redshift": (0.674, 0.676),
       "exponential": {"age": (0.1, 15.), "tau": (0.3, 10.), "massformed": (1., 15.), "metallicity": (0., 2.5)},
       "dust": {"type": "Calzetti", "Av": (0., 2.)}}

amv = np.load("results/aperture_correction_validated.npz", allow_pickle=True)
am = np.load("results/aperture_matched_photometry.npz", allow_pickle=True)
tab = np.load("results/sersic_parameter_table.npz", allow_pickle=True)
pivot = amv["pivot"]
ORDER = list(amv["filter_names"])
mtot_sersic = np.array([float({k: v for k, v in tab[f"{b}_cen"]}["m_tot"]) for b in ORDER])

def fit_mags(mags, instr, run):
    flam = mag_to_flam(mags, pivot); flam_err = 0.10 * flam
    g = pipes.galaxy(run, lambda ID: np.array([flam, flam_err]).T,
                     spectrum_exists=False, filt_list=FILT, phot_units="ergscma")
    f = pipes.fit(g, instr, run=run)
    try:
        f.fit(verbose=False, sampler="multinest", n_live=400)
    except (AttributeError, OSError):
        f.fit(verbose=False, sampler="nautilus", n_live=400)
    f.posterior.get_advanced_quantities()
    return f.posterior.samples, f

def med3(s, k):
    a = np.asarray(s[k], float); p = np.percentile(a, [16, 50, 84]); return p[1], p[1]-p[0], p[2]-p[1]

out = {"filter_names": np.array(ORDER), "pivot": pivot}

# --- estimator ladder under the QUIESCENT (headline) prior ---
print("QUIESCENT-prior fits (headline):")
for kind, mags in [("raw", amv["mag_raw_2"]), ("raw_apcorr", amv["mag_raw_apcorr_2"]),
                   ("filled", amv["mag_filled_2"]), ("total", amv["mag_total_2"]),
                   ("sersic", mtot_sersic)]:
    s, f = fit_mags(mags, QUI, f"AGEL0206_QUI_{kind}")
    p = np.percentile(np.asarray(s["stellar_mass"], float), [16, 50, 84])
    out[f"logM_{kind}_qui"] = p
    if kind == "total":
        out["samples_logM_total_qui"] = np.asarray(s["stellar_mass"], float)
        for nm, key in [("mass_weighted_age", "mass_weighted_age"), ("formed_mass", "formed_mass"),
                        ("sfr", "sfr"), ("ssfr", "ssfr"), ("metallicity", "exponential:metallicity"),
                        ("dust_Av", "dust:Av"), ("tau", "exponential:tau"), ("age_form", "exponential:age")]:
            if key in s:
                m, lo, hi = med3(s, key); out[f"prop_{nm}"] = np.array([m, lo, hi])
        # SED model spectrum + filter-convolved photometry for Fig 2 (quiescent, headline-consistent)
        zp50 = float(np.median(s["redshift"]))
        out["sed_z_p50"] = zp50
        out["sed_wav_obs"] = f.posterior.model_galaxy.wavelengths * (1.0 + zp50)
        sp16, sp50, sp84 = np.percentile(s["spectrum_full"], [16, 50, 84], axis=0)
        out["sed_spec_p16"], out["sed_spec_p50"], out["sed_spec_p84"] = sp16, sp50, sp84
        mp16, mp50, mp84 = np.percentile(s["photometry"], [16, 50, 84], axis=0)
        out["sed_model_phot_p16"], out["sed_model_phot_p50"], out["sed_model_phot_p84"] = mp16, mp50, mp84
        _fl = mag_to_flam(mags, pivot)
        out["sed_flam"], out["sed_flam_err"] = _fl, 0.10 * _fl
    print(f"  {kind:11s}: logM* = {p[1]:.3f} (-{p[1]-p[0]:.3f}/+{p[2]-p[1]:.3f})")

# --- auto-fit total under QUIESCENT prior (for apcorr-model systematic, prior-consistent) ---
s_auto, _ = fit_mags(am["mag_total_2"], QUI, "AGEL0206_QUI_autototal")
out["logM_autototal_qui"] = np.percentile(np.asarray(s_auto["stellar_mass"], float), [16, 50, 84])

# --- fiducial validated total (for SFH-prior systematic) ---
s_fid, _ = fit_mags(amv["mag_total_2"], FID, "AGEL0206_FID_total")
out["logM_total_fid"] = np.percentile(np.asarray(s_fid["stellar_mass"], float), [16, 50, 84])

# --- systematics (all on the validated total, prior-consistent where possible) ---
lm_qui = float(out["logM_total_qui"][1])
out["sys_apcorr_model_dex"] = abs(float(out["logM_autototal_qui"][1]) - lm_qui)   # auto vs validated (quiescent)
out["sys_sfh_prior_dex"] = abs(float(out["logM_total_fid"][1]) - lm_qui)          # fiducial vs quiescent
print(f"\nHEADLINE (quiescent) logM* (H0=70) = {lm_qui:.3f}")
print(f"apcorr-model sys = {out['sys_apcorr_model_dex']:.3f} dex   SFH-prior sys = {out['sys_sfh_prior_dex']:.3f} dex")
print("stellar props (quiescent):",
      {k.replace('prop_',''): round(float(out[k][0]),3) for k in out if k.startswith('prop_')})

np.savez("results/mstar_headline_quiescent.npz", **out)
print("Saved -> results/mstar_headline_quiescent.npz")
