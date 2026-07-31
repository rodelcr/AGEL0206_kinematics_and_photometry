"""Extract the full set of Bagpipes stellar-population properties from the FINAL fit
(validated-fit total-light run, AGEL0206_aperVAL_total_2Re) — the SAME fit that sets the
headline log M*. Reloads the cached posterior (no re-fit), reports + saves all derived
quantities. Run from repo root:  python -m scripts.sed_properties_final
"""
import numpy as np
import bagpipes as pipes
import os
from scripts.bagpipes_sersic_refit import FIT_INSTRUCTIONS, FILT_LIST, mag_to_flam

RUN = "AGEL0206_aperVAL_total_2Re"
v = np.load("results/aperture_correction_validated.npz", allow_pickle=True)
mags = v["mag_total_2"]; pivot = v["pivot"]
flam = mag_to_flam(mags, pivot); flam_err = 0.10 * flam

def load_data(ID):
    return np.array([flam, flam_err]).T

galaxy = pipes.galaxy(RUN, load_data, spectrum_exists=False,
                      filt_list=[os.path.abspath(p) for p in FILT_LIST], phot_units="ergscma")
fit = pipes.fit(galaxy, FIT_INSTRUCTIONS, run=RUN)
fit.fit(verbose=False, sampler="multinest", n_live=400)   # cached -> reloads instantly
fit.posterior.get_advanced_quantities()
s = fit.posterior.samples

# friendly name -> (sample key, unit, log?)
WANT = [
    ("stellar_mass",     "stellar_mass",            "log10(M/Msun)"),
    ("formed_mass",      "formed_mass",             "log10(M/Msun)"),
    ("mass_weighted_age","mass_weighted_age",       "Gyr"),
    ("sfr",              "sfr",                     "Msun/yr"),
    ("ssfr",             "ssfr",                    "log10(1/yr)"),
    ("metallicity",      "exponential:metallicity", "Zsun"),
    ("dust_Av",          "dust:Av",                 "mag"),
    ("tau",              "exponential:tau",         "Gyr"),
    ("age_form",         "exponential:age",         "Gyr (time since SF onset)"),
    ("redshift",         "redshift",                ""),
]
out = {}
print(f"{'property':18s}{'p50':>10s}{'-1sig':>9s}{'+1sig':>9s}   unit")
print("-" * 60)
for name, key, unit in WANT:
    if key not in s:
        print(f"{name:18s}  (key '{key}' not in samples)"); continue
    arr = np.asarray(s[key], float)
    p16, p50, p84 = np.percentile(arr, [16, 50, 84])
    out[name] = (p50, p50 - p16, p84 - p50)
    print(f"{name:18s}{p50:10.3f}{p50-p16:9.3f}{p84-p50:9.3f}   {unit}")

np.savez("results/sed_properties_final.npz",
         names=np.array(list(out.keys())),
         p50=np.array([out[k][0] for k in out]),
         elo=np.array([out[k][1] for k in out]),
         ehi=np.array([out[k][2] for k in out]),
         run=RUN)
print("\nSaved -> results/sed_properties_final.npz")
