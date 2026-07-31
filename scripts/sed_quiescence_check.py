"""Double-check the SED star-forming result against the spectrum (audit 2026-06-15).

The KCWI spectrum shows strong Ca H&K / G-band / Mg b / Fe absorption and NO emission
lines -> an OLD, passive population. The broadband SED fit (4 points) instead returns
SFR~57, sSFR~-9.6, age~2.6 Gyr. This script tests whether that is just the young+dusty
branch of the age-dust-M/L degeneracy by re-fitting the SAME validated total photometry
with an OLD-restricted SFH prior, and comparing fit quality (lnZ, chi^2 to the 4 points),
M*, age, and SFR. If the old solution fits comparably, the photometry is degenerate and
the spectrum (old) is decisive -> the SF result is an artifact.

Run:  python -m scripts.sed_quiescence_check
"""
import numpy as np, os
import bagpipes as pipes
from scripts.bagpipes_sersic_refit import FILT_LIST, mag_to_flam

v = np.load("results/aperture_correction_validated.npz", allow_pickle=True)
mags = v["mag_total_2"]; pivot = v["pivot"]
flam = mag_to_flam(mags, pivot); flam_err = 0.10 * flam
filt = [os.path.abspath(p) for p in FILT_LIST]

def load_data(ID):
    return np.array([flam, flam_err]).T

FITS = {
    "fiducial_exp": {              # current headline prior
        "redshift": (0.674, 0.676),
        "exponential": {"age": (0.1, 15.), "tau": (0.3, 10.), "massformed": (1., 15.), "metallicity": (0., 2.5)},
        "dust": {"type": "Calzetti", "Av": (0., 2.)},
    },
    "old_quiescent": {             # OLD-restricted: formed early, declined; little dust
        "redshift": (0.674, 0.676),
        "exponential": {"age": (4., 15.), "tau": (0.1, 1.5), "massformed": (1., 15.), "metallicity": (0., 2.5)},
        "dust": {"type": "Calzetti", "Av": (0., 0.6)},
    },
}

def chi2(model_phot):
    return float(np.sum(((flam - model_phot) / flam_err) ** 2))

print(f"{'fit':16s}{'logM*':>8s}{'age_mw':>8s}{'SFR':>9s}{'sSFR':>8s}{'Av':>6s}{'lnZ':>9s}{'chi2/4':>8s}")
print("-" * 72)
res = {}
for name, instr in FITS.items():
    g = pipes.galaxy(name + "_QCHK", load_data, spectrum_exists=False, filt_list=filt, phot_units="ergscma")
    f = pipes.fit(g, instr, run=name + "_QCHK")
    try:
        f.fit(verbose=False, sampler="multinest", n_live=400)
    except (AttributeError, OSError):
        f.fit(verbose=False, sampler="nautilus", n_live=400)
    f.posterior.get_advanced_quantities()
    s = f.posterior.samples
    lm = np.median(s["stellar_mass"]); age = np.median(s["mass_weighted_age"])
    sfr = np.median(s["sfr"]); ssfr = np.median(s["ssfr"]); av = np.median(s["dust:Av"])
    lnz = float(f.results["lnz"])
    c2 = np.median([chi2(mp) for mp in s["photometry"]])
    res[name] = dict(lm=lm, age=age, sfr=sfr, ssfr=ssfr, av=av, lnz=lnz, c2=c2)
    print(f"{name:16s}{lm:8.3f}{age:8.2f}{sfr:9.1f}{ssfr:8.2f}{av:6.2f}{lnz:9.2f}{c2:8.2f}")

dlnz = res["old_quiescent"]["lnz"] - res["fiducial_exp"]["lnz"]
dlm = res["old_quiescent"]["lm"] - res["fiducial_exp"]["lm"]
print(f"\nΔlnZ (old - fiducial) = {dlnz:+.2f}   (|Δ|<~2.5 => photometry cannot distinguish them)")
print(f"Δlog M* (old - fiducial) = {dlm:+.3f} dex  (old population => higher M/L => higher M*)")
print("If ΔlnZ is small, the 4-band SED is degenerate; the absorption-line spectrum (old) is decisive.")
