"""Sensitivity test: M* with the HST F200LP band DROPPED (2026-06-17).

Question (user, 2026-06-17): what stellar mass do we recover from JUST F140W +
the two JWST bands (F150W2, F322W2), i.e. without the bluest band F200LP?

F200LP (pivot ~4923 Å) is the only band blueward of the rest-frame 4000 Å break
at z=0.676 (rest ~2940 Å — samples the near-UV / young-star + dust-sensitive
side). Dropping it removes the lever arm on recent SF / dust, so this is also a
check on how much the blue point drives the M*/L (and whether the young+dusty
degeneracy the spectrum rules out was being pinned by F200LP).

Refits the VALIDATED-TOTAL photometry (results/aperture_correction_validated.npz,
mag_total_2) with a 3-band filter set, under BOTH the headline QUIESCENT prior and
the FIDUCIAL flat-age prior, and applies the same +0.0282 dex Planck-2015 cosmology
shift used by scripts/paper_values.py so the numbers are directly comparable to the
4-band headline log(M*/Msun) = 11.46 (registry, Planck 2015).

Run from repo root:  python -m scripts.mstar_drop_f200lp
"""
import os
import numpy as np
import bagpipes as pipes
from scripts.bagpipes_sersic_refit import FILT_LIST, mag_to_flam

DLOGM_COSMO = 0.0282                      # H0=70 -> Planck 2015 D_L^2 rescale (paper_values.py)
ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]

QUI = {"redshift": (0.674, 0.676),
       "exponential": {"age": (4., 15.), "tau": (0.1, 1.5), "massformed": (1., 15.), "metallicity": (0., 2.5)},
       "dust": {"type": "Calzetti", "Av": (0., 0.6)}}
FID = {"redshift": (0.674, 0.676),
       "exponential": {"age": (0.1, 15.), "tau": (0.3, 10.), "massformed": (1., 15.), "metallicity": (0., 2.5)},
       "dust": {"type": "Calzetti", "Av": (0., 2.)}}


def fit_mags(mags, pivot, filt_list, instr, run):
    flam = mag_to_flam(mags, pivot); flam_err = 0.10 * flam
    g = pipes.galaxy(run, lambda ID: np.array([flam, flam_err]).T,
                     spectrum_exists=False, filt_list=filt_list, phot_units="ergscma")
    f = pipes.fit(g, instr, run=run)
    try:
        f.fit(verbose=False, sampler="multinest", n_live=400)
    except (AttributeError, OSError):
        f.fit(verbose=False, sampler="nautilus", n_live=400)
    f.posterior.get_advanced_quantities()
    return f.posterior.samples


def med3(s, k):
    a = np.asarray(s[k], float); p = np.percentile(a, [16, 50, 84]); return p[1], p[1] - p[0], p[2] - p[1]


def main():
    amv = np.load("results/aperture_correction_validated.npz", allow_pickle=True)
    mag_total = amv["mag_total_2"]
    pivot_all = amv["pivot"]
    filt_all = [os.path.abspath(p) for p in FILT_LIST]

    keep = [i for i, b in enumerate(ORDER) if b != "F200LP"]   # F140W, F150W2, F322W2
    mags3 = mag_total[keep]; pivot3 = pivot_all[keep]; filt3 = [filt_all[i] for i in keep]
    print("Dropping F200LP. Bands used:", [ORDER[i] for i in keep])
    for i in keep:
        print(f"  {ORDER[i]:7s} pivot={pivot_all[i]:8.0f} AA  m_total={mag_total[i]:.3f}")

    out = {"bands_used": np.array([ORDER[i] for i in keep]), "pivot": pivot3,
           "mag_total_3band": mags3, "dlogm_cosmo": DLOGM_COSMO}

    print("\n--- 3-band (no F200LP) fits on validated-total photometry ---")
    for tag, prior in [("qui", QUI), ("fid", FID)]:
        s = fit_mags(mags3, pivot3, filt3, prior, f"AGEL0206_noF200LP_{tag}")
        p = np.percentile(np.asarray(s["stellar_mass"], float), [16, 50, 84])   # H0=70
        pP = p + DLOGM_COSMO                                                     # Planck 2015
        out[f"logM_3band_{tag}_H070"] = p
        out[f"logM_3band_{tag}_planck"] = pP
        age = med3(s, "mass_weighted_age") if "mass_weighted_age" in s else (np.nan, np.nan, np.nan)
        sfr = med3(s, "sfr") if "sfr" in s else (np.nan, np.nan, np.nan)
        av = med3(s, "dust:Av") if "dust:Av" in s else (np.nan, np.nan, np.nan)
        out[f"age_3band_{tag}"] = np.array(age); out[f"sfr_3band_{tag}"] = np.array(sfr)
        out[f"Av_3band_{tag}"] = np.array(av)
        label = "QUIESCENT (headline prior)" if tag == "qui" else "FIDUCIAL (flat-age prior)"
        print(f"\n{label}:")
        print(f"  logM* (H0=70)    = {p[1]:.3f}  (-{p[1]-p[0]:.3f}/+{p[2]-p[1]:.3f})")
        print(f"  logM* (Planck15) = {pP[1]:.3f}  (-{pP[1]-pP[0]:.3f}/+{pP[2]-pP[1]:.3f})   [4-band headline = 11.46]")
        print(f"  mass-wtd age = {age[0]:.2f} Gyr   SFR = {sfr[0]:.2f} Msun/yr   Av = {av[0]:.2f}")

    # 4-band headline reference: quiescent-prior validated-total (= registry M* source)
    hq = np.load("results/mstar_headline_quiescent.npz", allow_pickle=True)
    p4 = hq["logM_total_qui"] + DLOGM_COSMO
    out["logM_4band_qui_planck"] = p4
    print("\n(reference) 4-band QUIESCENT validated-total, Planck = %.3f (registry headline 11.46)" % p4[1])
    np.savez("results/mstar_drop_f200lp.npz", **out)
    print("\nSaved -> results/mstar_drop_f200lp.npz")


if __name__ == "__main__":
    main()
