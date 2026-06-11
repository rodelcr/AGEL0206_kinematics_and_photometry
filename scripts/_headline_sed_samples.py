"""Regenerate the Fig-2 SED npz (with SAMPLE arrays, matching bagpipes_sed_results.npz
keys) for the HEADLINE aperture-corrected-total photometry. -> bagpipes_sed_results_aperture.npz"""
import os, sys, numpy as np
sys.path.insert(0, os.getcwd())
import bagpipes as pipes
from scripts.bagpipes_sersic_refit import mag_to_flam, FILT_LIST, FIT_INSTRUCTIONS
d = np.load("results/aperture_matched_photometry.npz", allow_pickle=True)
ORDER = list(d["filter_names"]); pivot = d["pivot"]; mags = d["mag_total_2"]
flam = mag_to_flam(mags, pivot); flam_err = 0.1 * flam
filt = [os.path.abspath(p) for p in FILT_LIST]
g = pipes.galaxy("AGEL0206_aper_total_2Re", lambda ID: np.array([flam, flam_err]).T,
                 spectrum_exists=False, filt_list=filt, phot_units="ergscma")
fit = pipes.fit(g, FIT_INSTRUCTIONS, run="AGEL0206_aper_total_2Re")
fit.fit(verbose=False)
fit.posterior.get_advanced_quantities()
s = fit.posterior.samples
np.savez("results/bagpipes_sed_results_aperture.npz",
         filter_names=np.array(ORDER), pivot_wavelengths_AA=pivot,
         fluxes_lambda=flam, fluxes_lambda_error=flam_err, fluxes_lambda_error_fit=flam_err,
         model_photometry=s["photometry"], model_spectrum=s["spectrum_full"],
         model_wavelengths=fit.posterior.model_galaxy.wavelengths,
         redshift=s["redshift"], stellar_mass=s["stellar_mass"])
p = np.percentile(s["stellar_mass"], [16, 50, 84])
print(f"Fig2 SED npz: logM*={p[1]:.3f}  spec{s['spectrum_full'].shape} phot{s['photometry'].shape}")
print("Saved → results/bagpipes_sed_results_aperture.npz")
