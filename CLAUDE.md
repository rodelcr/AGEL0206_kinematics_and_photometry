# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis of the deflector galaxy in the AGEL J020613-011417 strong gravitational lens system. Two main analyses:
1. **Stellar velocity dispersion** from Keck/KCWI IFU spectroscopy using ppxf
2. **Aperture photometry + SED fitting** from HST and JWST imaging using Bagpipes

## Key Conventions

- All wavelengths in Angstroms unless otherwise noted
- Redshifts: deflector systemic `z = 0.67564` (from notebook 04 line-fitting; older notebooks use `0.67511` from AGEL DR2); source `z = 1.302`
- Coordinate system: ICRS, target at RA=31.55611, Dec=-1.23817
- IFU data cube axes: `[wavelength, y_spatial, x_spatial]` (shape ~[1273, 90, 90])
- Deflector spaxel region: `cube[:,45:64,45:55]` (19x10 spaxels) — used for the *integrated* spectrum in notebooks 01/03/03b/03c. This region is contaminated by the lensed arc; for the final σ measurements use circular apertures at R < Re/8, Re/2, Re instead (notebook 06).
- Noise/sky region: `cube[:,28:40,45:70]` (12x25 spaxels)
- Photometry uses AB magnitude system throughout
- Seeing FWHM = 1.27" (from KCWI `GUIDFWHM` header); IFU spaxel scale = 0.30"/spaxel; 1" ≈ 7.04 kpc at z=0.67564

## Notebook Organization

There are two versions of the base analysis notebooks (01, 02), plus a sequence of follow-up notebooks (03–06) that build on them:

- **Original exploratory notebooks** (`01_IFU_spectra_extraction_and_ppxf.ipynb`, `02_Bagpipes_SED_fitting.ipynb`): Full history of the analysis with iterative attempts. Kept for reference.
- **Streamlined notebooks** (`01_streamlined_IFU_ppxf.ipynb`, `02_streamlined_Bagpipes_SED.ipynb`): Clean versions with redundant cells removed. Use these for running the analysis.
- `03_bootstrap_ppxf_errors.ipynb` — wild-bootstrap error estimation on the 190-spaxel integrated spectrum at z=0.67511.
- `03b_bootstrap_ppxf_z067564.ipynb` — same bootstrap at z=0.67564 (the systemic redshift).
- `03c_bootstrap_summary.ipynb` — combines both input redshifts and all three SPS libraries (FSPS/EMILES/XSL); identifies the stable polynomial-degree range; reports the σ on the 190-spaxel integrated region (~301 ± 30 km/s — **NOT** the number to report in the paper).
- `04_redshift_verification.ipynb` — independent line-by-line Gaussian fitting → systemic z = 0.67564.
- `05_radial_sigma_profile.ipynb` — annular σ(R) profile with contamination-weighted masking.
- `05x_sigma_discrepancy_diagnostic.ipynb` — Tests 1–19 resolving the σ-vs-aperture discrepancy between NB01 (~301 km/s, 190-spaxel box) and NB05 (~210 km/s, inner circular aperture). **Resolution: σ depends on aperture; the 190-spaxel box is contaminated by the lensed arc and by a genuine outer rise.** Includes a failed Voronoi-binning attempt (Test 18) and a PowerBin attempt (Test 19) wired to `scripts/run_powerbin_test19.py`.
- `06_final_sigma_Re_apertures.ipynb` — **final paper σ values** at R < Re/8 (JFK95 σ_c), R < Re/2 (gradient diagnostic), and R < Re (σ_e, the Kormendy & Ho 2013 / Greene+2020 ARA&A definition). Applies 03c's FSPS+EMILES+XSL × multi-degree posterior combination to each circular aperture.

The streamlined IFU notebook uses veldis (degree=[4,30]) for the integrated spectrum and raw ppxf for per-spaxel and power-binned fitting. The `ppxf_per_spaxel()` function appends `best_fit` twice per degree iteration — account for this when indexing results (`best_fit_idx = deg_idx * 2`). This bug is **only** in notebook 01_streamlined; `scripts/bootstrap_ppxf.py` (which notebooks 03/03b/03c/06 build on) does not have it.

## Key Results

- **σ(<Re/8)**: *populated by notebook 06*
- **σ(<Re/2)**: *populated by notebook 06*
- **σ(<Re) = σ_e**: *populated by notebook 06* — **this is the primary number** for the M•–σ relation (Kormendy & Ho 2013 eq. 3, Greene+2020 fig. 5).
- σ on the 190-spaxel integrated region `cube[:,45:64,45:55]`: ~301 ± 30 km/s — **do NOT cite this in the paper**; it conflates the inner stellar kinematics with the lensed-arc contamination.
- Systemic redshift: z = 0.67564 (notebook 04); older results use z = 0.67511.
- Effective radius (IFU white-light, proper masking): Re = 2.61" = 18.4 kpc. Re/2 ≈ 1.30", Re/8 ≈ 0.33".
- log(M★/M☉) = 11.33 +0.07/−0.09 (notebook 02 Bagpipes).

## Scripts

All scripts live in `scripts/` and are importable as modules (`from scripts.bootstrap_ppxf import run_bootstrap`). Most also have a `__main__` CLI.

- `scripts/bootstrap_ppxf.py` — wild-bootstrap σ/V errors for ppxf fits. Hybrid Rademacher sign-flip × local-residual scaling (rolling-window 75 pix). Key functions:
  - `setup_ppxf_inputs(ifu_file, sps_name, z, ...)` — builds the ppxf_inputs dict from the cube+fixed 190-spaxel slice.
  - `setup_ppxf_inputs_from_spectrum(flux, noise, hdr, sps_name, z, ...)` — same dict but from a pre-extracted 1-D spectrum on the cube wavelength grid (used by notebook 06 for circular-aperture spectra).
  - `run_bootstrap_single_degree(ppxf_inputs, degree, best_fit_spectrum, n_bootstrap, window, seed)` — inner bootstrap loop.
  - `run_bootstrap(ifu_file, sps_name, ...)` — full pipeline for one SPS library; writes `results/ppxf_bootstrap_errors_{sps}[_<save_suffix>].npz`.
  - CLI: `python scripts/bootstrap_ppxf.py --sps_name {fsps,emiles,xsl,all} --n_bootstrap 500 [--save_suffix z067564]`. ~12 min per SPS library at N=500.
- `scripts/redshift_verify.py` — importable module. Provides `ABSORPTION_LINES` and `EMISSION_LINES` dictionaries (air wavelengths, NIST-verified), Gaussian single-line and multi-line fitters, and Ca H+K and [OII] doublet fitters with shared parameters. Used by notebook 04.
- `scripts/measure_Re.py` — multi-source effective-radius measurement (IFU white-light, HST F140W, HST F200LP) with four masking strategies (unmasked, zeroed, proper, PSF-convolved). CLI: `python scripts/measure_Re.py [--plot]`. Writes `results/Re_measurements.npz` with 10 labeled values; the paper uses the `IFU_whitelight_proper` row (Re = 2.61").
- `scripts/run_powerbin_test19.py` — standalone PowerBin + ppxf runner referenced from notebook 05x Test 19. Writes `results/test19_powerbin_results.npz` and PNG diagnostics.
- `scripts/photometry_masking_HST.py`, `scripts/photometry_masking_JWST.py` — interactive aperture photometry tools (see "Photometry Scripts" below).

## Root-level Python modules

`__init__.py`, `coords.py`, and `datafuncs.py` at the repo root are a **vendored copy of C. Fassnacht's `cdfutils`** (also mirrored at `notebooks/cdfutils/`). Treat as read-only third-party code; do not modify. They provide WCS helpers, sigma-clipping, and generic data utilities used by the photometry scripts.

## Results directory

`results/` contains the `.npz` outputs that notebooks and scripts consume downstream. Never regenerate these from scratch unless you have to — the bootstrap runs are expensive (~12 min per SPS library at N=500).

Key schemas:
- `ppxf_integrated_spectrum_results_{fsps,emiles,xsl}.npz` — best-fit arrays on the 190-spaxel integrated spectrum, 30 degrees each. Keys: `galaxy`, `best_fit` (30, n_pix), `degrees`, `vel_dis`, `mean_vel`, `fit_chi2`, `goodpixels`, `lam_gal_rest`, `noise_spec`.
- `ppxf_bootstrap_errors_{sps}[_z067564].npz` — full wild-bootstrap arrays shape `(30, 500)` for V, σ, χ², z, plus summary percentiles `{sigma,V,z}_p16/50/84`, asymmetric errors `{sigma,V}_boot_err_lo/hi`, originals `{sigma,V,chi2}_original`, metadata `degrees`, `sps_name`, `n_bootstrap`, `window`, `seed`, optional `z_input`.
- `ppxf_bootstrap_errors_R_le_{Re8,Re2,Re}_{sps}[_z067564].npz` — same schema as above, produced by notebook 06 for the circular apertures.
- `final_sigma_Re_apertures.npz` — concatenated posterior samples at the three apertures (produced by notebook 06). Keys: `sigma_samples_{Re8,Re2,Re}`, `sigma_{Re8,Re2,Re}_p16/50/84`, `aperture_arcsec`, `aperture_kpc`, `n_spaxels_{Re8,Re2,Re}`, `degree_range_used`, `sps_libs_used`, `z_inputs_used`, `n_bootstrap`. **This is the file Figure 2 of the ApJL loads.**
- `Re_measurements.npz`, `radial_sigma_profile.npz`, `radial_sigma_combined_posterior.npz`, `bagpipes_sed_results.npz`, `test19_powerbin_results.npz`.

## Running the Code

- Use the **ISMGas** conda environment: `conda activate ISMGas`
- Notebooks are in `notebooks/` and assume data files are in the repository root
- The photometry scripts are standalone interactive tools: `python scripts/photometry_masking_HST.py`
- Bagpipes caches results in `pipes/posterior/` as HDF5 files; delete the `.h5` file to re-run a fit
- For local development, symlink data files from the original directory (see Data Files section)

## Dependencies

See `requirements.txt`. Critical packages: ppxf, veldis, bagpipes, photutils, powerbin.
Bagpipes requires either MultiNest (C library) or Nautilus (pure Python) as a sampling backend.

## Data Files (Not in Repo)

Large FITS files must be obtained separately and placed in the repo root:
- **IFU cube:** `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` (KCWI, ~253 MB)
- **HST images:** `AGEL020613-011417A_F200LP_WFC3_*.fits`, `AGEL020613-011417A_F140W_WFC3_*.fits`
- **JWST images:** `jw05594-o101_t103_nircam_clear-f150w2_i2d.fits`, `jw05594-o101_t103_nircam_clear-f322w2_i2d.fits`
- **Stellar templates:** `TEXT/` directory (1,273 spectra for veldis), `spectra_*.npz` files (for ppxf)

## Photometry Scripts

`scripts/photometry_masking_HST.py` and `scripts/photometry_masking_JWST.py` share the same architecture:
- `run_photometry_math()` — core function for aperture stats, error propagation, magnitude conversion
- `PhotometryTool` class — interactive matplotlib GUI with sliders and click-to-mask
- Main script handles 3 execution modes: reload local params, import from another file, or interactive
- HST version uses `PHOTFLAM`/`PHOTPLAM` header keywords for AB zero point
- JWST version uses `PIXAR_SR` keyword: `AB_ZP = -6.10 - 2.5*log10(PIXAR_SR)`

## Photometric Calibration Reference

### AB Zeropoints from FITS Headers
For drizzled HST images (`_drc.fits`, `_drz.fits`) in units of e-/s:
```python
ZP_ST = -2.5 * log10(PHOTFLAM) - 21.10
ZP_AB = ZP_ST - 5.0 * log10(PHOTPLAM) + 18.6921
mag_AB = -2.5 * log10(count_rate_e_per_s) + ZP_AB
```
**Do NOT hardcode zeropoints** — read `PHOTFLAM` and `PHOTPLAM` from each image header.

### Known Zeropoints (for quick reference, verify against headers)
| Instrument/Filter | ZP_AB | PHOTFLAM | PHOTPLAM |
|-------------------|-------|----------|----------|
| ACS/WFC F606W | 26.483 | 7.9197e-20 | 5921.89 |
| WFC3/UVIS F200LP | 27.344 | 5.1851e-20 | 4923.48 |
| WFC3/IR F140W | ~26.45 | 1.4737e-20 | 13922.9 |

### Surface Brightness
```
SB (mag/arcsec^2) = -2.5 * log10(flux_per_pix / pix_area_arcsec2) + ZP_AB
```
where `pix_area_arcsec2 = pixel_scale^2`.

### Pixel Scales
| Detector | Pixel scale |
|----------|-------------|
| ACS/WFC | 0.05"/pix |
| WFC3/UVIS | 0.04"/pix |
| WFC3/IR | 0.13"/pix |

### Drizzled Image Properties
- `_drc.fits` / `_drz.fits` are in **e-/s** (count rate), already combined
- `NDRIZIM` = number of input images (exposures × chips)
- `EXPTIME` = total effective exposure time
- No pixel area map (PAM) correction needed for drizzled images
- **Drizzle pixel correlation**: noise in drizzled images is correlated between adjacent pixels. Empirical noise is ~1.3-2x higher than the idealized CCD equation predicts. The STScI ETC does NOT account for this.

## Exposure Time Calculator

A local ETC script is at `2026_HST_Proposal/Figure_generation/hst_etc.py` — note this lives in a **sibling directory**, not inside this repo.
Implements the WFC3 IHB CCD equation:
```
SNR = C*t / sqrt(C*t + (Bsky + Bdark)*t + Nread*RN^2)
```
(per pixel, for extended sources where C = count rate per pixel)

Validated against WFC3 IHB examples 9.9.1 and 9.9.3 (within 2%).

### Key Detector Parameters for ETC
| Parameter | WFC3/UVIS | WFC3/IR | ACS/WFC |
|-----------|-----------|---------|---------|
| Dark current | 0.0027 e-/s/pix | 0.05 e-/s/pix | 0.0067 e-/s/pix |
| Read noise | 3.1 e-/read | 12.5 e-/read | 3.7 e-/read |
| Sky (avg, F606W) | 0.040 e-/s/pix | -- | 0.042 e-/s/pix |
| Sky (avg, F200LP) | 0.094 e-/s/pix | -- | -- |
| Sky (avg, F140W) | -- | 1.17 e-/s/pix | -- |

### Color Corrections
When normalizing in a different band than the observation filter, apply the AB_nu color term for the source SED. For an E0 galaxy in F140W: AB_nu = -1.41 (effective mag is ~1.4 mag brighter than V-band). Our arc SB measurements are already in the observation filter, so no correction needed.
