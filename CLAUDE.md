# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis of the deflector galaxy in the AGEL J020613-011417 strong gravitational lens system. Two main analyses:
1. **Stellar velocity dispersion** from Keck/KCWI IFU spectroscopy using ppxf
2. **Aperture photometry + SED fitting** from HST and JWST imaging using Bagpipes

## Key Conventions

- All wavelengths in Angstroms unless otherwise noted
- Redshifts: deflector z=0.675, source z=1.302
- Coordinate system: ICRS, target at RA=31.55611, Dec=-1.23817
- IFU data cube axes: `[wavelength, y_spatial, x_spatial]` (shape ~[1273, 90, 90])
- Deflector spaxel region: `cube[:,45:64,45:55]` (19x10 spaxels)
- Noise/sky region: `cube[:,28:40,45:70]` (12x25 spaxels)
- Photometry uses AB magnitude system throughout

## Notebook Organization

There are two versions of each analysis notebook:

- **Original exploratory notebooks** (`01_IFU_spectra_extraction_and_ppxf.ipynb`, `02_Bagpipes_SED_fitting.ipynb`): Full history of the analysis with iterative attempts. Kept for reference.
- **Streamlined notebooks** (`01_streamlined_IFU_ppxf.ipynb`, `02_streamlined_Bagpipes_SED.ipynb`): Clean versions with redundant cells removed. Use these for running the analysis.

The streamlined IFU notebook uses veldis (degree=[4,30]) for the integrated spectrum and raw ppxf for per-spaxel and power-binned fitting. The `ppxf_per_spaxel()` function appends `best_fit` twice per degree iteration — account for this when indexing results (`best_fit_idx = deg_idx * 2`).

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
