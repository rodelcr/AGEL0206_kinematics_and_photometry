# AGEL0206 Deflector Galaxy: Stellar Kinematics and Photometric Analysis

Analysis of the deflector galaxy in the **AGEL J020613-011417** strong gravitational lens system. This repository contains tools and notebooks for measuring the stellar velocity dispersion from Keck/KCWI integral field unit (IFU) spectroscopy, performing aperture photometry on HST and JWST imaging, and estimating the stellar mass via SED fitting with [Bagpipes](https://bagpipes.readthedocs.io/).

## Target

| Property | Value |
|----------|-------|
| Name | AGEL J020613-011417 (AGEL0206) |
| RA, Dec (ICRS) | 31.55611, -1.23817 |
| Deflector redshift | z = 0.675 |
| Source redshift | z = 1.302 |

## Repository Structure

```
.
├── notebooks/
│   ├── 01_IFU_spectra_extraction_and_ppxf.ipynb   # IFU spectral extraction + ppxf kinematics
│   ├── 02_Bagpipes_SED_fitting.ipynb               # Broadband SED fitting with Bagpipes
│   └── 03_Prospector_SED_fitting.ipynb              # Exploratory SED fitting (alternative)
├── scripts/
│   ├── photometry_masking_HST.py                    # Interactive HST aperture photometry tool
│   └── photometry_masking_JWST.py                   # Interactive JWST aperture photometry tool
├── example_outputs/                                  # Sample parameter JSON files from photometry
├── figures/                                          # ppxf fitting result plots
├── *.dat                                             # Filter transmission curves (HST + JWST)
├── requirements.txt
└── README.md
```

## Data Requirements

The following data files are **not included** in this repository (too large for git). Place them in the repository root directory before running the notebooks.

### IFU Data
- `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` — KCWI IFU data cube (~253 MB). Obtained from Keck Observatory Archive.

### HST Imaging (from [MAST](https://mast.stsci.edu/))
- `AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits` — WFC3/UVIS F200LP cutout
- `AGEL020613-011417A_F200LP_WFC3_drc_sci.fits` — Full-field science image
- `AGEL020613-011417A_F200LP_WFC3_drc_wht.fits` — Weight map
- `AGEL020613-011417A_F140W_WFC3_cutout_L3.fits` — WFC3/IR F140W cutout
- `AGEL020613-011417A_F140W_WFC3_drz_sci.fits` — Full-field science image
- `AGEL020613-011417A_F140W_WFC3_drz_wht.fits` — Weight map

### JWST Imaging (from [MAST](https://mast.stsci.edu/), Program 05594)
- `jw05594-o101_t103_nircam_clear-f150w2_i2d.fits` — NIRCam F150W2
- `jw05594-o101_t103_nircam_clear-f322w2_i2d.fits` — NIRCam F322W2

### Stellar Template Libraries
- `TEXT/` directory — 1,273 reference stellar spectra for veldis/ppxf template fitting
- `spectra_emiles_9.0.npz`, `spectra_fsps_9.0.npz`, `spectra_xsl_9.0.npz` — Stellar population synthesis model libraries for ppxf

## Notebooks

### 01 - IFU Spectra Extraction and ppxf Fitting

Loads the KCWI IFU data cube, extracts 1D spectra from the deflector galaxy region, and fits stellar kinematics using both [veldis](https://github.com/kvgc153/veldis) (a ppxf wrapper) and [ppxf](https://pypi.org/project/ppxf/) directly. Includes:
- Integrated deflector spectrum fitting
- Per-spaxel velocity dispersion mapping
- Power-binned (adaptive S/N binning) kinematic analysis

### 02 - Bagpipes SED Fitting

Fits broadband photometry (HST F200LP, F140W + JWST F150W2, F322W2) using Bagpipes with an exponential-tau star formation history and Calzetti dust attenuation. Key result: **log(M\*/M_sun) = 11.34 +0.06/-0.08**.

### 03 - Prospector SED Fitting (Exploratory)

Early/alternative SED fitting attempts. Retained for reference; see notebook 02 for the primary analysis.

## Photometry Scripts

The `scripts/` directory contains interactive aperture photometry tools for HST and JWST images. Both scripts provide:
- Matplotlib slider GUI for adjusting elliptical aperture and background annulus geometry
- Click-to-mask interface for contaminating sources
- Noise map support with full error propagation
- AB magnitude and flux density output (Jy, erg/s/cm2/Hz)
- Save/load of aperture parameters (JSON) and masks (FITS)

**Usage:**
```bash
# Edit the 'fname' variable near the top of the script to select your target FITS file
python scripts/photometry_masking_HST.py
python scripts/photometry_masking_JWST.py
```

Three execution modes are supported:
1. **Reload** — Load existing parameter/mask files and recalculate photometry
2. **Import** — Import parameters/masks from another filter and reproject to match
3. **Interactive** — Launch the GUI to define aperture and mask from scratch

## Dependencies

```bash
pip install -r requirements.txt
```

Bagpipes requires a Bayesian sampling backend. Install one of:
- [MultiNest](https://github.com/JohannesBuchner/PyMultiNest) (recommended, requires C library)
- [Nautilus](https://github.com/johannesulf/nautilus) (`pip install nautilus-sampler`, pure Python fallback)

## Known Issues

- The Bagpipes notebook references `'JWST_NIRCAM.F322W2.dat'` (uppercase NIRCAM) while the actual filename is `JWST_NIRCam.F322W2.dat`. This works on macOS (case-insensitive filesystem) but will break on Linux. Match the case if deploying on a case-sensitive system.
