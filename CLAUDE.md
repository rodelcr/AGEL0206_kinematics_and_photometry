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

## Key Results (σ_e pipeline, N=500 production)

**Primary paper number** — σ_e(<Re) on the NEW `_mtwdo_` reduction, wide arc-masked window ppxf, **with He I 3819 source-emission mask + M10 full sky-line audit (2026-05-27)**:

- **σ_e(<Re) = 269.62 ± 13.27 km/s** (asymmetric −13.45 / +13.10) at `wR3800_5400_arcmask` (rest 3800–5400 Å, Ca H+K through Mg b/Fe5335 + z=1.302 source-emission masks incl. He I 3819 + 35-entry BAD_PIXELS_REST sky/CR catalog) — **paper headline** for the M•–σ relation (Kormendy & Ho 2013 eq. 3, Greene+2020 fig. 5). **2026-06-08: was ±11.77; the D7 R_e-source systematic (±6.13) was folded into the budget (Task 5, nb15).** All numbers in this block are emitted by `scripts/paper_values.py` → `results/PAPER_VALUES.json` (single source of truth — do not hand-edit)
- Source: `scripts/run_wide_sigma_e.py --cube new_clean_hei`, caches at `results/run_wide_sigma_e/new_clean_hei/wR3800_5400_arcmask_{fsps,emiles,xsl}_T*_N500.npz`. Cube: `raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`.
- Pipeline: bad-pixel mask (M10-audited, 35 entries — 26 original CR residuals + 9 OH airglow/sky bands added via full sky-noise audit, all verified non-source via arc spectrum) + no-Balmer mask (Hδ, Hγ, Hβ kept in the fit) + ARC_MASKS_REST catalog with He I 3819 source emission added (5237–5253 Å). See `_clean_hei` presets in `scripts/run_wide_sigma_e.py`.

Systematic-budget components on the NEW reduction (post-M10, 2026-05-27):

| component | ± km/s | notes |
|-----------|--------|-------|
| stat (N=500) | 4.6 (asym −5.16/+4.13) | wild-bootstrap pool across FSPS+EMILES+XSL × 15 polynomial degrees on cleaned + He-I + M10-sky-masked new cube. **This pooled width already marginalizes over the SPS-library and polynomial-degree systematics** — between-SPS spread is ±2.04 (FSPS 271.3/EMILES 267.2/XSL 270.4), within-SPS bootstrap ±4.22, quadrature 4.73 ≈ pooled 4.64. SPS collapsed from ~26 km/s (narrow) to ~4 km/s (wide) per J5, so no separate SPS line is needed |
| I-shape (10 shapes × N=250, **M11**) | **2.27** | **rigorous re-derivation 2026-05-27** on NEW cube + M10 masks; peak-to-peak/2 of 10-shape spread (266.83–271.37). Was ±1.5 carried from old-cube sweep |
| F200 mask (3 weights × N=500, **M11**) | **6.65** | **rigorous re-derivation 2026-05-27** on NEW cube + M10 masks; (w00−w100)/2 = (269.69−256.39)/2. Bigger than the ±3.8 carried value because the cleaned NEW cube + Balmer-unmasked fit is more sensitive to arc dilution (more signal pixels → more leverage) |
| frame (vac/air, carried) | 5.0 | from prior frame-fix work; structural choice (per-SPS native frame) |
| centering (HST WCS, carried) | 4.0 | HST-derived; reduction-independent |
| fit-window (3 windows × N=500, **M11**) | **3.82** | **rigorous re-derivation 2026-05-27** on NEW cube + M10 masks. wR3800_5400_arcmask 269.66, wR4000_5400_arcmask 268.58, w6500_7500 276.22 → peak-to-peak/2 = 3.82. **Dramatically smaller than the ±15.0 carried value** because the cleaned + He-I-masked NEW cube agrees across windows to ±4 km/s (vs ±15 on the OLD cube). The narrow w6500_7500 still has large stat error (±30) but its central value sits near the wide-window numbers now |
| **reduction-pass (refined 2026-05-27 post-M10)** | **3.45** | half-Δ between cleaned + He-I + M10 new and old reductions = (269.62−262.72)/2 = 3.45. Was ±3.86 pre-M10; M10 sky masks tightened the inter-reduction gap further. Flag: only 2 reductions; refine if 3rd lands |
| **R_e-source (D7 wide, 4 estimators × N=500, 2026-06-08)** | **6.13** | **NEW — Task 5, nb15, `scripts/run_sigma_e_Re_systematic_wide.py`.** Wide-window re-measurement of D7 (was ±8.45 at narrow, never folded in). Full peak-to-peak/2 across mean 2.305"(269.62)/F140W 2.168"(267.44)/F200LP 2.441"(272.44)/CaHK+G 2.902"(279.69). User decision 2026-06-08: fold the FULL spread (CaHK+G's +10 km/s at 2.90" reflects the rising σ(R) gradient). Light-R_e-family-only spread is ±2.50 (reported as a cross-check) |
| **TOTAL (symmetric)** | **13.27** | quadrature sum (was ±11.77 pre-R_e-source; R_e-source ±6.13 folded in 2026-06-08). Computed by `scripts/paper_values.py` |
| **TOTAL (asymmetric)** | **−13.45 / +13.10** | preserves stat-side skew |

> **σ_e arc-masking-approach cross-check (Task 1, nb13, 2026-06-08):** the masking-*definition* systematic (expert/sersic/perband/global reprojected to the IFU grid) = **±5.85 km/s** — the spectroscopic analogue of the M★ ±0.16 dex term. It **overlaps the F200-mask w-sweep term (±6.65)** and is NOT added separately (larger-of-two = 6.65 kept). σ_e robust to arc-mask definition at this level. `scripts/run_sigma_e_mask_systematic.py`.

Cross-checks:
- **OLD cube cleaned + He I + M10 (`--cube headline_clean_hei`)**: σ_e(<Re) = **262.72 ± 17.99 km/s** (asym −18.10/+17.88). 6.90 km/s Δ to the headline → source of the refined ±3.45 reduction-pass systematic.
- **Pre-M10 cross-checks (cleaned + He I only)**: NEW cube 271.87 (now −2.25 above headline), OLD cube 264.16. Δ = +7.71 km/s. Adding M10 sky masks shrinks the gap to +6.90.
- **Pre-He-I cross-checks (clean only)**: NEW 268.98, OLD 260.44, Δ = +8.54.
- **Legacy (un-cleaned, no He I, no M10)**: NEW 265.76, OLD 254.85, Δ = +10.91. Full cleanup brings the cubes ~4 km/s closer together.
- **OLD cube narrow Ca H+K + G** (notebook 09d): σ_e(<Re) = 267.95 ± 30.10 km/s at `w6500_7500`. The cleaned + He-I + M10 NEW headline (269.62) is within 2 km/s of this narrow-window value.
- **M9 (DO NOT MASK)**: visible +5-7% bump at def-rest 5193–5204 Å is the Mg b LOSVD wing (NOT source emission per arc spectrum, NOT sky residual per noise spectrum). Smoke-masking it drops σ_e by 7 km/s (largest single mask effect) → it's signal, not contamination. Documented to avoid future second-guessing.

Hδ may still need targeted masking — flagged for revisit; see `METHODS_AND_SYSTEMATICS.md` Part III.5 item #0 and the inline TODO in `scripts/bootstrap_ppxf.py:_determine_goodpixels_no_balmer`.

Prior cross-checks (notebook 07c Gültekin pipeline, narrow window only — superseded by nb09 but still valid):
- §6cum cumulative I-weighted aperture ppxf: **σ_e(<Re) = 267.32 ± 24 km/s**, σ_e(<Re/2) ≈ 225.78 ± 18, σ_e(<Re/8) ≈ 209.18 ± 20
- §7 discrete Gültekin annular sum (arc-filtered to R < Re_safe = 3Re/4 = 1.72"): **256.17 −13.0/+12.7 km/s**
- §7b flat-σ extrapolation into outer annulus: 274.37 −16.2/+17.4 km/s
- nb07e arc-spectrum subtraction: matches §6cum bit-identically — residual arc dilution sub-dominant at N=500
- F200LP mask sensitivity at narrow (§6cum-nomask): σ_e(<R_e) = 250.96 ± 23 km/s (Δ = −16.36 vs §6cum headline) — see nb09 §10.1 for the equivalent test at the wide arc-masked window where the systematic shrinks to ±3.8
- Do NOT include ann5 in §7 (unfiltered gives 386 km/s, unphysical) — use §6cum or §7b

### Method choice — wide arc-masked window vs cumulative vs annular

**nb09 wide arc-masked** (current paper headline) — single ppxf fit on the I-band-weighted R<R_e aperture spectrum across rest 3800–5400 Å with explicit z=1.302 source-emission masking (Mg II 2796/2803, [O II] 3727/3729, **He I 3819**, [Ne III] 3869). Pros over nb07c §6cum: 4× more spectral pixels (2161 vs 555 good pixels) → ±6 stat vs ±24 stat; orthogonal feature-set cross-check via `wR4000_5400_arcmask` (Hβ + Mg b, no Ca H+K) gives 15 km/s window-spread systematic. Cons: needs explicit catalog of source-emission lines mapped into deflector rest frame (one-time cost; reusable for other AGEL targets).

§6cum (cumulative I-weighted aperture ppxf, nb07c, narrow window only) is preserved as the **method cross-check**. See `~/.claude/.../memory/reference_cumulative_vs_annular_sigma_e.md` for the full case. Short version of why §6cum was preferred over §7 (and why both are now superseded by nb09 for the headline):
- Single ppxf fit on the I-weighted aperture spectrum at R<R_max — matches what KH13 / Greene+20 / SAURON / ATLAS3D / MaNGA actually compute (Cappellari+2006 eq. 1)
- No binning to defend (no equal-r vs equal-N debate)
- Single LOSVD fit preserves line-shape information that §7's moment-pooling discards
- Bright-center I-weighting auto-suppresses arc contamination (verified by nb07e: < 0.1 km/s shift)
- No per-SPS V_sys anchoring choice (§7 needs one; can shift σ_e by 5–15 km/s)
- Robust to non-axisymmetry and asymmetric arc-mask gaps

§7 discrete annular sum is for:
- The σ(R) profile (radial physics, elliptical-bulge check)
- Per-annulus systematics inspection (S/N, EW, σ posterior)
- Cross-check on the cumulative number — never quote §7 as the headline

Annular binning (only relevant for §7 / σ(R)): **equal-N inside R_safe = 3R_e/4 + 1 outer flagged bin** (nb07c) is the right choice — balanced bootstrap variance per bin, narrow outer bin cleanly diagnoses arc contamination. Equal-width annuli (nb07) are kept as a sensitivity test only; the 257 (equal-width) vs 267 (equal-N) shift is at the SPS-systematic level (±24).

Historical / superseded:
- σ on 190-spaxel integrated region `cube[:,45:64,45:55]`: ~301 ± 30 km/s — **do NOT cite**; arc-contaminated
- Notebook 06 aperture posteriors at R<Re/8, Re/2, Re (the earlier pre-Gültekin path) are still in `results/final_sigma_Re_apertures.npz` but superseded by the Gültekin numbers above

Other headlines:
- Systemic redshift: z = 0.67564 (notebook 04); older results use z = 0.67511
- Effective radius (F140W + F200LP masked CoG mean): **Re = 2.305" = 16.23 kpc** (paper headline — supersedes older Re = 2.61" IFU-only value)
- log(M★/M☉) = **11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors [**11.04 ± 0.14 +0.32** at 20%] — NEW headline (2026-05-29) from principled F200LP-located + IR-extended arc masking, raw aperture photometry. One-sided +sys = Sérsic fill-in of deflector light under the arc, reaching log M★ ≈ 11.46. **Supersedes 11.33 +0.07/−0.09** (older expert-aperture, smaller masks — now a mid-range cross-check). See arc-mask section below + `NOTES_photometry_mask_systematics_2026-05-29.md`.

### Arc-mask verification (photometry side, 2026-05-29 — notebook 12)

The hand-painted F200LP arc mask is **reproducible from objective criteria** and the photometry
headline is **invariant to the masking method**. Two independent selectors in
`scripts/arc_mask_verification.py` (→ `results/arc_mask_verification.npz`, `arc_mask_bagpipes.npz`,
`results/figures/arcmask_*.png`): (A) color m_F200LP−m_F140W (arc bluer than the red deflector),
(B) Sersic-residual (subtract the 2D Sersic deflector model, flag positive excess). Results:
- F200LP **no-mask is 0.42 mag brighter than expert** (real arc contamination), but the
  **Sersic-residual mask reproduces expert to +0.016 mag**, R_e to 3%, and **Δlog M★ = +0.034**
  (≪ ±0.08 posterior). → Sersic-residual is the candidate **standard for F200LP**.
- An **S/N-regime sweep** (SNR∈{2..20}, k∈{2..8}) quantifies over-masking via the masked-flux
  fraction from the smooth model: F200LP stays clean (0.03–0.10); **F140W Sersic over-masks**
  (arc sits on bright IR deflector light) → use the **color method on F140W**.
- Caveat flagged (not fixed): `measure_Re.hst_Re` hardcodes 0.08″/pix for F200LP (cutout is 0.05″)
  — ΔR_e between masks cancels, but absolute F200LP R_e from that script is biased. See
  `NOTES_arc_mask_verification_2026-05-29.md`. **Later TODO:** spectroscopic invariance.

**4-band extension + M★ budget (2026-05-29, `scripts/photometry_systematics.py`, `NOTES_photometry_mask_systematics_2026-05-29.md`):**
- **Principled recipe:** F200LP locates the arc (best source contrast); reproject its footprint to
  every band; the deeper IR bands **extend** it (region-grow into 2-component-Sérsic-residual source
  pixels contiguous with the arc). Deflector model = **2-component (bulge+disk) Sérsic** — single
  Sérsic under-fits the bright IR galaxy (median residual +1.0) and over-masks 0.2–0.4 mag; do NOT
  use independent per-band Sérsic masks on F140W/JWST. No PSF (env lacks webbpsf) — flagged.
- **raw vs filled:** the IR-extended masks are large (F150W2 +9119 px), so **raw photutils is biased
  LOW and mask-size-dependent** (discards under-arc deflector light), while the **Sérsic fill-in is
  mask-definition-independent** (per-band vs global agree 0.01 dex). Fill-in correction +0.18–0.96 mag.
- **M★ budget** (`results/photometry_systematics_Mstar.npz`, fig `Mstar_budget.png`): headline =
  **raw, one-sided-up systematic**: log M★ = 11.16 ± 0.08 (stat,10%) +0.31 (sys) [11.04 ± 0.14 +0.32
  at 20%]; fill-in upper reach 11.46/11.36. 10%→20% shifts central ~0.10 dex. **per-band vs global**
  negligible for filled, large for raw. Scripts: `principled_mask_photometry.py` (single-Sérsic
  audit), `mask_method_comparison.py` (expert/HST-reproj/Sérsic × raw/filled/total).
- **Explicit masking-approach systematic on M★ = ±0.16 dex** (`results/Mstar_masking_systematic.npz`,
  fig `Mstar_masking_systematic.png`): peak-to-peak/2 across all approaches; dominated by under-arc
  light (raw↔filled ±0.15), mask-definition (per-band↔global) negligible ±0.004, mask-extent 0.18.
  **TODO:** the analogue masking systematic on σ_e (reproject approaches to IFU, re-run `run_wide_sigma_e.py`).
- **PSF effect on the fill — quantified (`scripts/psf_fill_model.py`, `results/psf_fill_model.npz`):**
  PSF-convolved 2-comp Sérsic (Gaussian at instrument FWHM; env lacks webbpsf) shifts the filled mag
  by **≤0.004 mag** in every band → ΔM★ ≪0.01 dex, negligible. The fill is PSF-robust (arc is outside
  the PSF core). No longer a flagged unknown.
- **R_e pixscale fix:** `measure_Re.hst_Re` now reads pix scale from the WCS (F200LP was wrongly
  0.08″, is 0.05″) → diagnostic F200LP R_e 3.05→1.91″. **Headline R_e=2.305″ unaffected** (from
  `final_sigma_e.py`, already WCS-correct). Flag: the two CoG algorithms differ (measure_Re 1.91 vs
  final_sigma_e 2.52 for F200LP) — separate pre-existing methodology gap.

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
