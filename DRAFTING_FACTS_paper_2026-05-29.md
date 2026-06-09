---
title: "AGEL J020613−011417 — Drafting Fact Sheet for the Methods & Results Sections"
subtitle: "Facts, decisions, and provenance keyed to the paper outline — last revised 2026-06-08 (post-M12)"
author: "Compiled for R. Córdova Rosado from project notes, memory, METHODS_AND_SYSTEMATICS.md, and TESTS_AND_DIAGNOSTICS.md"
date: "2026-06-08"
---

# How to use this document

This is a **drafting reference**, not prose for the paper. Every bullet is a fact or
decision we have already made, with its source file in `[brackets]` so you can verify it.
Where a number is superseded, the current value is given and the stale one flagged.
**Three gaps** are flagged up front — they are not in this repo and you must source them
elsewhere before drafting those paragraphs:

- **G1 — Lens modeling (Ferrami et al.).** No BPL / point-mass / M(<r_break) / joint-posterior /
  PSF / shear / subhalo content exists in this kinematics+photometry repo. Source from the
  lens-modeling repo (`~/Documents/AGEL/202509_DESJ0206_modeling/` or
  `20251112_DESJ0206_Pyauto_PRONTO/`) or the Ferrami draft directly.
- **G2 — X-ray / "truly quiescent" statement.** No X-ray, Chandra, or explicit quiescence
  data anywhere in the read files. The galaxy is only ever described as a "passive elliptical."
  You need an external source if you want the X-ray quiescence claim in §Results.
- **G3 — "We do not account for peculiar velocities."** This sentence is not currently written
  in any methods doc. It is a true statement about our analysis (we use the line-fit systemic z
  directly) — just add it; nothing to reconcile.

**Headline numbers (current, 2026-06-08, post-M12; emitted by `scripts/paper_values.py` → `results/PAPER_VALUES.json`):**

- σ_e(<R_e) = **269.62 ± 13.27 km/s** (asym −13.45 / +13.10); often rounded **270 ± 13 km/s**. **2026-06-08: was ±11.77; the D7 R_e-source systematic (±6.13) was folded into the budget (M12, Task 5, nb15).**
- log(M⋆/M☉) = **11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors [**11.04 ± 0.14 (stat) +0.32 (sys)** at 20%] — **NEW principled arc masking** (F200LP-located + IR-extended, raw aperture). The one-sided +sys is the Sérsic fill-in of deflector light hidden under the arc, reaching log M⋆ ≈ 11.46. **Supersedes 11.33 +0.07/−0.09** (older expert-aperture value, smaller masks; now a mid-range cross-check). See §2.1.1b / §3.2 and `NOTES_photometry_mask_systematics_2026-05-29.md`.
- R_e = **2.305″ = 16.23 kpc**.
- z_deflector = **0.67564**; z_source = **1.302**.

**Cosmology:** FlatLambdaCDM, H₀ = 70 km/s/Mpc, Ω_m = 0.3; 1″ ≈ 7.04 kpc (7043.5 pc/arcsec),
D_A = 1452 Mpc at z = 0.67564. [reference_hst_jwst_data_properties.md]

> **Staleness warning:** the *body* of `METHODS_AND_SYSTEMATICS.md` (dated 2026-05-18) still
> prints pre-M11 numbers (254.85 / 268.98, total ±18) in its narrative prose. Its own §0 table
> and `CLAUDE.md` carry the current headline **269.62 ± 13.27** (post-M12, 2026-06-08). Always
> quote the current value (canonical source: `results/PAPER_VALUES.json`).

---

# §1 Data and observations

*All values below are read directly from the FITS headers of the analysis files (verified
2026-06-08), not from memory. Target: **AGEL J020613−011417** = DES J0206−0114, ICRS
RA = 31.55611°, Dec = −1.23817° (HST `RA_TARG`/`DEC_TARG`); deflector z = 0.67564, source z = 1.302.
The field is the massive cluster **ACT-CL J0206.2−0114** (ACT DR5, Hilton et al. 2021) — the JWST
pointing target.*

| Facility | Inst. / mode | Band(s) | Program (PI) | UT date | Exp. | Pixel scale |
|---|---|---|---|---|---|---|
| Keck II | KCWI — KCRM red (RL) **+** KCB blue (BL), Medium slicer, 2×2 | RL 5625–8941 Å (σ_e); BL λc 4500 Å | **U204 + K409 + U002** (3-night combine) | 2024-12 → 2025-11 | **≈180 min RED + ≈198 min BLUE** (36 × 300 s RED) | 0.300″/spaxel |
| HST | WFC3/UVIS | F200LP | **16773 (Glazebrook)** | 2022-07-14 | 600.0 s (NDRIZIM 4) | 0.050″/pix |
| HST | WFC3/IR | F140W | **16773 (Glazebrook)** | 2022-07-14 | 597.7 s (NDRIZIM 3) | 0.080″/pix |
| JWST | NIRCam SW, module B | F150W2 (CLEAR) | **05594 (Mahler)** | 2024-08-27 | 1836.0 s | 0.03075″/pix |
| JWST | NIRCam LW, module B | F322W2 (CLEAR) | **05594 (Mahler)** | 2024-08-27 | 1836.0 s | 0.06301″/pix |

### §1.1 Keck/KCWI IFU spectroscopy (the σ_e data)

- **Keck II / KCWI**, KCRM **red** arm, **RL (Red-Low)** grating (`RGRATNAM`), **Medium** slicer
  (`IFUNAM`), **2×2** binning (`BINNING`/`CCDSUM`); nod-and-shuffle **off** (`RNASNAM='Open'`).
- **The headline cube is a 3-NIGHT, 3-PROGRAM combine** (see §2.2.1), NOT a single dataset — the
  cube header (`DATE-OBS` 2025-11-17, `PROGPI` Jones, `XPOSURE` 300 s) describes only the **first
  input frame**, so cite the combine, not the header:

  | Night (UT) | Program (PI) | RED (RL) | BLUE (BL) |
  |---|---|---|---|
  | 2024 Dec 29 | U204 | 4 × 300 s | 1 × 1320 s *(excluded by reducer)* |
  | 2025 Aug 30 | K409 (PI TBD) | 12 × 300 s | 4 × 990 s |
  | 2025 Nov 17 | U002 (Jones) | 20 × 300 s | 5 × 1320 s |

  → **36 × 300 s ≈ 180 min RED + ≈ 198 min BLUE on-source.** The NEW **`_mtwdo_`** reduction
  (Master-Twilight + Dome hybrid flats, K. R. Gupta), `Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`;
  dithered (PA 0°/45°/90°); CR rejection via **astroscrappy** (`fsmode=median`, `psffwhm=2.50`).
  `OBSERVER` = Glazebrook, Gupta, Jones, Kacprzak, Alcorn, Rhoades, Barone, Tran, Chen, Vasa.
  *(K409 PI + Aug-30 DIMM seeing still TBD for acknowledgements.)*
- **Both KCWI arms recorded simultaneously.** The **red** arm (**KCRM**, RL grating, central λ
  7150 Å `RCWAVE` — the cube below) is the σ_e dataset. The **blue** arm (**KCWI/KCB**, **BL
  (Blue-Low)** grating `BGRATNAM`, central λ **4500 Å** `BCWAVE`, **KBlue** filter `BFILTNAM`, N&S
  off) is the ≈198 min BLUE above; it is **not used for σ_e** (the deflector's Ca H+K…Mg b absorption
  falls in the red) but supplies the **source-z = 1.302 cross-check** (Fe II UV multiplet; §2.2.2).
  **The blue combined cube is not in this repo** (only the RL cube) — quote the blue arm from the
  observing log, not a local file.
- **Conditions (Nov-17 frame, header-literal):** airmass **1.165** (`AIRMASS`); **guide-star FWHM
  1.271″** (`GUIDFWHM`) — the adopted seeing; parallactic angle 15.6°, rotator 0°. (Per-night
  airmass/seeing for the other two nights not carried in the combined cube.)
- **Cube:** 100 × 100 spaxels × 3317 λ-pixels; **0.300″/spaxel** (square; from WCS), **30″ × 30″**
  FOV; wavelength **5625–8941 Å** at **1.0 Å/pix** (`CRVAL3`/`CD3_3`), red central λ 7150 Å
  (`RCWAVE`). Spectral resolution **R ≈ 10000** at the fit band → **σ_inst ≈ 12.6 km/s** (well below
  σ_e ≈ 270; the LSF is carried into ppxf). [reference_kcwi_data_properties.md]
- Pointing `RA`/`DEC` = 02:06:13.38 / −01:14:20.8; target `TARGRA`/`TARGDEC` = 02:06:13.58 / −01:14:19.8.

### §1.2 HST/WFC3 imaging (F200LP, F140W)

- **Program 16773, PI Karl Glazebrook** (`PROPOSID`/`PR_INV_*`); `TARGNAME` = DESJ0206-0114; both
  bands **UT 2022-07-14** (`EXPSTART` MJD 59774.30–59774.31), ORIENTAT 111.9°, drizzled, units
  **ELECTRONS/S**.
- **F200LP** — WFC3/**UVIS**, aperture UVIS2: EXPTIME **600.0 s**, NDRIZIM 4, NCOMBINE 2, drizzle
  scale **0.050″/pix**; `PHOTFLAM` = 5.1851×10⁻²⁰, `PHOTPLAM` = 4923.48 Å → **ZP_AB = 27.344**.
- **F140W** — WFC3/**IR**: EXPTIME **597.7 s**, NDRIZIM 3, NCOMBINE 3, drizzle scale **0.080″/pix**;
  `PHOTFLAM` = 1.4829×10⁻²⁰, `PHOTPLAM` = 13922.91 Å → **ZP_AB = 26.446**.
- AB zeropoints are recomputed per-image from `PHOTFLAM`/`PHOTPLAM` (never hardcoded);
  `ZP_AB = −2.5·log10(PHOTFLAM) − 21.10 − 5·log10(PHOTPLAM) + 18.6921`. **Drizzle pixel correlation:**
  empirical noise ≈1.3–2× the CCD-equation value → motivates the 10% flux floor (§2.1.4).

### §1.3 JWST/NIRCam imaging (F150W2, F322W2)

- **Program 05594** — *"JWST Cluster SLICE — Strong Lensing and Cluster Evolution"*, **PI Guillaume
  Mahler** (SURVEY category); `TARGPROP` = ACT-CL_J0206.2−0114 (the cluster field; the deflector
  sits in it). Both bands **UT 2024-08-27** (`DATE-BEG` 06:02–06:40), NIRCam **module B**, units
  **MJy/sr** (i2d), NINTS 1 / NGROUPS 4.
- **F150W2** (SHORT channel, CLEAR pupil): `EFFEXPTM`/`XPOSURE` **1836.0 s**; `PIXAR_SR` =
  2.2222×10⁻¹⁴ sr → **0.03075″/pix** (PIXAR_A2 = 9.454×10⁻⁴ ″²) → **ZP_AB = 28.033**.
- **F322W2** (LONG channel, CLEAR pupil): `EFFEXPTM`/`XPOSURE` **1836.0 s**; `PIXAR_SR` =
  9.3307×10⁻¹⁴ sr → **0.06301″/pix** (PIXAR_A2 = 3.970×10⁻³ ″²) → **ZP_AB = 26.475**.
- JWST AB from area: `ZP_AB = −6.10 − 2.5·log10(PIXAR_SR)`. No PAM correction (drizzled/resampled).

---

# §2 Methods

## §2.1 Photometric Analysis

### §2.1.1 Masking

- **Use the F200LP `_mask.fits` as the arc mask — NOT the F140W one.** [feedback_masking_strategy.md; reference_hst_jwst_data_properties.md; METHODS §II.4]
  - F200LP mask (2512 HST pixels) is hand-tuned to cover **the arc only** and leaves the
    **deflector core unmasked** — it is the original arc-tuned segmentation mask.
  - F140W mask (1121 HST pixels) covers **the arc PLUS the deflector core** (the photometry
    tool's aperture geometry pulled the core in); the two masks have **zero spatial overlap**.
  - Consequence: F140W is valued for its *image* (sharper Sérsic fit), never for its mask.
    Using the F140W mask as an arc mask would destroy the deflector-core signal.
### §2.1.1b Principled (reproducible) arc masking — NEW 2026-05-29

The hand-painted masks are now backed by an objective, reproducible recipe (so a referee can
reproduce the arc selection), validated + applied to all four bands.
[`NOTES_arc_mask_verification_2026-05-29.md`, `NOTES_photometry_mask_systematics_2026-05-29.md`;
`scripts/arc_mask_verification.py`, `mask_method_comparison.py`, `photometry_systematics.py`]

- **Locator = F200LP** (highest contrast on the blue z=1.302 source). Two independent objective
  selectors reproduce the F200LP hand mask: (i) a color cut m_F200LP−m_F140W (arc bluer than the
  red deflector); (ii) a **Sérsic-residual** mask (subtract a 2D deflector model, flag positive
  excess > 3σ_sky). On F200LP the Sérsic-residual mask reproduces the expert photometry to
  **0.016 mag**, R_e to 3%.
- **Transfer to other bands by reprojecting the F200LP footprint** (WCS + `map_coordinates`).
  Reproduces the expert JWST aperture mags to **0.01–0.02 mag**. Do **NOT** use an independent
  per-band Sérsic-residual mask on F140W/JWST — a single Sérsic under-fits the bright extended
  IR galaxy (median residual +1.0) and over-masks by 0.2–0.4 mag.
- **IR extends the source footprint:** the deeper/sharper JWST bands see the arc's fuller extent
  (F150W2 mask grows ~6× vs HST). We grow the F200LP arc into each band's 2-component-Sérsic
  residual source pixels that are contiguous with the arc (region growing; cannot grab the core
  pedestal or field companions).
- **Deflector model = 2-component (bulge+disk) Sérsic** — single Sérsic under-fits JWST; 2-comp
  cuts galaxy-body RMS 3.3–3.5×. **PSF effect on the fill quantified** (`scripts/psf_fill_model.py`):
  a PSF-convolved 2-comp Sérsic (Gaussian at instrument FWHM; env lacks `webbpsf`) shifts the filled
  mag by **≤0.004 mag** → ΔM⋆ ≪0.01 dex, negligible. The fill is PSF-robust (arc is outside the PSF core).
- **raw vs fill-in:** raw photutils (masked) discards the deflector light under the arc; the
  Sérsic fill-in (replace masked pixels with the model) recovers it. With the large IR-extended
  masks the fill-in correction is **+0.18 to +0.96 mag** per band. **filled photometry is
  mask-definition-independent** (per-band vs global agree to 0.01 dex); **raw is mask-size-
  dependent** (up to 0.78 mag) — see §3.2.
- **per-band vs global mask** (union of all per-band masks): negligible for filled (0.01 dex M⋆),
  large for raw — confirming filled is the robust quantity, raw the conservative lower bound.

- **Files** (in sibling dir `../velocity_dispersion_from_IFU/`):
  `AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits` (500×500, 25″×25″);
  the mask was revised 2026-02-28 (`_mask_20260228.fits`); original kept as `_mask_OLD.fits`.
- **JWST bands:** same arc geometry; the HST F200LP mask is re-projected from the HST WCS
  onto each JWST frame. [METHODS §II.4]
- **For the IFU/σ_e application** (context): F200LP mask reprojected to the 0.30″ IFU grid via
  `scipy.ndimage.map_coordinates(order=0)` nearest-neighbor; ~38/184 spaxels inside R_e flagged
  (~21% by count, ~27% by I-weight). [feedback_masking_strategy.md]

### §2.1.2 Estimating R_e

- **Headline R_e = 2.305″ = 16.23 kpc** at z = 0.67564. [CLAUDE.md; METHODS §II.3/A3]
- **Method:** mean of the **F140W and F200LP masked curves-of-growth**
  (`scripts/final_sigma_e.py:curve_of_growth`, the PRODUCTION script):
  - F140W masked CoG = **2.168″**; F200LP masked CoG = **2.441″**; mean = **2.305″**.
  - Masked CoG zeros mask-flagged HST pixels before the radial flux integral, so the lensed
    arc / companions do not bias the half-light point.
  - Centers: HST-mean `photutils.centroid_2dg` centers (same as the σ_e pipeline), both bands.
- **Do NOT confuse with `scripts/measure_Re.py`** → `results/Re_measurements.npz`. That script
  uses the image-geometric center and produces **10 variants** (3 sources × 4 masking strategies:
  `unmasked`, `zeroed`, `proper`, `PSF-convolved`); its "proper-mask" mean (2.633″) is **not the
  headline** — early-reference only. (`measure_Re.hst_Re` previously hardcoded pixscale=0.08 for
  both bands; **fixed 2026-06-08 (A3c/nb14) to read the WCS pixscale + `r_max_arcsec=6.0`, now
  reconciling with the production CoG to <0.04″** — but `final_sigma_e` remains the headline path.)
- **Superseded prior value:** IFU white-light R_e = 2.61″ (nb05/06 era) — superseded because the
  HST-based 2.305″ has the F140W+F200LP arc masks built in. A Ca H+K+G depth-map gives 2.866″.
- **R_e as a σ_e systematic (test D7):** four R_e definitions shift σ_e (rises monotonically with
  R_e). Re-measured at the wide window (M12, nb15): spread **±6.13 km/s**, now **folded into the
  formal budget** (§3.1, §2.4.3); the narrow-window 16.9 km/s figure is superseded. [CLAUDE.md]

### §2.1.3 Measuring Fluxes

- **Tools:** `photutils` aperture photometry (interactive masking GUIs
  `scripts/photometry_masking_{HST,JWST}.py`). **`synphot` is NOT currently used** — flux→AB
  is done directly from header zeropoints. *(If the outline wants synphot, that is aspirational;
  flag with co-authors.)* [METHODS §II.4]
- **AB zeropoints — read per-image from FITS headers, never hardcoded.** [CLAUDE.md; reference_hst_jwst_data_properties.md]
  - HST (drizzled e⁻/s):
    `ZP_AB = −2.5·log10(PHOTFLAM) − 21.10 − 5·log10(PHOTPLAM) + 18.6921`.
  - JWST (i2d, MJy/sr): `ZP_AB = −6.10 − 2.5·log10(PIXAR_SR)`.
- **Four bands / instruments / programs:**

| Telescope | Instrument | Filter | Program (PI) | Drizzled scale | EXPTIME | ZP_AB | AB (aperture) |
|---|---|---|---|---|---|---|---|
| HST | WFC3/UVIS | F200LP | 16773 (Glazebrook) | 0.0500″/pix | 600.0 s | 27.344 | 22.613 |
| HST | WFC3/IR | F140W | 16773 (Glazebrook) | 0.0800″/pix | 597.7 s | 26.446 | 19.134 |
| JWST | NIRCam SW | F150W2 | 05594 (Mahler) | 0.03075″/pix | 1836.0 s | 28.033 | 18.942 |
| JWST | NIRCam LW | F322W2 | 05594 (Mahler) | 0.06301″/pix | 1836.0 s | 26.475 | 18.604 |

  [reference_hst_jwst_data_properties.md]
- JWST pointing comes from the cluster program ACT-CL J0206.2−0114 (ACT DR5, Hilton et al. 2021);
  the deflector sits in that field. No PAM correction (drizzled). **Drizzle pixel correlation:**
  empirical noise ~1.3–2× the idealized CCD-equation value — motivates the 10% flux floor (§2.1.4).
- **2D Sérsic verification given masking (the morphology cross-check; notebook 08):** single-component
  2D Sérsic per band (`scripts/measure_Re.py` framework + the **2026-05-11 bound-fix** n∈[1.0,6.0],
  ellip∈[0.0,0.6], 3-init grid — prevented an n=0.30 flat-disk escape on F200LP), extrapolated to total
  flux:

| Filter | AB aperture | AB Sérsic-total | Δmag | n | r_eff (″) |
|---|---|---|---|---|---|
| F200LP | 22.613 | 20.672 | −1.94 | 1.42 | 2.49 |
| F140W  | 19.134 | 18.282 | −0.85 | 1.54 | 1.88 |
| F150W2 | 18.942 | 18.154 | −0.79 | 1.40 | 1.72 |
| F322W2 | 18.604 | 17.633 | −0.97 | 1.97 | 2.03 |

  Aperture under-counts outer-profile flux by 0.8–1.9 mag, but the integrated-mass effect is only
  +0.07 dex (see §2.1.4 / Results). Caches: `results/sersic_total_photometry.npz`,
  `results/bagpipes_sersic_refit.npz`. [METHODS §II.6]

### §2.1.4 SED Modeling with `bagpipes`

- **Tool:** `bagpipes` (Carnall et al. 2018, 2019); BPASS+CLOUDY templates (Bagpipes default).
  Notebook `02_streamlined_Bagpipes_SED.ipynb`. [METHODS §II.5; reference_sed_fitting.md]
- **Model:** exponential-τ SFH; Calzetti dust (A_V free, 0–2); free metallicity (0–2.5 Z☉);
  redshift prior (0.674, 0.676) bracketing the line-fit z. Sampler **Nautilus** (n_live=400),
  500-sample posterior.

  ```python
  fit_instructions = {
    "redshift":   (0.674, 0.676),
    "exponential": {"age": (0.1,15.), "tau": (0.3,10.),
                    "massformed": (1.,15.), "metallicity": (0.,2.5)},
    "dust": {"type": "Calzetti", "Av": (0., 2.)},
  }
  ```
- **Why a 10% fractional error floor on each band:** a conservative, uniform systematic budget
  that covers both the formal photon noise **and** the drizzle pixel-correlation underestimate
  (empirical noise ~1.3–2× the idealized CCD-equation value). [METHODS §II.4/II.7; reference_sed_fitting.md]
- **Systematic test to fill in masked pixels = the 2D-Sérsic-total refit:** because the aperture
  excludes arc-masked pixels and truncates the outer profile, we fit a 2D Sérsic per band
  (§2.1.3), extrapolate to total flux, and re-run Bagpipes with identical priors. This recovers
  the flux missing from masked/truncated regions and bounds the mass bias.
  Result: log(M⋆/M☉) = **11.40 +0.11/−0.15** (Sérsic-total) vs **11.33 +0.07/−0.09** (aperture)
  → aperture is biased low by ~0.07 dex (+16% in M⋆), within the quoted uncertainty. [METHODS §II.6]
  **NOTE (2026-05-29):** this analytic-total refit is now superseded by the per-aperture 2-component
  Sérsic **fill-in** under the principled IR-extended masks (§2.1.1b / §3.2): the headline is the raw
  value 11.16 (10%) with a one-sided +0.31 dex fill-in systematic reaching 11.46. The 11.40
  analytic-total sits inside that bracket.
- **Posteriors (older expert-aperture, now a cross-check):** log(M⋆/M☉)=11.33 +0.07/−0.09; mass-weighted age 5.10 Gyr
  [3.69–6.11]; τ = 0.68 Gyr; Z = 0.89 Z☉; A_V = 0.82 mag; z_phot = 0.675 [0.674–0.676].
  Per-band residuals: F200LP +4.7%, F140W −11.3%, F150W2 −8.0%, F322W2 +9.6% (mild single-τ-SFH
  NIR-slope tension; not catastrophic, within ±0.07–0.15 dex). [METHODS §II.5]
- **Caches / figure:** `results/bagpipes_sed_results.npz` (the file Figure 2 loads),
  `pipes/posterior/AGEL0206/0206_real.h5` (delete to refit), `results/bagpipes_sersic_refit.npz`.

---

## §2.2 Spectroscopic Analysis

### §2.2.1 Calibration — @Keerthi / Kaustubh

- **Instrument:** Keck II / KCWI, **KCRM red arm**, Red-Low (RL) grating, Medium slicer, 2×2
  binning; dithered at PA 0°/45°/90°. [METHODS §I.1; reference_kcwi_data_properties.md]
- **Cube (headline):** `Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits` — the NEW
  **`_mtwdo_`** reduction (Master Twilight + Dome hybrid flats), K. R. Gupta's latest re-reduction.
  The earlier no-`_mtwdo_` cube is the OLD reduction (second-reduction cross-check; co-aligned
  to <0.002″, differs only in flat-fielding; σ_e differs ~7 km/s — the "reduction-pass" systematic).
- **Reduction chain:** raw frames → `kcwidrp`; stacked with `ISMgas/kcwi/kcwiReduction` by
  **K. R. Gupta (Kaustubh)**; Dec-29 RED frames pre-reduced by Yuguang Chen. *(Calibration prose
  is Keerthi/Kaustubh's to write — these are the facts on hand.)*
- **Cube geometry:** shape (3317, 100, 100); **0.300″ spaxels**; ~30″×30″ FoV; λ = 5625–8941 Å at
  Δλ = 1.0 Å/pix; BUNIT = FLAM16/arcsec².
- **Three confirmed nights** (totals ≈ **180 min RED + 198 min BLUE** on source):

| Night (UT) | Program (PI) | RED | BLUE |
|---|---|---|---|
| 2024 Dec 29 | U204 | 4×300 s | 1×1320 s (excluded by reducer) |
| 2025 Aug 30 | K409 (PI TBD) | 12×300 s | 4×990 s |
| 2025 Nov 17 | U002 (Jones) | 20×300 s | 5×1320 s |

  **FLAG:** the cube header (DATE-OBS 2025-11-17, PROGPI Jones) describes only the first input
  frame — cite the multi-night combine, not the header. K409 PI and Aug-30 DIMM seeing still TBD.
- **Spectral resolution / LSF:** FWHM_inst = 2.355 × DISPSCAL(0.2940) × 1.0 Å = **0.692 Å**
  → **R ≈ 10000** at the fit band; **σ_v,inst ≈ 12.6 km/s observed (7.5 km/s rest-frame)** — the
  ~270 km/s dispersion is resolved by ~20×. LSF robustness ≤ 0.83 km/s over DISPSCAL ×0.5–×2.0.
- **Seeing:** guide-camera FWHM = **1.27″** (`GUIDFWHM`, Nov-17; conservative).

### §2.2.2 Redshift Measurement

- **Method — single/multi-line Gaussian fits on absorption (+emission) features** in the
  integrated spectrum, on **NIST air rest wavelengths** (`scripts/redshift_verify.py`, notebook 04):
  - Ca K 3933.66, Ca H 3968.47 (fit as a doublet with shared width), G-band 4304.40, Hδ, Hβ,
    plus the [O II] 3726/3729 doublet (shared z, σ, flux ratio). Per-line z inverse-variance combined.
    Hγ excluded (not trusted); Mg I / Fe / Na D out of range or blended.
  - **Deflector systemic z = 0.67564** (supersedes AGEL DR2 z = 0.67511; ~95 km/s offset, used only
    in old notebooks). [METHODS §I.2; reference_redshift_verification.md]
- **Double verification with ppxf:** the Gaussian line-fit z is independent of, and consistent with,
  the ppxf V_sys (ppxf returns z via z = (1+z₀)·exp(V/c) − 1; Cappellari 2023 eq. 5c).
- **Source z = 1.302** (AGEL DR2), cross-checked via Fe II UV resonance multiplet + Mg II 2796/2803
  in the KCRM blue-arm cube. *(Measured z_source ± uncertainty + its cache are a `[TBD]` placeholder
  — confirm before quoting a source-z error bar.)*
- **Peculiar velocities:** we **do not** correct for peculiar velocities — we use the line-fit
  systemic z directly (see gap G3; this sentence needs to be added to the draft).
- **Wavelength convention:** AIR rest wavelengths throughout (NIST-verified); KCWI native vacuum is
  converted to air. Do not mix SDSS vacuum values (~1 Å ≈ 80 km/s systematic at 4000 Å).
  [feedback_wavelength_convention.md]

### §2.2.3 Stellar Kinematics & Velocity Dispersion

**Spatial regime / masking / R_e**

- Headline σ_e is measured on a **single I(r)-weighted 1-D aperture spectrum at R < R_e** — the
  literature-standard σ_e (Cappellari et al. 2006 eq. 1; the KH13 / Greene+20 / SAURON / ATLAS3D /
  MaNGA convention). **R_e = 2.305″ = 16.23 kpc** (§2.1.2).
- **Center:** HST-mean `photutils.centroid_2dg` of F140W+F200LP cores, propagated through the KCWI
  WCS; adopted RA = 31.55613°, Dec = −1.23819° (ICRS). F140W↔F200LP offset 0.36″ (F200LP shifted by
  the UV-bright arc; F140W is the clean bulge centroid).
- **Spatial arc mask:** F200LP segmentation mask reprojected to the IFU grid (nearest-neighbor);
  38/184 spaxels inside R_e flagged (~27% by I-weight). F140W mask not used (covers the core).
- **I-weight map:** per-spaxel collapse of the IFU 6500–7500 Å white-light band; spectra **summed**
  (not averaged) so the wild bootstrap preserves per-spaxel noise.
- R < R_e/8 (JFK95 σ_c) was **dropped** — at 0.288″ it sits inside seeing FWHM/2 = 0.64″ (~3 spaxels).
  No central σ is quoted; σ_e(<R_e/2) ≈ 225 km/s is reported as a gradient diagnostic only.

**ppxf setup**

- **ppxf v9.x** (Cappellari & Emsellem 2004; Cappellari 2017, 2023), **moments = 2** (first two
  velocity moments V, σ). [METHODS §I.5/I.6; scripts/bootstrap_ppxf.py]
- **Polynomial degrees:** additive Legendre `degree` swept over **15–29 (15 values)** and pooled;
  `mdegree = 0`. σ-vs-degree is flat within the bootstrap envelope (no LOSVD absorption by the
  polynomial). *(Note: METHODS §I.5 prose mis-labels these "multiplicative" — the code passes
  `degree=…, mdegree=0`, i.e. additive. Describe them as additive.)*
- **Three SPS libraries, pooled at the posterior level:** FSPS, EMILES, XSL (8–9 templates each over
  rest 3500–5000 Å).
- **Per-SPS native wavelength frame:** FSPS = vacuum, EMILES = air, XSL = air. Galaxy frame matched
  per-SPS (scalar-median vac↔air via Ciddor 1996, following Cappellari's pattern). Why: the legacy
  "convert galaxy to air for all three" produced a spurious −90 km/s FSPS V_sys; the fix collapsed
  the V_sys split-track ~110 → ~15–18 km/s (σ shift ≤2 km/s). [reference_ppxf_vacair_handling.md]
- **Per-SPS V_sys subtracted before pooling** (σ is V-invariant; V_sys medians FSPS −19, EMILES −4,
  XSL −1 km/s). A single global V_sys would inflate σ_e by ~10–15 km/s in the discrete cross-check.
- **Emission-line handling (sigma clipping / removing emission lines):**
  - Custom **no-Balmer goodpixels** (`_determine_goodpixels_no_balmer`): masks only the forbidden
    lines ([O II], [O III], [O I], [N II], [S II]) and **keeps Hδ, Hγ, Hβ IN the fit** (they are
    absorption in this passive deflector; masking them discards ~120 stellar pixels). **Hδ: keep
    unmasked — RESOLVED (M12/nb16):** Hδ is well-fit (local-MAD not an outlier), and masking it
    shifts σ_e by +6–8 km/s = LOSVD *information*, not contamination (the M9 "it's signal" pattern);
    the `bootstrap_ppxf.py` TODO is closed.
  - **z = 1.302 source-emission masks** (`ARC_MASKS_REST`, mapped by (1+z_s)/(1+z_l)=1.374):
    Mg II 2796/2803 (def-rest 3835–3855 Å), [O II] 3727/3729 (5115–5135), [Ne III] 3869 (5260–5340),
    **He I 3819** (5237–5253, added 2026-05-27 as test M8).
  - **Tellurics — handled by masking, not correction.** One band is masked: the **O₂ A-band
    leading-edge residual** at **obs 7593–7626 Å = def-rest 4525–4545 Å** (`ARC_MASKS_REST`).
    Confirmed telluric, not z=1.302 source, on 2026-05-18 by the same-observed-λ test (the spike sits
    at the *same obs λ* in the deflector aperture, the sky box, and a 4–8″ off-deflector spectrum). No
    full telluric correction is applied; the M10 OH/sky bands (below) catch sky-emission residuals.
  - **M10 sky-line audit:** bad-pixel catalog `BAD_PIXELS_REST` = 35 entries (26 CR residuals +
    9 OH airglow/sky bands), all verified non-source via the arc spectrum.
- **Wild bootstrap errors:** hybrid **Rademacher sign-flip × local-residual scaling** (75-pixel
  rolling-window residual rescale, clipped [0.2,5.0]); **N = 500 production** per (SPS × degree),
  N = 50 smoke (agree within 1 km/s). Errors = 16/84 percentiles (asymmetric). joblib parallel,
  **BLAS pinned to 1 thread/worker**.

**σ_e definition & wavelength window**

- **σ_e definition — analytic (Gültekin) and discrete (Cappellari) forms are the SAME quantity.**
  (Gültekin et al. 2009 eq. 1; Kormendy & Ho 2013 eq. 3; 2D-IFU form Cappellari et al. 2006 eq. 1.)
  - **Analytic / luminosity-weighted integral** (Gültekin 1D longslit, generalised to 2D — the 2πr
    is the axisymmetric area element dA = 2πr dr):

    $$\sigma_e^2 = \frac{\int_0^{R_e} (V^2 + \sigma^2)\, I(r)\, 2\pi r\, dr}{\int_0^{R_e} I(r)\, 2\pi r\, dr}.$$
  - **Discrete Cappellari (2006) eq. 1 spaxel/bin sum** — the *same* moment evaluated over spaxels
    (or Voronoi bins) n, with $F_n$ = the I-band flux (luminosity weight) of spaxel n; the 2πr factor
    is captured implicitly because the number of spaxels per annulus scales with area:

    $$\sigma_e^2 = \frac{\sum_{n:\, R_n \le R_e} F_n\,[(V_n - V_\mathrm{sys})^2 + \sigma_n^2]}{\sum_{n:\, R_n \le R_e} F_n}.$$

    This is the form SAURON IV / ATLAS3D III / MaNGA-DAP (Westfall+2019) use on Voronoi bins.
- **These two forms are equivalent and BOTH are already implemented and cross-checked in this project**
  (do not re-derive): the **analytic integral is the headline path** — the cumulative I-weighted
  single-ppxf fit on the R<R_e aperture spectrum (preserves LOSVD shape; what KH13/Greene+20 compute) —
  and the **discrete Cappellari/Gültekin spaxel-sum is the §7 / nb07c annular implementation**
  (`run_gultekin_mc`, F_j = Σ I over each unmasked annulus). The two agree at **<1σ**: integral
  (§6cum, narrow) 267.32 ± 24 vs discrete annular (§7, arc-filtered to R<R_safe=3R_e/4) 256.17 −13.0/+12.7.
  Full implementation + subtleties in `reference_gultekin_implementation.md`.
- **Why the integral path is the headline** (not the discrete sum): the discrete spaxel-sum is (i)
  per-spaxel-S/N-limited and (ii) sensitive to arc contamination in the outer spaxels — an *unfiltered*
  discrete sum over the resolved PowerBins out to ~3.4″ is dominated by the z=1.302 arc (outer bins
  carry σ ≳ 300–440 km/s, V swings of ±150–270; the σ_e blows up to ~350). The I-weighted aperture
  auto-suppresses the arc (bright deflector core dominates the weight) and avoids the S/N floor, so it
  implements the same Cappellari definition robustly. The discrete §7 sum is therefore arc-filtered to
  R<R_safe and used only as the architecture-independent cross-check / for the σ(R) profile.
- **Standard wavelength range (headline): rest 3800–5400 Å** (`wR3800_5400_arcmask`) — Ca H+K through
  Mg b / Fe5270 / Fe5335; 2161 good pixels (~2.6× the narrow window). The fit is constrained by the
  highest-S/N features — **Ca H, Ca K, G-band** (3800–4400) — and is **consistent across windows**.
  - Blue edge 3800 Å = XSL template floor (def-rest 3675 Å) + 5% velocity padding + retaining Ca H+K;
    red edge set by KCWI coverage. **Provenance is data-driven, NOT a specific paper.** Closest
    precedent: LEGA-C (Bezanson et al. 2018; van der Wel et al. 2021) at z~0.7.
    **DO NOT cite "Cappellari 2017 recommends extending blueward" — that was a hallucination, retracted.**
  - Cross-check window: **rest 4000–5400 Å** (`wR4000_5400_arcmask`, Hβ + Mg b, **no Ca H+K**) —
    orthogonal Lick-tradition feature set.
  - Narrow legacy window: obs 6500–7500 Å = rest 3879–4476 Å (Ca H+K + G); σ_e = 267.95 ± 30.10.

**A few systematics tests to describe** (full budget in Results §3.1):
I-shape (10 I-weight-map shapes, ±2.27), F200 mask weight (w∈{0,0.5,1}, ±6.65), fit-window
(3 windows, ±3.82), frame (vac/air, ±5.0), centering (5-center ±0.4″ sweep, ±4.0), reduction-pass
(NEW vs OLD cube, ±3.45). Plus the **M9 "DO NOT MASK" finding**: the +5–7% bump at def-rest
5193–5204 Å is the Mg b LOSVD wing (signal, not contamination); masking it drops σ_e by 7 km/s.

---

## §2.3 Lens Modeling

> **GAP G1 — not in this repo.** The companion-paper (Ferrami et al.) lens-modeling content
> (multiple lens models, PSF construction, joint shear+subhalo constraints, potential-correction
> triple-checks, BPL vs point-mass differentiating models) is **not present in this kinematics+
> photometry repo**. Source it from the lens-modeling repo
> (`~/Documents/AGEL/202509_DESJ0206_modeling/`, `20251112_DESJ0206_Pyauto_PRONTO/`) or the
> Ferrami draft. The only on-hand pointer: PROJECT_BRIEF.md notes the lens-model enclosed mass
> combines with (σ, M⋆, z) to give a 4D constraint on the central BH + bulge scaling relations.

System identity for the summary sentence: **AGEL J020613−011417** (DES J0206−01), AGEL DR2;
source z = 1.302 (spiral), deflector z = 0.67564.

---

## §2.4 End-to-end procedures (ordered walkthroughs)

*This section spells out the actual step order of each pipeline — what runs, in what sequence, with
which script — so the Methods section and a referee can reconstruct every headline number. Numeric
results are the M12 headline (`results/PAPER_VALUES.json`).*

### §2.4.1 Spectroscopic σ_e pipeline — step by step

Driver for the headline: `scripts/run_wide_sigma_e.py --cube new_clean_hei --n_bootstrap 500`
(calls `final_sigma_e.load_setup` → `extract_aperture_spectrum` → `bootstrap_ppxf` machinery).

1. **Reduction (upstream).** Keck II / KCWI, KCRM red arm, Red-Low grating, Medium slicer, 2×2
   binning, dithered PA 0/45/90°. The headline cube is the NEW **`_mtwdo_`** reduction
   (Master-Twilight+Dome hybrid flats, K. R. Gupta), `Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`.
   The pre-`_mtwdo_` cube is the OLD reduction (second-reduction cross-check). Cube axes
   `[λ, y, x]`; 0.30″/spaxel; seeing FWHM 1.27″ (`GUIDFWHM`); σ_inst ≈ 12.6 km/s; R ≈ 10000.
2. **Geometry + R_e.** HST-mean centre = `photutils.centroid_2dg` of the F140W+F200LP cores
   (offset 0.36″; F200LP pulled by the UV arc, so F140W is the clean bulge centroid), propagated
   through the KCWI WCS to a sub-pixel IFU centre. R_e = mean of the F140W (2.168″) + F200LP
   (2.441″) **masked** curves-of-growth = **2.305″** (`final_sigma_e.curve_of_growth`, r_max=6″).
3. **Spatial arc mask → IFU grid.** The F200LP `_mask.fits` (2512 HST px, arc-only) is reprojected
   onto the 0.30″ IFU grid with `scipy.ndimage.map_coordinates(order=0)` (nearest-neighbour) →
   `arc_spax_mask` (~38/184 spaxels inside R_e flagged, ~27% by I-weight).
4. **I-weighted aperture extraction.** Build the 6500–7500 Å white-light I-map; the R<R_e aperture
   spectrum = Σ_spaxel cube · w, w = I-weight with arc spaxels **hard-masked (mask_weight=0)** and
   re-normalised. Spectra are summed (not averaged) so the bootstrap preserves per-spaxel noise.
   (`extract_aperture_spectrum`, `final_sigma_e.py`.)
5. **Spectral masking (goodpixels).** (a) **no-Balmer** mask — forbidden lines only ([O II], [O III],
   [O I], [N II], [S II]); **Hδ, Hγ, Hβ kept** (absorption; M12/nb16 confirmed keep-unmasked).
   (b) **z=1.302 source-emission** masks `ARC_MASKS_REST` (Mg II, [O II], [Ne III], **He I 3819**) +
   the O₂ A-band telluric. (c) **`BAD_PIXELS_REST`** = 35 entries (26 CR residuals + 9 M10 OH/sky
   bands). All in deflector rest-frame Å.
   - **Outlier rejection / sigma-clipping decision:** ppxf's in-fit **`clean=False`** — `clean=True`
     rejects 0 pixels here (the KCWI noise array is overestimated vs the actual residual scatter), so
     it is inert. Instead, cosmic-ray/bad pixels are flagged by a **3σ local-MAD clip** (rolling
     75-pixel window, |residual|/local-MAD > 3σ) on the canonical residuals → 52 CR-like spikes,
     curated into the 26 CR entries of `BAD_PIXELS_REST`; the 9 M10 sky bands were flagged separately
     at >2.5× median sky-noise. **Caveats (carried as a stated sensitivity, not a budget term):** the
     3σ threshold is not strongly motivated vs alternatives, and cleaning shifts σ_e by only +1.4 km/s
     (below stat 1σ); the same pixels replicate the shift on the OLD cube (M6) → real CR residuals.
6. **Per-SPS frame-aware ppxf.** For each of FSPS / EMILES / XSL: match the galaxy frame to the
   SPS native frame (FSPS=vacuum, EMILES/XSL=air; scalar-median vac↔air via Ciddor 1996), run
   `ppxf` (Cappellari) with **moments=2**, **additive Legendre `degree` swept 15–29 (15 values)**,
   `mdegree=0`, over the wide window (rest 3800–5400 Å). Subtract the per-SPS V_sys before pooling.
7. **Wild-bootstrap errors.** Per (SPS × degree), **N=500** hybrid wild bootstrap
   (`bootstrap_ppxf.compute_local_residual_scaling`): residuals r = galaxy − bestfit; a rolling
   **75-pixel** window estimates the local (wavelength-dependent) noise scale s(λ) (clipped to
   **[0.2, 5.0]**); each iteration applies a **Rademacher** ±1 sign-flip to the locally-rescaled
   residual, adds it back to the best fit, and re-fits → one σ draw. Seeds are deterministic
   (`BOOT_SEED=42 + 50000 + 10000·W_IDX + 100·sps_idx + degree_i`) so caches bit-reproduce.
   joblib `loky`, **BLAS pinned to 1 thread/worker**; N=50 smoke before N=500 production.
8. **Pool + budget.** Concatenate all σ samples (3 SPS × 15 degrees × 500) → 16/50/84 percentiles
   = central + asymmetric stat. The pooled width already marginalises SPS + polynomial degree. Add
   the §2.4.3 systematic components in quadrature → **σ_e = 269.62 −13.45/+13.10 (sym ±13.27) km/s**.

### §2.4.2 Photometry / M⋆ pipeline — step by step

Drivers: `scripts/photometry_systematics.py` (masking + photometry), `bagpipes_sersic_refit.py`
(SED), orchestrated/displayed in notebook 12. Recipe = **F200LP locates, IR extends, fill-in recovers**.

1. **Locate the arc on F200LP** (best source contrast): a 2D Sérsic deflector model is subtracted
   and positive residual > 3σ_sky flags arc pixels (`arc_mask_verification.py`). This **objective
   Sérsic-residual mask reproduces the expert hand mask to 0.016 mag**, R_e to 3% (k=3 is the clean
   S/N saturation point from a SNR∈{2..20}×k∈{2..8} sweep).
2. **Reproject the F200LP footprint to every band** (WCS + `map_coordinates`) — reproduces the
   expert JWST aperture mags to 0.01–0.02 mag. Do **not** fit an independent per-band single Sérsic
   (under-fits the bright IR galaxy, over-masks 0.2–0.4 mag).
3. **IR-extend.** In the deeper/sharper JWST bands, region-grow the arc into 2-component-Sérsic
   residual source pixels **contiguous** with the arc seed (cannot grab the core pedestal or field
   companions); F150W2 mask grows ~6× vs HST.
4. **Deflector model = 2-component (bulge+disk) Sérsic** (cuts JWST galaxy-body RMS 3.3–3.5× vs
   single Sérsic). PSF-convolved fill quantified at ≤0.004 mag (`psf_fill_model.py`) — negligible.
5. **Two photometry estimators per band/mask:** **raw** `photutils` aperture (masked → discards
   under-arc deflector light → biased low, mask-size-dependent) vs **Sérsic fill-in** (replace
   masked pixels with the model → recovers it → mask-definition-independent, per-band vs global
   agree 0.01 dex). Fill-in correction +0.18–0.96 mag with the large IR-extended masks.
6. **AB magnitudes** from per-image header zeropoints (never hardcoded): HST `PHOTFLAM/PHOTPLAM`,
   JWST `PIXAR_SR`. Four bands: F200LP, F140W (HST WFC3), F150W2, F322W2 (JWST NIRCam).
7. **Bagpipes SED** (`02_streamlined_Bagpipes_SED`): exponential-τ SFH, Calzetti dust (A_V 0–2),
   free Z (0–2.5 Z☉), z prior (0.674,0.676); **Nautilus** sampler n_live=400, 500-sample posterior.
   A **10% fractional flux floor** per band covers photon noise + the drizzle pixel-correlation
   under-estimate (empirical noise ≈1.3–2× CCD-equation); **20%** is run as a conservative variant.
8. **M⋆ budget** = 8 fits {per-band, global} × {raw, filled} × {10%, 20%} → headline (empirical,
   raw-central, one-sided +sys, user choice): **log M⋆ = 11.16 ± 0.08 (stat) +0.31 (sys) at 10%**
   [11.04 ± 0.14 at 20%], fill-in reach 11.46 (the asymmetric form used in the figures is
   +0.32/−0.08, the stat⊕sys quadrature); explicit masking-approach systematic ±0.16 dex
   (`Mstar_masking_systematic.npz`).

### §2.4.3 Systematic-audit procedures (how each budget term was measured)

Each σ_e systematic is a **dedicated sweep** re-using the §2.4.1 machinery, varying ONE axis and
quoting **peak-to-peak/2** (or half-Δ for two-point axes), all on the headline cube + masks:

- **I-shape (±2.27):** 10 alternative I-weight maps (`run_isource_shape_sweep.py`) × 3 SPS × 15 deg
  × N=250; peak-to-peak/2 of the 10 central σ_e (266.83–271.37).
- **F200 mask-weight (±6.65):** down-weight the arc spaxels at w∈{0,0.5,1} × 3 SPS × N=500;
  (w00 269.69 − w100 256.39)/2. *(This is the mask-**strength** axis. The mask-**definition** axis —
  expert/sersic/perband/global reprojected to the IFU grid, `run_sigma_e_mask_systematic.py` — gives
  ±5.85, which overlaps this; larger-of-two is kept, not added.)*
- **Fit-window (±3.82):** 3 windows (wR3800_5400 / wR4000_5400 / w6500_7500) × 3 SPS × N=500;
  peak-to-peak/2 (269.66 / 268.58 / 276.22).
- **Frame (±5.0):** structural — vac/air per-SPS native-frame choice.
- **Centering (±4.0):** 5 perturbed HST-mean centres (±0.4″ sweep; `NOTES_centering_investigation`).
- **Reduction-pass (±3.45):** half-Δ between the NEW and OLD cleaned cubes (269.62 vs 262.72); 2
  reductions only — refine if a 3rd lands.
- **R_e-source / D7 (±6.13):** the 4 R_e estimators (mean/F140W/F200LP/CaHK+G, 2.168–2.902″;
  `run_sigma_e_Re_systematic_wide.py`); full peak-to-peak/2.

Quote the budget table (§3.1) as the result; the iterative history that *led* to this masking +
component set is the **internal audit trail (M8→M12)** at the end of this document — not for the
manuscript Methods section.

### §2.4.4 Reproducibility / infrastructure

- **Architecturally-independent estimators** (must not share intermediate state): **§6cum** cumulative
  I-weighted single-ppxf aperture (the analytic Gültekin path = headline) · **§7** discrete
  Cappellari (2006) annular spaxel-sum (`run_gultekin_mc`, arc-filtered to R<R_safe) · **nb07e**
  arc-spectrum subtraction. The three agree at <1σ (267.32 integral vs 256.17 discrete, narrow).
- **Caching + seed discipline.** Every expensive bootstrap writes a `.npz` keyed by its parameters
  (`results/<sweep>/<params>_N{N}.npz`); a cached result is **never** re-run. Seeds are deterministic
  (§2.4.1 step 7) so re-runs bit-reproduce; **N=50 smoke before N=500 production**; joblib with
  BLAS=1/worker to avoid oversubscription.
- **Deterministic values registry (single source of truth).** `scripts/paper_values.py` recomputes
  **every** headline number from the result caches (with per-entry provenance + formula) → emits
  `results/PAPER_VALUES.json`. The ApJL figures **load** that JSON (no hard-coded literals); the
  summary docs are validated by `paper_values.py --check <files>`, an anchored, precision-aware
  **drift linter** (parses only the canonical `σ_e = N ± M` / `log M⋆ = N` statements; dated
  snapshots carry a `pv-skip-file` marker). Regenerate + verify:
  `python scripts/paper_values.py && python scripts/paper_values.py --check *.md`.

---

# §3 Results

## Figure 2 (left panel) — the integrated spectrum + SED

- **Left:** KCWI/KCRM Red-Low integrated, I(r)-convolved spectrum of the deflector elliptical, with
  the ppxf best fit overlaid (and inline residuals — Rodrigo's standard nested-gridspec style).
  Loads from the wide-window fit + bootstrap pool.
- **Right:** HST+JWST four-band SED + Bagpipes posterior model. Loads `results/bagpipes_sed_results.npz`.
- Figure file: `AGEL0206_sigma_e_SED_final_wide.pdf` (per outline) / `results/AGEL0206_spectra_SED_fit.pdf`.

## §3.1 Velocity Dispersion Measurement

- **Headline: σ_e(<R_e) = 269.62 ± 13.27 km/s** (asym −13.45/+13.10), via the Gültekin (2009)
  luminosity-weighted formalism evaluated with the Cappellari ppxf implementation on the
  I(r)-weighted R<R_e aperture spectrum, plus the wild-bootstrap error pool.
- **Reproduce:** `scripts/run_wide_sigma_e.py --cube new_clean_hei --n_bootstrap 500`.
  Caches: `results/run_wide_sigma_e/new_clean_hei/wR3800_5400_arcmask_{fsps,emiles,xsl}_T*_N500.npz`.

**R_e measurement and checks** (for the §3.1 paragraph): R_e = 2.305″ = 16.23 kpc, the F140W (2.168″)
+ F200LP (2.441″) masked-CoG mean; supersedes IFU-only 2.61″; masked CoG removes arc/companion bias.

**Systematic error budget (post-M12, 2026-06-08; computed by `scripts/paper_values.py`):**

| Component | ± km/s | Why |
|---|---|---|
| stat (N=500) | 4.6 (−5.16/+4.13) | wild-bootstrap pool over 3 SPS × 15 degrees; already marginalizes SPS (between-SPS ±2.04) + degree spread |
| I-shape (10 shapes) | 2.27 | peak-to-peak/2 of 10 I-weight-map shapes (266.83–271.37) |
| F200 mask (3 weights) | 6.65 | (w00 269.69 − w100 256.39)/2; cleaned cube has more signal pixels → more arc-dilution leverage. (Arc-mask-*definition* cross-check ±5.85 overlaps this, not added — M12/nb13) |
| frame (vac/air) | 5.0 | structural per-SPS native-frame choice |
| centering (HST WCS) | 4.0 | 5-center ±0.4″ sweep; reduction-independent |
| fit-window (3 windows) | 3.82 | wR3800_5400 269.66 / wR4000_5400 268.58 / w6500_7500 276.22; **dropped from carried ±15** on the cleaned cube |
| reduction-pass | 3.45 | (269.62 − 262.72)/2, NEW vs OLD cleaned cube; **only 2 reductions — refine if a 3rd lands** |
| **R_e-source (D7 wide)** | **6.13** | **NEW M12 (2026-06-08, nb15):** peak-to-peak/2 over 4 R_e estimators (mean 2.305″/F140W 2.168″/F200LP 2.441″/CaHK+G 2.902″ → 269.62/267.44/272.44/279.69). User chose full spread; CaHK+G's +10 is the rising σ(R) gradient. Light-family-only = ±2.50 |
| **TOTAL (sym)** | **13.27** | quadrature |
| **TOTAL (asym)** | **−13.45 / +13.10** | preserves stat skew |

- **No separate SPS line** is folded in: the SPS-library spread **collapsed from ~26 km/s at the
  narrow window** (FSPS 253.0 < EMILES 267.5 < XSL 279.7) **to ~4 km/s at the wide window** (FSPS
  271.3 / EMILES 267.2 / XSL 270.4; `PAPER_VALUES.json:per_sps_p50`) — a key argument for the wide
  window; it now sits inside the pooled stat width (between-SPS ±2.04).
- Quadrature is deliberately conservative (components not strictly independent).
- **R_e-source systematic IS now folded in (±6.13, M12 2026-06-08):** re-measured at the wide window
  (nb15) and included in the budget above — the narrow-window ±16.9 is superseded. The full
  4-estimator spread was adopted (user decision); the CaHK+G 2.90″ estimator drives it (+10 km/s =
  rising σ(R), see below), while the light-R_e family alone gives ±2.50.

**Systematics paragraph order** (per outline): first R_e and I(r) (I-shape ±2.27, R_e-source flagged),
then the spectrum (mask ±6.65, frame ±5.0, window ±3.82, reduction ±3.45, centering ±4.0).

**Cross-checks** (narrow window, architecturally independent, superseded for headline but valid):

- §6cum cumulative I-weighted aperture ppxf: 267.32 ± 24; σ_e(<R_e/2) ≈ 225.78; σ_e(<R_e/8) ≈ 209.18.
- §7 discrete Gültekin annular sum (equal-N inside R_safe=3R_e/4): 256.17 −13.0/+12.7.
- §7b flat-σ extrapolation: 274.37 −16.2/+17.4.
- nb07e arc-spectrum subtraction matches §6cum bit-identically (<1 km/s).
- Narrow Ca H+K + G window (nb09d): 267.95 ± 30.10 — headline within ~2 km/s.

**Headline evolution** (mask-audit arc): 268.98 (clean) → 271.87 (+He I, M8) → 269.62 (+M10 sky audit)
→ ±11.77 (M11 re-derivation) → **±13.27 (M12, +R_e-source ±6.13, 2026-06-08)**. He I pushes σ up ~2.9,
the 9 sky bands pull it down ~2.2 (near-cancel); central value 269.62 stable across M10→M12 (only the error grew).

**Compare with Greene (2016/2020) and KH13 choices:** σ_e(<R_e) here uses the Cappellari et al. 2006
eq. 1 luminosity-weighted single-aperture definition, exactly as adopted by Kormendy & Ho 2013
(M•–σ eq. 3) and Greene et al. 2020 ARA&A (fig. 5). We dropped R<R_e/8 (JFK95 σ_c) because seeing
under-resolves it. For a pure elliptical, galaxy R_e ≈ bulge R_e, so KH13's "bulge R_e" equals the
measured total R_e — no bulge/disk decomposition needed.

**σ_e physical-meaning discussion:** the deflector is a dispersion-supported **elliptical bulge**;
σ(r) should *decrease* gently with radius — any rising σ(r) in raw output is a continuum-dilution
artifact from the z=1.302 arc, not a physical gradient. Resolved kinematics (nb11): 17 central
spaxels at S/N≥5, σ ∈ [144, 251] km/s (median 201), declining with R, no coherent rotation above noise.

## §3.2 Stellar Mass Estimate

- **log(M⋆/M☉) = 11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors [**11.04 ± 0.14 +0.32** at 20%]
  — NEW headline from the principled F200LP-located + IR-extended arc masking with **raw** aperture
  photometry (`scripts/photometry_systematics.py`; `results/photometry_systematics_Mstar.npz`;
  fig `results/figures/Mstar_budget.png`). **Empirical / conservative:** raw masked photometry
  discards the deflector light under the (now IR-extended) arc, so the central is biased LOW and
  the systematic is **one-sided UP**.
- **The systematic is intentionally uneven (could be higher):** the Sérsic fill-in — which restores
  the deflector light under the masked arc and is **mask-definition-independent** (per-band vs global
  agree to 0.01 dex) — gives log(M⋆/M☉) = **11.46 (10%) / 11.36 (20%)**, i.e. the +0.31 dex upper
  reach. The mask-definition term is negligible (±0.003–0.034 dex).
- **Flux-error choice matters at ~0.10 dex:** 10%→20% lowers the central (11.16→11.04, looser
  likelihood lets the prior pull mass down) and widens stat (±0.08→±0.14). Both quoted.
- **Explicit masking-approach systematic on M⋆ = ±0.16 dex** (10%: ±0.160, 20%: ±0.162) —
  peak-to-peak/2 across all masking approaches {expert-aperture, per-band/global × raw/filled}
  (`results/Mstar_masking_systematic.npz`, fig `Mstar_masking_systematic.png`). Decomposition:
  **under-arc light (raw↔filled) ±0.15** (dominant; the adopted one-sided +0.31), **mask-definition
  (per-band↔global) ±0.004** (negligible — filled is mask-def-independent), **mask-extent
  (expert↔IR-extended) 0.18**. This ±0.16 is the symmetric statement of the same effect quoted
  one-sided (+0.31) in the headline.
- **Superseded values (cross-checks, smaller masks):** older expert-aperture log(M⋆/M☉)=11.33
  +0.07/−0.09 (sits mid-range); 2D-Sérsic-total refit 11.40 +0.11/−0.15. Consistent with the new
  raw↔filled bracket [11.16, 11.46].
- **"Typical elliptical at z~0.7":** M⋆ ≈ 1.4 × 10¹¹ M☉ (raw headline) to ≈ 2.9 × 10¹¹ (fill-in upper) at z = 0.6756, passive/quiescent SFH
  (mass-weighted age ~5 Gyr, low sSFR). Place on the z≈0.5–1.0 mass-size relation (van der Wel et al.
  2014; Mowla et al. 2019) and the strong-lens-deflector analog sample (Sonnenfeld et al. 2015,
  SL2S V — best z-match). [notebook 10 verified citations]
- **X-ray / quiescence:** **GAP G2 — no X-ray data in repo.** Source externally if you want the
  "truly quiescent" claim; otherwise frame from the passive SED only.

## §3.3 Mass Estimates from Lens Modeling

> **GAP G1 — see §2.3.** The range of lens-model results, the joint-posterior mass + uncertainty, the
> BPL vs point-mass differentiating models, M(<r_break), the pedagogical-figure intuition, the small
> table reproduced from Ferrami et al., the comparison to local galaxies, and the "shape of the mass
> profile matters" argument are all **in the lens-modeling repo / Ferrami draft, not here.**
> The only in-repo pointer is the 4D-constraint sentence in PROJECT_BRIEF.md. Pull these numbers from
> the companion paper before drafting §3.3.

---

# Key citations (verified in notebook 10)

- **ppxf:** Cappellari & Emsellem (2004); Cappellari (2017, MNRAS 466, 798); Cappellari (2023).
- **σ_e aperture definition:** Cappellari et al. (2006, MNRAS 366, 1126) eq. 1;
  **Gültekin et al. (2009)** eq. 1; **Kormendy & Ho (2013, ARA&A 51, 511)** eq. 3 (scatter 0.29 dex);
  **Greene, Strader & Ho (2020, ARA&A 58, 257)** fig. 5 **[⚠ AUDIT FLAG — see below]**; Saglia et al. (2016).
- **z=0 IFS calibration:** Cappellari et al. (2013, ATLAS3D XV & XX); Krajnović et al. (2013, XVII);
  Zhu et al. (2024, MaNGA DynPop III).
- **Mass–size / analogs:** van der Wel et al. (2014); Mowla et al. (2019); Ge et al. (2021);
  Sonnenfeld et al. (2015, SL2S V) **[⚠ AUDIT FLAG — see below]**.
- **Vacuum↔air:** Ciddor (1996). **Wide-window precedent:** Bezanson et al. (2018); van der Wel et al. (2021).
- **Source-z catalog:** AGEL DR2 (Carleton et al., in prep). **Cluster:** Hilton et al. (2021, ACT DR5).

# Citation audit + referee challenges (2026-06-08)

**Citation audit** (vs the user's Zotero library; "not in library" ≠ fabricated — most are real
foundational papers simply not yet added to Zotero, flagged so the .bib can be completed):

- **Verified in-library, metadata exact:** Cappellari (2017, MNRAS 466, 798) ✓; Kormendy & Ho
  (2013, ARA&A 51, 511) ✓; Greene, Strader & Ho (2020, ARA&A 58, 257) ✓ *(metadata only — see flag)*;
  Hilton et al. (2021, ApJS 253, 3 = ACT SZ-cluster catalog) ✓; Cappellari (2023, MNRAS 526, 3273) ✓.
  Also in-library and relevant: **McConnell & Ma (2013)** — a massive-BH M•–σ / M•–M_bulge anchor.
- **⚠ FLAG 1 — Greene, Strader & Ho (2020, ARA&A 58, 257):** bibliographic metadata is CORRECT, but
  this review is titled **"Intermediate-Mass Black Holes."** It does provide M•–σ / M•–M★ relations
  (legitimately citable as a *comparison* relation, as Fig. 4 uses it), but for a ~10⁹ M☉ SMBH the
  *primary* anchors should be KH13 / McConnell & Ma 2013 / Saglia+2016, with Greene+20 as one of the
  overplotted comparison relations. **ACTION: verify the "fig. 5" number** (the scaling-relation
  figure in Greene+20 may not be fig. 5) and ensure the text frames it as a comparison, not the
  primary massive-BH relation. *(Not a fabrication — a framing/figure-number check.)*
- **⚠ FLAG 2 — Sonnenfeld "SL2S V" (2015):** the library contains SL2S **IV** (Sonnenfeld et al.
  **2013**, ApJ 777, 98), not paper **V** (2015, ApJ 800, 94). **ACTION: decide which you mean** and
  either correct to "IV (2013)" or add the 2015 V paper to Zotero. Do not let the 2013 record stand
  in for a 2015 citation.
- **⚠ FLAG 3 — Bezanson et al. (2018):** no Bezanson-first-author 2018 item in the library. If you
  mean the LEGA-C velocity-dispersion paper (Bezanson et al. 2018, ApJ 858, 60), **add it** / verify.
- **To add to Zotero before the .bib is complete** (real papers, just absent; metadata NOT
  independently verified here — confirm on add): Cappellari & Emsellem (2004), Cappellari et al.
  (2006, SAURON IV), Gültekin et al. (2009), Saglia et al. (2016), Cappellari et al. (2013, ATLAS3D),
  Krajnović et al. (2013), Zhu et al. (2024, MaNGA DynPop III), van der Wel et al. (2014),
  Mowla et al. (2019), Ge et al. (2021), Ciddor (1996), van der Wel et al. (2021, LEGA-C DR3).

**Referee challenges to the current doc** (anticipate / address in drafting):

1. **R_e is not converged (A3c, nb14).** The headline R_e = 2.305″ is a *raw* curve-of-growth at
   r_max=6″; the CoG never flattens (no sky pedestal subtraction), so R_e climbs with r_max and sits
   at the top of a 2.1–2.5″ method family (sky-sub 2.14″, Sérsic 2.18″). **Strongest soft spot.** A
   referee will ask for a Sérsic / sky-subtracted / PSF-deconvolved R_e. Currently flagged, not fixed.
   Knock-on: R_e sets both the σ_e aperture and (weakly) M★.
2. **R_e-source budget term uses the FULL 4-estimator spread (±6.13).** Including the CaHK+G 2.90″
   estimator — which is larger than the light-R_e family and whose +10 km/s is the rising-σ(R)
   gradient — arguably conflates aperture-uncertainty with σ(R)-physics. Defensible (conservative),
   but a referee could argue for the light-family ±2.50 (→ ±12.0 total). Documented both ways.
3. **Rising outer σ(r) attributed to arc dilution.** The claim "σ(r) should decrease; any rise is the
   z=1.302 arc" is the right physical call, but it is *also* what lets us drop the high-σ outer bins;
   make the arc-contamination evidence explicit (EW / continuum-dilution diagnostic, nb07e) so it does
   not read as motivated reasoning.
4. **He I 3819 source-emission mask** rests on a 3-pixel residual cluster — defensible (consistent
   across reductions, matches z=1.302), but small; a referee may probe it. The M8/M9/M10 audit trail
   covers this.
5. **σ_inst = 12.6 km/s vs σ_e ≈ 270:** well-resolved, not a challenge — but state the LSF treatment.

# Flagged open items (carry into drafting)

1. ~~R_e-source systematic not folded in~~ **RESOLVED 2026-06-08 (M12, nb15):** re-measured at the
   wide window (±6.13) and folded into the budget → headline ±11.77 → **±13.27**.
2. Reduction-pass ±3.45 rests on only 2 reductions — refine if a 3rd lands.
3. ~~Hδ may still need targeted masking~~ **RESOLVED 2026-06-08 (M12, nb16):** decision = **keep
   unmasked** (Hδ is well-fit, not a contaminant; masking it would discard LOSVD information).
   TODO in `bootstrap_ppxf.py` closed.
3b. **NEW (A3c, nb14):** raw CoG R_e is r_max-non-convergent; headline 2.305″ is the top of a
   2.1–2.5″ family (**R_e method systematic ±0.08″**), not folded into σ_e/M★ — flagged for a decision.
4. Asymmetric-error tables: older 4-corner/M10 tables print −17.92/+17.65 (pre-M11), M11 was
   −11.98/+11.57 (±11.77); **current (M12) is −13.45/+13.10 (±13.27). Quote the M12 value**
   (canonical: `results/PAPER_VALUES.json`).
5. No X-ray / quiescent data in repo (G2). 6. No lens-model content in repo (G1).
7. "We do not account for peculiar velocities" needs to be written in (G3 — true, just add it).
8. Source-z (1.302) error bar + its cache are a `[TBD]` placeholder.
9. K409 (Aug-30) PI and per-night DIMM seeing still TBD for acknowledgements.

---

# Internal appendix — σ_e measurement audit trail (M8→M12) — NOT for the manuscript

*This is the iterative history of how the final §2.4.3 masking + systematic-component set was reached.
The Methods section should present the **final** procedure (§2.4) and the budget (§3.1), not this
trail. Kept here as internal provenance / referee-rebuttal backup. Each step is also a `TESTS` row
with its own cache.*

- **M8 — He I 3819 (added as a source-emission mask).** A 3-pixel positive-residual cluster at
  def-rest 5244–5248 Å = source-rest 3818–3820 Å matches He I 3819.6 from the z=1.302 source;
  consistent across both reductions → astrophysical, not a CR. Masked. σ_e (NEW) 268.98 → 271.87.
- **M9 — 5193–5204 Å bump (deliberately NOT masked).** A +5–7% feature; checked against the arc
  spectrum (no source counterpart) and the noise spectrum (no sky counterpart) → it is the **Mg b
  LOSVD wing** (signal). Smoke-masking it drops σ_e by 7 km/s → **DO NOT mask.** (Same "it's signal,
  not contamination" logic later applied to Hδ in M12.)
- **M10 — full sky-line audit.** Flag every band with cube `noise_sky` > 2.5× median across the fit
  window, and cross-check each against the arc spectrum to confirm no source counterpart → 9 OH/sky
  bands added to `BAD_PIXELS_REST` (26 → 35 entries). σ_e (NEW) 271.87 → 269.62.
- **M11 — cube-matched systematic re-derivation.** Re-ran the I-shape / mask-weight / fit-window
  sweeps on the NEW cube + M10 masks (replacing carried OLD-cube values). The fit-window term
  collapsed ±15 → ±3.82 on the cleaned cube; total sys ±17.16 → ±10.81 (total ±11.77).
- **M12 (2026-06-08) — R_e-source fold-in + two cross-checks.** Folded the wide-window R_e-source
  systematic (±6.13) into the budget → sys ±12.43, **total ±13.27**. Added the arc-mask-**definition**
  cross-check (±5.85; overlaps the F200-mask term, not double-counted) and resolved the Hδ open item
  (**keep unmasked** — well-fit, masking it discards LOSVD information). Central σ_e 269.62 unchanged
  across M10→M12; only the error grew.

**Net effect of the trail:** central σ_e settled at 269.62 km/s by M10 and has not moved since; the
error bar evolved ±17.78 (pre-M11) → ±11.77 (M11) → ±13.27 (M12). The mask set (He I + 35-entry
sky/CR + no-Balmer) and the 7-component budget that the final §2.4 procedure uses are the M12 result.
