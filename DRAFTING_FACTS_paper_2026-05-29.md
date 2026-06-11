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
The **`<!-- PV:auto:… -->` blocks** (the headline numbers below + the §3.1 budget table) are
**generated from `results/PAPER_VALUES.json`** — do not hand-edit inside the markers; refresh with
`python scripts/paper_values.py --render DRAFTING_FACTS_paper_2026-05-29.md`.
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

<!-- PV:auto:headline -->
**Headline numbers** — *generated from `results/PAPER_VALUES.json` by `scripts/paper_values.py --render`; do not hand-edit inside the markers.*

- σ_e(<R_e) = **267.31 ± 12.79 km/s** (asym −12.98 / +12.61); often rounded **267 ± 13 km/s**. Measured at the best-mask R_e=2.097″ aperture; the R_e-source systematic (±5.13, best-mask CoG family) is folded into the budget (§3.1).
- log(M⋆/M☉) = **11.47 +0.09/−0.14 (stat, 10%)** — aperture-corrected total at the matched 2 R_e = 4.19″ elliptical aperture (empirical aperture + single-Sérsic model wings; GAMA/Taylor+2011 fluxscale). Five estimators: raw **11.18** (empirical) · raw+apcorr **11.35** · filled **11.36** · **total 11.47 (headline)** · Sérsic-total **11.40**. Sérsic-only systematic ±0.19 dex; masking-approach ±0.09 dex. See §2.1.1b / §3.2.
- R_e = **2.097″ = 14.76 kpc** (best-mask F140W+F200LP CoG; method systematic ±0.10″). **Supersedes 2.305″** (expert mask).
- z_deflector = **0.67564**; z_source = **1.302**.
<!-- /PV:auto:headline -->

**Cosmology:** FlatLambdaCDM, H₀ = 70 km/s/Mpc, Ω_m = 0.3; 1″ ≈ 7.04 kpc (7043.5 pc/arcsec),
D_A = 1452 Mpc at z = 0.67564. [reference_hst_jwst_data_properties.md]

> **Staleness warning:** older docs (`METHODS_AND_SYSTEMATICS.md` 2026-05-18 body; pre-2026-06-11
> handoffs) print superseded σ_e numbers (254.85 / 268.98 / 269.62, totals ±18/±13.27). The current
> headline is **267.31 ± 12.79** (post "best mask throughout", 2026-06-11, R_e=2.097″). Always
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
  red deflector); (ii) a **Sérsic-residual** mask (subtract a 2D deflector model, flag pixels whose
  positive residual exceeds **k·σ_sky**, where **k is the detection threshold in units of the sky
  noise σ_sky** and σ_sky is the robust background scatter; we adopt **k = 3**). On F200LP the
  Sérsic-residual mask reproduces the expert photometry to **0.016 mag**, R_e to 3%.
- **Transfer to other bands by reprojecting the F200LP footprint** (WCS + `map_coordinates`).
  Reproduces the expert JWST aperture mags to **0.01–0.02 mag**.
- **WCS registration (NEW 2026-06-10).** A raw WCS reprojection leaves a **0.09–0.18″** residual
  offset between the F200LP arc and the JWST/F140W frames (HST UVIS↔IR and HST↔JWST astrometry; the
  JWST i2d GWCS differs from the FITS-header WCS astropy reads). On the thin arc that mis-covers one
  edge and lets a companion slip past the mask. We **cross-correlation-register** the reprojected
  F200LP arc (and the F140W color layer) to each band's own image near the deflector
  (`photometry_systematics.register_shift`), which drives the residual offset to **0.00″** (verified
  F140W 0.179→0, F150W2 0.154→0, F322W2 0.089→0). Field-galaxy contamination check after
  registration: **no unmasked source (>8σ) falls inside the photometry aperture** in any band.
  *Do not* anchor on the F200LP deflector centroid instead — it is arc-biased (the faint blue
  deflector centroid is pulled ~0.26″ by the arc), so the arc cross-correlation is the right fiducial.
- **Deflector model = single-component Sérsic** (per user 2026-06-10). The deflector is an
  elliptical, so a single Sérsic is the physically-motivated light profile — we do **not** fit a
  bulge+disk decomposition (an earlier 2-component model was only a numerical kludge to make the
  residual-based mask growth behave; see the mask-evolution note below). The single Sérsic is fit
  per band with masked pixels excluded, capturing **66–87%** of the box flux (F140W is the weakest,
  with a localized central residual from the absent bulge/disk split, but that residual sits inside
  the arc-free core and does not enter the mask).
  - *Implementation gotcha (documented so it isn't reintroduced):* `astropy.modeling.Sersic2D`'s
    `amplitude` parameter is **I(r_eff), not the central peak**. Since I(0)/I(r_eff)=exp(b_n)≈2150
    at n=4, initialising `amplitude` at the data peak over-predicts the centre ~2000× and LevMar
    collapses the amplitude to ≈0 (the whole galaxy then reads as positive residual → catastrophic
    over-masking). Initialise `amplitude ≈ peak·exp(−b_n)` (~peak·0.05 for n≈3) and **bound it to
    (peak·1e-3, peak·1.5)**, exactly as `sersic_total_photometry.fit_sersic2d`.
- **IR extends the source footprint — grown band-specifically (HST color gate / deep-JWST
  morphology guard):** the deeper/sharper JWST bands see the arc's fuller extent, so we region-grow
  the F200LP arc seed into the IR. The "is this pixel source, not deflector body?" test is gated two
  ways depending on the band's depth (a candidate must also be a >k·σ_sky residual, within R<5″, and
  spatially contiguous with the arc seed — region growing via `scipy.ndimage.label`):
  1. **HST bands (F200LP, F140W) — COLOR GATE:** the pixel must be **bluer than the deflector core**
     by ≥ **ΔC = 0.5 mag** in **F200LP−F140W** (`DCOLOR`), at **k = 3** (`K_EXT`).
  2. **Deep JWST bands (F150W2, F322W2) — MORPHOLOGY GUARD:** the pixel must be **source-dominated**,
     `data − sky − model > MORPH_FRAC·model` (excess beats the single-Sérsic deflector model → it is
     OUTSIDE the Sérsic-dominated region), at the lower **k = 2** (`K_EXT_JWST`). **No HST color is
     required**, so the deep JWST data captures the **diffuse** arc emission HST is too shallow to
     color-confirm. Because the single Sérsic already accounts for 66–87% of the deflector flux,
     `resid > model` rejects the deflector body (model-dominated) while admitting the diffuse source
     (excess-dominated) — verified: raw mag stable (no deflector loss), 0 unmasked >8σ sources in the
     aperture, the deflector core <0.4″ never masked. F150W2 IR-extension grows 1049→5288 px (the
     diffuse arc), F322W2 415→1230 px.
  - **Why the color gate is needed (and is the right physics).** A *single* Sérsic, being the
    correct elliptical model, still under-fits the bright extended IR galaxy body slightly, leaving
    a low-level **positive** residual over the *red* galaxy. Residual-only growth (criterion 1 alone,
    the old 2-component recipe) therefore sweeps red galaxy-body light — and even neighbouring
    companions — into the "arc" mask (verified: F150W2 mask balloons and raw mags fall 1–2 mag).
    The z=1.302 source is a **blue, star-forming** lensed galaxy; the deflector is a **red, passive**
    elliptical. Requiring the grown pixels to be **bluer than the deflector core** (criterion 2)
    exploits this intrinsic 4000-Å-break contrast to keep the blue source extension while rejecting
    the red body. The mask thus becomes a property of the **source color**, not of the (imperfect)
    deflector residual. This is the same physical contrast used by the validated F200LP color
    selector (i) above; we are simply applying it per-pixel during the IR growth.
  - **Color construction.** Per target-band grid, reproject the F200LP and F140W surface brightness
    (`scripts.arc_mask_verification.reproject_intensity`, WCS + bilinear `map_coordinates`) and form
    `−2.5·log10(SB_F200LP/SB_F140W)` (arc → more negative). Constant ZP offsets **cancel** because
    we threshold on the contrast vs the deflector-core color (median in a 0.4–0.8″ annulus), not on
    an absolute color.
  - **Caveats.** (a) *(RESOLVED 2026-06-10)* the HST color gate is conservative where HST is shallow
    and was missing faint diffuse JWST arc emission (user catch via the deep-stretch
    `_diffuse_mask_check`); the **deep-JWST morphology guard** above now captures it using JWST's own
    sensitivity + the source-dominance test, without over-masking the deflector (raw mag stable). (b)
    ΔC = 0.5 mag and MORPH_FRAC = 1.0 are tunables; the resulting masks reproduce the expert hand-mask
    M⋆ (11.33), our calibration that they neither over- nor under-mask. (c) The morphology guard
    assumes the single Sérsic captures the deflector majority — true here (66–87%); if a future target
    has a worse single-Sérsic fit, re-check that `resid>model` still protects the body.
- **raw vs fill-in:** raw photutils (masked) discards the deflector light under the arc; the
  single-Sérsic fill-in (replace masked pixels with the model) recovers it. Under the color-gated
  masks the fill-in correction is **−0.19 to −0.27 mag** per band (smaller than the old 2-component
  +0.18 to +0.96 mag, because the cleaner masks remove less deflector light). The **headline is the
  empirical raw value with the best (per-band, color-gated) mask**; the fill-in is carried as a
  **one-sided +sys** (deflector light the mask removes). **PSF effect on the fill** stays negligible
  (`scripts/psf_fill_model.py`: ≤0.004 mag → ΔM⋆ ≪0.01 dex; arc is outside the PSF core).
- **per-band vs global mask** (union of all per-band masks): masking systematic **±0.086 dex**
  (10%); the raw↔filled (under-arc) term dominates (±0.057), the per-band↔global mask-definition
  term ±0.028. See §3.2.

- **Mask-evolution note (why we settled on the single-Sérsic + color-gated mask).** Three attempts,
  documented in `results/figures/mask_attempts_comparison.png` + notebook 12 §8:
  - **(A) 2-component (bulge+disk) Sérsic + residual-only IR growth** — over-masks: residual-only
    growth grabs red galaxy body + a companion (F150W2 15768 px, F322W2 4258 px) → raw biased low
    (11.16) and a large +0.31 fill-in. The bulge+disk was never physical for this elliptical; it was
    a kludge to tame the residual growth.
  - **(B) single-Sérsic with a broken fit (amplitude→0) + color gate** — DISCARDED artifact: the
    near-zero model made the residual term meaningless and collapsed the fill-in to ~0.01 mag
    (apparent 11.34 ± 0.013). Caught via the black model panel in the diagnostic figures.
  - **(C) single-Sérsic (fixed init) + color-gated IR growth + WCS registration** — FINAL: tighter,
    arc-focused, correctly-registered masks (F150W2 7698 px, F322W2 1997 px; companion now masked),
    raw central 11.36 **consistent with the validated expert hand-mask 11.33**, fill real but smaller
    (−0.19 to −0.28), masking systematic ±0.086 (vs ±0.16 for 2-comp). Physically motivated (single
    Sérsic) + physically-gated (source color) + astrometrically registered (§2.1.1b).

- **Files** (in sibling dir `../velocity_dispersion_from_IFU/`):
  `AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits` (500×500, 25″×25″);
  the mask was revised 2026-02-28 (`_mask_20260228.fits`); original kept as `_mask_OLD.fits`.
- **JWST bands:** same arc geometry; the HST F200LP mask is re-projected from the HST WCS
  onto each JWST frame. [METHODS §II.4]
- **For the IFU/σ_e application** (context): F200LP mask reprojected to the 0.30″ IFU grid via
  `scipy.ndimage.map_coordinates(order=0)` nearest-neighbor; **35/154 = 23% of spaxels (28% by
  I-weight)** flagged inside R_e=2.097″ (a geometric LOWER bound — the 1.27″ seeing spreads the arc
  much wider; see §2.4.1 step 3 + `scripts/psf_mask_fraction.py`). [feedback_masking_strategy.md]

### §2.1.2 Estimating R_e

- **Headline R_e = 2.097″ = 14.76 kpc** at z = 0.67564 (best-mask F140W+F200LP CoG). [CLAUDE.md; METHODS §II.3/A3]
- **"Best mask throughout" (2026-06-11):** the headline R_e is now the **best-mask** (single-Sérsic + color/morph gate + WCS registration) curve-of-growth = **2.097″**, superseding the expert-mask 2.305″. The best masks genuinely remove outer diffuse-arc + interloper light that slightly inflated the expert-mask R_e; the model-based Sérsic r_eff (1.897″) and sky-subtracted CoG (1.922″) both confirm the inward shift. The HST 2-R_e companion mask is empty (companions live only in deep JWST), so it does not enter the HST R_e — the global color/morph mask alone gives 2.097″ (self-consistent fixed point). `scripts/re_mask_sensitivity.py`, `scripts/Re_bestmask_reconciliation.py`.
- **photutils validation (`scripts/validate_Re_photutils.py`):** our azimuthally-averaged masked CoG is independently reproduced by photutils `RadialProfile` (integrated) to **±0.002″**, and unmasked matches photutils `CurveOfGrowth` (direct 2D aperture sum) to ±0.004″. A naive masked aperture-*sum* (photutils `CurveOfGrowth` dropping masked pixels) would bias R_e **+0.25–0.41″ high** → our azimuthal-fill is the correct masked-CoG treatment. `centroid_2dg` is robust; `centroid_quadratic` fails on extended light (expected — a point-source tool).
- **Method systematic (best mask) = ±0.100″** (raw CoG 2.097 / sky-sub CoG 1.922 / Sérsic r_eff 1.897, peak-to-peak/2; `results/Re_cog_reconciliation_bestmask.npz`).
- **Method:** mean of the **F140W and F200LP best-mask curves-of-growth**
  (`scripts/final_sigma_e.py:curve_of_growth`, the PRODUCTION script):
  - best-mask F140W CoG = **1.912″**; F200LP CoG = **2.281″**; mean = **2.097″** (expert-mask values were 2.168/2.441/2.305).
  - Masked CoG zeros mask-flagged HST pixels before the radial flux integral, so the lensed
    arc / companions do not bias the half-light point.
  - Centers: HST-mean `photutils.centroid_2dg` centers (same as the σ_e pipeline), both bands.
- **Do NOT confuse with `scripts/measure_Re.py`** → `results/Re_measurements.npz`. That script
  uses the image-geometric center and produces **10 variants** (3 sources × 4 masking strategies:
  `unmasked`, `zeroed`, `proper`, `PSF-convolved`); its "proper-mask" mean (2.633″) is **not the
  headline** — early-reference only. (`measure_Re.hst_Re` previously hardcoded pixscale=0.08 for
  both bands; **fixed 2026-06-08 (A3c/nb14) to read the WCS pixscale + `r_max_arcsec=6.0`, now
  reconciling with the production CoG to <0.04″** — but `final_sigma_e` remains the headline path.)
- **Superseded prior values:** IFU white-light R_e = 2.61″ (nb05/06 era); expert-mask HST CoG = 2.305″ (pre-2026-06-11). A Ca H+K+G depth-map gives 2.90″ (kinematic-tracer estimator, not the photometric R_e).
- **R_e as a σ_e systematic (test D7 → best-mask grid):** σ_e rises monotonically with R_e. Re-measured at the wide window over a 7-point R_e grid bracketing the 2.097″ headline (`scripts/run_sigma_e_Re_grid.py`); the **adopted** R_e-source systematic = **±5.13 km/s** (best-mask CoG family {1.912/2.097/2.281″}), folded into the budget (§3.1, §2.4.3). The full grid incl. CaHK+G 2.90″ (±9.98) and CaHK+G's +12.4 km/s deviation are cross-checks, NOT folded. [CLAUDE.md]

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
  **NOTE (updated 2026-06-10):** this analytic-total refit is superseded by the per-aperture
  **single-Sérsic fill-in** under the principled color-gated IR-extended masks (§2.1.1b / §3.2): the
  headline is the empirical raw value **11.36 (10%)** with a one-sided **+0.14** dex fill-in
  systematic reaching **11.47**. The 11.40 analytic-total sits just inside that bracket and the
  central is consistent with the validated expert hand-mask 11.33. The enhanced comparison figure
  `results/figures/nb08_sersic_vs_aperture.png` now overlays, for both the empirical-aperture and
  Sérsic-total choices, the best-fit Bagpipes model SED, the measured photometry, and the
  filter-convolved model photometry.
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

- **Method — single/multi-line Gaussian fits** to the integrated-spectrum features on **NIST air
  rest wavelengths**, per-line z inverse-variance combined (`scripts/redshift_verify.py`, notebook 04).

**Deflector absorption features (systemic z = 0.67564 ± 0.00033; 6 fits):**

| Feature | Ion/species | rest λ_air (Å) | obs λ_air (Å) | z_line | σ (Å) |
|---|---|---|---|---|---|
| Ca K | Ca II | 3933.66 | 6591.79 | 0.67574 ± 0.00096 | 11.5 |
| Ca H | Ca II | 3968.47 | 6649.61 | 0.67561 ± 0.00097 | 16.6 |
| Hδ | H I | 4101.74 | 6871.80 | 0.67534 ± 0.00058 | 3.3 |
| G-band | CH (band) | 4304.40 | 7209.31 | 0.67487 ± 0.00160 | 11.3 |
| Hβ | H I | 4861.33 | 8149.30 | 0.67635 ± 0.00082 | 2.8 |
| Ca H+K | Ca II doublet | 3933.66 / 3968.47 | — | 0.67567 ± 0.00069 | 11.2 |
| **Weighted** | | | | **0.67564 ± 0.00033** | |

**Source emission features (lensed arc, z = 1.30263 ± 0.00003; [O II] doublet, red cube):**

| Feature | rest λ_air (Å) | obs λ_air (Å) | z_line |
|---|---|---|---|
| [O II] 3726 | 3726.03 | 8579.2 | 1.30357 ± 0.00005 |
| [O II] 3729 | 3728.82 | 8585.6 | 1.30185 ± 0.00005 |
| **[O II] doublet** | | | **1.30251 ± 0.00004** |

  - **Deflector systemic z = 0.67564** (supersedes AGEL DR2 z = 0.67511; ~95 km/s offset, used only
    in old notebooks). [METHODS §I.2; reference_redshift_verification.md]
  - Other absorption features in `ABSORPTION_LINES` (Hγ 4340.47, Fe4383 4383.55, Mg I/Mg b 5175.27,
    Fe I 5270.40, Na D 5893.00) are NIST-air but were NOT used in the systemic fit (out of the red-cube
    range at z=0.676, or blended); listed for completeness in `scripts/redshift_verify.py`.

> **Air↔vacuum audit (2026-06-11 — PASS).** (1) Every rest wavelength is the **NIST air** value,
> verified to **<0.01 Å** (Ca K 3933.66 vs NIST 3933.663, Hβ 4861.33 vs 4861.333, [O II] 3726.03 vs
> 3726.032, etc.). **G-band 4304.40 is the CH molecular bandhead** (a blend centroid, not a single
> atomic line) → it is the most discrepant/uncertain line (z=0.67487, ±0.0016) and is down-weighted
> by its error. (2) The KCWI cube is **vacuum** (FITS `CTYPE3='WAVE'`; `AWAV` would denote air) and is
> converted to **air** before fitting via the **Ciddor (1996)/Morton (2000) IAU-standard** index of
> refraction n = 1 + 8.34254×10⁻⁵ + 2.406147×10⁻²/(130−s²) + 1.5998×10⁻⁴/(38.9−s²), s = 10⁴/λ_vac —
> **the same formula the σ_e ppxf pipeline uses** (`bootstrap_ppxf.py`; Ciddor 1996). (3) Thus
> z = λ_obs,air/λ_rest,air − 1 combines **both wavelengths in the air frame** → internally consistent,
> z **unbiased**. Had we left the cube in vacuum and used air rest lines, z would be biased high by
> ≈ (1+z)·(n−1) ≈ +5×10⁻⁴ (≈ +90 km/s); the conversion removes this.
- **Double verification with ppxf:** the Gaussian line-fit z is independent of, and consistent with,
  the ppxf V_sys (ppxf returns z via z = (1+z₀)·exp(V/c) − 1; Cappellari 2023, the ppxf
  relativistic velocity convention — *eq. number not yet verified, do not assert "eq. 5c"*).
- **Source z = 1.30263 ± 0.00003** — measured from the **[O II] λλ3726,3729 doublet** in the red cube
  (notebook 04; obs 8579.2/8585.6 Å), consistent with AGEL DR2 z = 1.302. (The ±0.00003 is the
  formal doublet-fit error; the per-line spread 1.30185–1.30357 is the realistic uncertainty.)
- **Source rest-UV resonance lines (arc ISM; `scripts/source_uv_redshift.py`, 2026-06-11).** The
  bright star-forming arc shows the classic Mg II / Fe II resonance absorption:
  - **Mg II λλ2796.35,2803.53 (red arm):** strongly detected (**16.7σ / 18.4σ**), z(absorption) =
    **1.30113 ± 0.00003**, i.e. a **−195 km/s** centroid offset vs the systemic [O II] z=1.30263 (both
    in the *same* red cube, so it is a real *relative* offset, not cross-arm calibration). **We do NOT
    interpret this as an outflow** without further work — for a resonance line on a *lensed* arc the
    offset is degenerate among: (i) Mg II **emission infill** filling the red absorption wing (P-Cygni;
    pulls the centroid blueward with no bulk motion), (ii) **differential lensing / aperture** — the
    Mg II absorption (continuum-weighted) and [O II] emission (H II-region-weighted) sample different
    source regions that the lens positions/magnifies differently, and (iii) the blended-doublet
    single-Gaussian centroid on an asymmetric profile. Report it as an **observed Mg II–[O II] velocity
    offset**, with outflow as one (unconfirmed) possibility.
  - **Fe II λλ2344.21,2382.77 (blue arm):** detected (7.9σ / 3.0σ; Fe II 2374 only 1.2σ), but at z ≈
    1.306, ~+470 km/s from the red systemic. Because Fe II sits in the **separate blue cube** (an older
    `AWAV`/air reduction, ≠ the red `WAVE`/vacuum mtwdo cube), this offset is dominated by a **blue↔red
    wavelength-calibration difference**, so the Fe II *velocity* is unreliable cross-arm — the detection
    is real, the velocity is a `[TODO: cross-calibrate the blue arm]`. (Rest λ are Morton 2003 vacuum;
    the blue cube is converted air→vac for the fit. NOT a clean independent z — quote Mg II / [O II].)
- **Companion ellipticals — a galaxy GROUP at the deflector redshift** (notebook 04a,
  `scripts/redshift_verify.bootstrap_redshift`, N=200 wild bootstrap on the red cube). Two early-type
  companions flank the primary deflector, **both at z ≈ 0.676 = the deflector z** (physically associated,
  not interlopers):

| Companion | position | sep. from deflector | z (absorption) | z (emission) |
|---|---|---|---|---|
| **NE** ("2 o'clock", larger) | IFU (x,y)=(61,64) | ≈ 4.5″ ≈ 31 kpc NE | 0.67580 −0.00067/+0.00062 | 0.67565 −0.00045/+0.00060 |
| **SW** ("7 o'clock", compact) | IFU (x,y)=(38,45) | ≈ 4.6″ ≈ 32 kpc SW | 0.67588 −0.00034/+0.00047 | 0.67558 −0.00048/+0.00071 |

  Both companions sit within \|Δcz\| ≲ 60 km/s of the deflector (0.67564) → AGEL J020613−011417 is the
  central galaxy of a **compact group/triple** at z ≈ 0.676. (These are the field ellipticals masked in
  the aperture photometry; their group membership confirms they are not background interlopers.)
- **Peculiar velocities:** we **do not** correct for peculiar velocities — we use the line-fit
  systemic z directly (see gap G3; this sentence needs to be added to the draft).
- **Wavelength convention:** AIR rest wavelengths throughout (NIST-verified to <0.01 Å); KCWI native
  vacuum (`CTYPE3='WAVE'`) converted to air via Ciddor 1996 — see the air↔vacuum audit box above. Do
  not mix SDSS vacuum values (~1 Å ≈ 80 km/s systematic at 4000 Å). [feedback_wavelength_convention.md]

### §2.2.3 Stellar Kinematics & Velocity Dispersion

**Spatial regime / masking / R_e**

- Headline σ_e is measured on a **single I(r)-weighted 1-D aperture spectrum at R < R_e** — the
  co-add-then-fit method of **Cappellari et al. 2006 (SAURON IV, §2.3)**, the literature-standard σ_e
  (KH13 / Greene+20 / SAURON / ATLAS3D / MaNGA convention). **R_e = 2.097″ = 14.76 kpc** (best mask; §2.1.2).
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
- **Per-SPS V_sys subtracted before pooling** (σ is V-invariant; a single global V_sys would inflate
  σ_e by ~10–15 km/s in the discrete cross-check). **All three libraries now converge to the same
  redshift** — at the headline R_e=2.097″ aperture the per-SPS V_sys are FSPS −15.7 / EMILES −19.1 /
  XSL −14.9 km/s → z = 0.67555 / 0.67553 / 0.67556 (a ~4 km/s, Δz≈3×10⁻⁵ spread; mean **z_ppxf =
  0.67555 ± 0.00001**). This convergence is *the* signature that the air↔vacuum frame fix worked: the
  legacy "convert all three to air" left FSPS (the vacuum-native library) ~90–110 km/s discrepant from
  EMILES/XSL; per-SPS native framing collapses it to the ~4 km/s seen here (residual folded into the
  ±5 km/s frame systematic).
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

- **σ_e definition — the luminosity-weighted second moment of the LOSVD within R_e.** (Defined as the
  radial integral in **Gültekin et al. 2009, eq. (1)** [verified 2026-06-12]; the relation we place it on
  is **Kormendy & Ho 2013 eq. 3**.) Two mathematically-equivalent representations:
  - **Analytic / luminosity-weighted integral — Gültekin et al. 2009 eq. (1)** (their 1D form;
    generalised to the 2D axisymmetric area element dA = 2πr dr):

    $$\sigma_e^2 = \frac{\int_0^{R_e} (V^2 + \sigma^2)\, I(r)\, 2\pi r\, dr}{\int_0^{R_e} I(r)\, 2\pi r\, dr}.$$
  - **Discrete spaxel/bin sum** — the same moment evaluated over spaxels (or Voronoi bins) n, with
    $F_n$ = the I-band flux (luminosity weight) of spaxel n; the 2πr factor is captured implicitly
    because the number of spaxels per annulus scales with area:

    $$\sigma_e^2 = \frac{\sum_{n:\, R_n \le R_e} F_n\,[(V_n - V_\mathrm{sys})^2 + \sigma_n^2]}{\sum_{n:\, R_n \le R_e} F_n}.$$

    This discrete spaxel/bin sum **is Cappellari et al. 2013 (ATLAS3D XV, MNRAS 432, 1709) eq. 29**
    — $\langle v_{\rm rms}^2\rangle_e \equiv \langle V^2+\sigma^2\rangle_e \equiv \sum_n F_n(V_n^2+\sigma_n^2)/\sum_n F_n$
    — equivalently the discretised form of the Gültekin 2009 eq. 1 integral. (Their eq. 29 also carries
    an **inclination/projection element** in the dynamical-model framework; our §7 uses the
    *observed-plane* projected $V_n,\sigma_n$ directly and omits that term.) **It is NOT "Cappellari
    2006 eq. 1"** — a prior draft mis-cited it there; Cappellari 2006 eq. 1 is the aperture correction
    $(\sigma_R/\sigma_e)=(R/R_e)^{-0.066}$, and its σ_e (SAURON IV §2.3) is the co-add-then-fit method
    below. [eq. numbers verified 2026-06-12: Cappellari 2006 eq. 1, Cappellari 2013 eq. 29, Gültekin 2009 eq. 1]
- **Both forms are implemented and cross-checked here** (do not re-derive). The **headline** (§6cum/nb09)
  is the **co-add-then-fit method of Cappellari et al. 2006 (SAURON IV, §2.3)** — co-add the I-weighted
  R<R_e spectra into one effective spectrum and fit a single LOSVD with pPXF; the width of that co-added
  spectrum *equals* the luminosity-weighted second moment above (this is what KH13/Greene+20/ATLAS$^{\rm 3D}$
  report). The **discrete annular sum (Gültekin 2009 eq. 1)** is the §7 / nb07c cross-check
  (`run_gultekin_mc`, F_j = Σ I over each unmasked annulus). The two agree at **<1σ**: co-add-fit
  (§6cum, narrow) 267.32 ± 24 vs discrete annular (§7, arc-filtered to R<R_safe=3R_e/4) 256.17 −13.0/+12.7.
  Full implementation + subtleties in `reference_gultekin_implementation.md`.
- **★ What the headline code ACTUALLY computes** (`final_sigma_e.extract_aperture_spectrum` → pPXF;
  this is the number we report, 267.31). **We do NOT evaluate either σ_e² sum or integral above.**
  We form the **I(r)-weighted co-added aperture spectrum**

  $$S_e(\lambda) \;=\; \sum_{n:\, r_n < R_e} w_n\, S_n(\lambda), \qquad
    w_n = \frac{I_n}{\sum_{m:\, r_m < R_e} I_m}, \qquad w_n \equiv 0 \ \text{for arc-masked spaxels (mask\_weight=0)},$$

  where $S_n(\lambda)$ is the spaxel spectrum and $I_n$ its 6500–7500 Å observed-band IFU flux (the
  luminosity weight). **The choice of weight map $I_n$ is itself a budgeted systematic — the "I-shape"
  term (±2.27 km/s; §2.4.3, `run_isource_shape_sweep.py`).** Holding the spaxel selection fixed, we
  re-measure σ_e with **10 alternative weight maps** — (1) the headline 6500–7500 Å IFU band, (2) full
  IFU white-light, (3–4) HST F140W / F200LP reprojected (raw), (5–6) the same arc-masked, (7–8) their
  1-D curve-of-growth annular means, (9–10) their 2-D Sérsic-model fits — and take peak-to-peak/2 of the
  spread (266.83–271.37 km/s on the NEW cube + M10). So $I_n$ = IFU white-light is the *headline* choice,
  not an assumption: σ_e moves ≤±2.27 km/s across this whole family of luminosity weightings.
  pPXF then fits a **single Gaussian LOSVD** (`moments=2`) to $S_e(\lambda)$,

  $$S_e(\lambda) \;\approx\; \Big[\textstyle\sum_k a_k\, T_k(\lambda)\Big] \otimes \mathcal{G}(v;\,V,\sigma)
    \;+\; \text{(additive Legendre poly, deg 15–29)},$$

  and **$\sigma_e \equiv \sigma$**, the fitted LOSVD dispersion. Because superposing spaxel spectra at
  different line-of-sight velocities $V_n$ broadens the co-add, the width $\sigma$ of that single fit
  *equals* the luminosity-weighted second moment $\sqrt{\langle (V-V_{\rm sys})^2 + \sigma_n^2\rangle}$
  — so it **measures** the Cappellari 2013 eq. 29 / Gültekin 2009 eq. 1 quantity **without ever forming
  the explicit sum** (and with $V_{\rm sys}$ absorbed automatically into the fitted $V$; §2.4.1). Code:
  `extract_aperture_spectrum` (the weighted sum) + `bootstrap_ppxf.setup_ppxf_inputs_from_spectrum` +
  `ppxf(..., moments=2)`. The explicit discrete sum $\sum_j F_j(V_j^2+\sigma_j^2)/\sum_j F_j$ is computed
  **only** in the §7 annular cross-check (`run_gultekin_mc`), never in the headline.
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
   through the KCWI WCS to a sub-pixel IFU centre. R_e = mean of the F140W (1.912″) + F200LP
   (2.281″) **best-mask** curves-of-growth = **2.097″** (`final_sigma_e.curve_of_growth`, r_max=6″;
   photutils-validated ±0.002″).
3. **Spatial arc mask → IFU grid.** The F200LP `_mask.fits` (2512 HST px, arc-only) is reprojected
   onto the 0.30″ IFU grid with `scipy.ndimage.map_coordinates(order=0)` (nearest-neighbour) →
   `arc_spax_mask`. **Mask coverage within R_e (PSF caveat; `scripts/psf_mask_fraction.py`):** this
   *geometric* mask flags **35/154 = 23% of spaxels (28% by I-weight)** inside R_e=2.097″. But the arc
   light is spread by the **1.27″ KCWI seeing**, so the geometric fraction is a **lower bound**: growing
   the mask by the PSF (binary dilation) raises the masked fraction to **71% / 81%** (by FWHM/2 = 0.64″)
   or **96% / 98%** (by the full FWHM = 1.27″) — i.e. a hard seeing-grown mask removes nearly the entire
   aperture and is **infeasible**. (Convolving the binary mask with the PSF and cutting at >25% / >10%
   arc-light contamination gives 42% / 75% by count; the >50% cut, 1%, is a thin-arc dilution artifact.)
   This is precisely why we instead (i) **I-weight** the aperture — the bright deflector core dominates,
   the arc-contaminated outskirts are down-weighted — and (ii) carry the **mask-weight systematic
   (±6.65 km/s, w=0→1; §2.4.3)** to bracket the residual arc dilution rather than attempt a hard
   PSF-grown mask.
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
   the §2.4.3 systematic components in quadrature → **σ_e = 267.31 −12.98/+12.61 (sym ±12.79) km/s**
   (at the best-mask R_e=2.097″ aperture; cache `results/run_wide_sigma_e/resys_best_mean/`).

### §2.4.2 Photometry / M⋆ pipeline — step by step

Drivers: `scripts/photometry_systematics.py` (masking + photometry), `bagpipes_sersic_refit.py`
(SED), orchestrated/displayed in notebook 12. Recipe = **F200LP locates, IR extends, fill-in recovers**.

1. **Locate the arc on F200LP** (best source contrast): a 2D Sérsic deflector model is subtracted
   and a pixel is flagged as arc where its positive residual exceeds **k·σ_sky** — i.e. **k is the
   detection threshold in units of the sky noise σ_sky** (so k=3 means residual > 3σ_sky), with a
   parallel S/N floor (the pixel's own S/N > SNR_thresh). We adopt **k = 3**
   (`arc_mask_verification.py`): this **objective Sérsic-residual mask reproduces the expert hand mask
   to 0.016 mag**, R_e to 3%. **k = 3 is the clean saturation point** of a 2-D sweep over the two
   thresholds (SNR_thresh ∈ {2..20} × k ∈ {2..8}) — below k≈3 it grabs noise, above it the masked-flux
   fraction stops changing.
2. **Reproject the F200LP footprint to every band** (WCS + `map_coordinates`) — reproduces the
   expert JWST aperture mags to 0.01–0.02 mag.
3. **IR-extend with a COLOR GATE.** In the deeper/sharper JWST bands, region-grow the arc seed into
   pixels that are **(a)** significant positive residual above the single-Sérsic deflector model
   (data−model > 3σ_sky), **(b)** bluer than the deflector core by ≥0.5 mag in F200LP−F140W, and
   **(c)** contiguous with the seed within R<5″. The **color gate (b) is essential**: a single Sérsic
   slightly under-fits the bright IR body, so residual-only growth would sweep the *red* galaxy body
   and companions into the mask; the *blue* z=1.302 source vs *red* passive deflector contrast keeps
   only genuine source pixels. See §2.1.1b for the detailed rationale + caveats.
4. **Deflector model = single-component Sérsic** — physically motivated for the elliptical (per user
   2026-06-10); not a bulge+disk decomposition. Fit gotcha: `Sersic2D.amplitude` is I(r_eff) not the
   peak — init ~peak·exp(−b_n), bound to (peak·1e-3, peak·1.5), else LevMar collapses to ≈0.
   PSF-convolved fill quantified at ≤0.004 mag (`psf_fill_model.py`) — negligible.
5. **Two photometry estimators per band/mask:** **raw** `photutils` aperture (masked → discards
   under-arc deflector light → biased low, mask-size-dependent) vs **single-Sérsic fill-in** (replace
   masked pixels with the model → recovers it). Fill-in correction **−0.19 to −0.27 mag** under the
   color-gated masks (smaller than the old 2-component +0.18–0.96 mag — cleaner masks remove less).
6. **AB magnitudes** from per-image header zeropoints (never hardcoded): HST `PHOTFLAM/PHOTPLAM`,
   JWST `PIXAR_SR`. Four bands: F200LP, F140W (HST WFC3), F150W2, F322W2 (JWST NIRCam).
7. **Bagpipes SED** (`02_streamlined_Bagpipes_SED`): exponential-τ SFH, Calzetti dust (A_V 0–2),
   free Z (0–2.5 Z☉), z prior (0.674,0.676); **Nautilus** sampler n_live=400, 500-sample posterior.
   A **10% fractional flux floor** per band covers photon noise + the drizzle pixel-correlation
   under-estimate (empirical noise ≈1.3–2× CCD-equation); **20%** is run as a conservative variant.
8. **M⋆ budget** = 8 fits {per-band, global} × {raw, filled} × {10%, 20%} (script-reproducible via
   `scripts/mstar_masking_budget.py` → `Mstar_masking_systematic.npz`, `photometry_systematics_Mstar.npz`)
   → headline (empirical, raw-central with the best per-band color-gated mask, one-sided +sys, user
   choice): **log M⋆ = 11.36 ± 0.08 (stat) +0.14 (sys) at 10%** [11.26 ± 0.14 +0.18 at 20%], fill-in
   reach 11.47; explicit masking-approach systematic **±0.086 dex** (was ±0.16 under the 2-component
   masks). Central consistent with the validated expert hand-mask 11.33.

### §2.4.3 Systematic-audit procedures (how each budget term was measured)

Each σ_e systematic is a **dedicated sweep** re-using the §2.4.1 machinery, varying ONE axis and
quoting **peak-to-peak/2** (or half-Δ for two-point axes), all on the headline cube + masks:

- **I-shape (±2.27):** 10 alternative I-weight maps (`run_isource_shape_sweep.py`) × 3 SPS × 15 deg
  × N=250; peak-to-peak/2 of the 10 central σ_e (266.83–271.37).
- **F200 mask-weight (±6.65):** down-weight the arc spaxels at w∈{0,0.5,1} × 3 SPS × N=500;
  (w00 269.69 − w100 256.39)/2. *(This is the mask-**strength** axis. The mask-**definition** axis —
  expert/sersic/perband/global reprojected to the IFU grid, `run_sigma_e_mask_systematic.py` — gives
  ±4.58, which overlaps this; larger-of-two is kept, not added.)*
- **Fit-window (±3.82):** 3 windows (wR3800_5400 / wR4000_5400 / w6500_7500) × 3 SPS × N=500;
  peak-to-peak/2 (269.66 / 268.58 / 276.22).
- **Frame (±5.0):** structural — vac/air per-SPS native-frame choice.
- **Centering (±4.0):** 5 perturbed HST-mean centres (±0.4″ sweep; `NOTES_centering_investigation`).
- **Reduction-pass (±3.45):** half-Δ between the NEW and OLD cleaned cubes (269.62 vs 262.72); 2
  reductions only — refine if a 3rd lands.
- **R_e-source (±5.13):** σ_e re-measured over a 7-point R_e grid bracketing the 2.097″ best-mask
  headline (`run_sigma_e_Re_grid.py`). **Adopted** = best-mask CoG family {1.912/2.097/2.281″}
  peak-to-peak/2 = (269.99−259.73)/2. The full grid incl. CaHK+G 2.90″ (±9.98) and CaHK+G's +12.4
  deviation are cross-checks, NOT folded (CaHK+G is a different I-map definition, not an error on
  the photometric R_e). Was ±6.13 at the old expert-mask headline.

Quote the budget table (§3.1) as the result; the iterative history that *led* to this masking +
component set is the **internal audit trail (M8→M12)** at the end of this document — not for the
manuscript Methods section.

### §2.4.4 Reproducibility / infrastructure

- **Architecturally-independent estimators** (must not share intermediate state): **§6cum** cumulative
  I-weighted single-ppxf aperture (the Cappellari 2006 co-add-then-fit method = headline) · **§7** discrete
  **Gültekin (2009) eq. 1** annular spaxel-sum (`run_gultekin_mc`, arc-filtered to R<R_safe) · **nb07e**
  arc-spectrum subtraction. The three agree at <1σ (267.32 co-add-fit vs 256.17 discrete, narrow).
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

- **Final redshifts (to quote in Results).** The deflector systemic redshift from **Gaussian line
  fits** (6 absorption features, §2.2.2) is **z = 0.67564 ± 0.00033**, in agreement with the
  independent **ppxf stellar-kinematic** redshift **z = 0.67555 ± 0.00002** (V_sys = −16 ± 4 km/s
  relative to the line-fit value; headline R_e=2.097″ aperture). **All three SPS libraries converge:**
  FSPS 0.67555 / EMILES 0.67553 / XSL 0.67556 (~4 km/s spread) — the signature of the air↔vacuum
  per-SPS frame fix (§2.4.1), which collapsed the legacy ~110 km/s FSPS-vs-air offset. The 16 km/s
  line-fit↔ppxf offset is well within the ±59 km/s line-fit uncertainty (and the ppxf
  formal ±4 km/s carries an additional ~±5 km/s per-SPS frame systematic). **We adopt z = 0.67564**
  (the line-fit value supersedes AGEL DR2 z = 0.67511). The lensed-**source** redshift is
  **z = 1.30263 ± 0.00003** (red-cube [O II] λλ3726,3729 doublet). [§2.2.2; `scripts/redshift_verify.py`,
  `run_sigma_e_Re_grid.py`]
- **Headline: σ_e(<R_e) = 267.31 ± 12.79 km/s** (asym −12.98/+12.61), the luminosity-weighted σ_e
  within R_e (Gültekin et al. 2009 eq. 1 definition) measured by the **Cappellari et al. 2006 (SAURON IV
  §2.3) co-add-then-fit method** — a single pPXF fit to the I(r)-weighted R<R_e aperture spectrum
  (R_e = 2.097″ best-mask CoG) — plus the wild-bootstrap error pool. **Supersedes 269.62 ± 13.27 at R_e=2.305″** — the "best mask throughout"
  cascade (2026-06-11) adopts the best-mask R_e=2.097″, lowering σ_e by 2.3 km/s along the
  rising σ(R) profile (well within errors).
- **Reproduce:** σ_e at the headline aperture is `results/run_wide_sigma_e/resys_best_mean/`
  (R_e=2.097″), produced by `scripts/run_sigma_e_Re_grid.py --n_bootstrap 500`; the R_e-source
  systematic (±5.13, best-mask CoG family) is from the same grid (`results/sigma_e_Re_grid_N500.npz`).

**R_e measurement and checks** (for the §3.1 paragraph): R_e = 2.097″ = 14.76 kpc, the F140W (1.912″)
+ F200LP (2.281″) **best-mask**-CoG mean (photutils-validated ±0.002″); supersedes expert-mask 2.305″
and IFU-only 2.61″; masked CoG removes arc/companion bias.

<!-- PV:auto:budget -->
**σ_e systematic budget** — *generated by `scripts/paper_values.py --render`; do not hand-edit inside the markers.*

| Component | ± km/s | Note |
|---|---|---|
| stat (N=500) | 4.51 | asym −5.04/+3.98; pooled 3 SPS × 15 deg (marginalizes SPS + degree) |
| I-shape | 2.27 | 10-shape peak-to-peak/2 on NEW cube + M10 |
| F200 mask | 6.65 | (w00-w100)/2; larger-of-two vs mask-approach 5.85 (no double-count) |
| frame (vac/air) | 5.00 | carried constant |
| centering | 4.00 | carried constant |
| fit-window | 3.82 | peak-to-peak/2 across 3 fit windows |
| reduction-pass | 3.45 | half-Δ between NEW and OLD reductions |
| R_e-source (best-mask grid) | 5.13 | peak-to-peak/2 across best-mask CoG family {F140W 1.91", mean 2.10", F200LP 2.28"} at the 2.097" headline (user 2026-06-11; CaHK+G & full grid are cross-checks, not folded) |
| **TOTAL (sym)** | **12.79** | quadrature (sys ±11.96 ⊕ stat) |
| **TOTAL (asym)** | **−12.98 / +12.61** | preserves stat-side skew |

Cross-checks (not added): arc-mask-definition ±4.58 (overlaps F200 mask) · full-grid R_e-source ±9.98 (incl. CaHK+G, conservative ceiling) · CaHK+G(2.90″) deviation +12.38 km/s.
<!-- /PV:auto:budget -->

*(Per-component sweep detail — the per-window / per-estimator values behind each ± — is in §2.4.3.)*

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

**Compare with Greene (2016/2020) and KH13 choices:** σ_e(<R_e) here uses the luminosity-weighted
single-aperture σ_e (Cappellari et al. 2006 SAURON IV §2.3 co-add-then-fit; Gültekin 2009 eq. 1 definition),
exactly as adopted by Kormendy & Ho 2013
(M•–σ eq. 3) and Greene et al. 2020 ARA&A (fig. 5). We dropped R<R_e/8 (JFK95 σ_c) because seeing
under-resolves it. For a pure elliptical, galaxy R_e ≈ bulge R_e, so KH13's "bulge R_e" equals the
measured total R_e — no bulge/disk decomposition needed.

**σ_e physical-meaning discussion:** the deflector is a dispersion-supported **elliptical bulge**;
σ(r) should *decrease* gently with radius — any rising σ(r) in raw output is a continuum-dilution
artifact from the z=1.302 arc, not a physical gradient. Resolved kinematics (nb11): 17 central
spaxels at S/N≥5, σ ∈ [144, 251] km/s (median 201), declining with R, no coherent rotation above noise.

## §3.2 Stellar Mass Estimate

- **HEADLINE: log(M⋆/M☉) = 11.47 +0.09/−0.15** (aperture-corrected total, **matched 2 R_e aperture**,
  10% flux floor; `scripts/aperture_matched_photometry.py`, `results/aperture_matched_photometry.npz`).
  Total-light mass: empirical flux in a **matched 2 R_e elliptical aperture** (same physical region all
  4 bands; single-Sérsic deflector model; color/morph-gated WCS-registered arc mask + companion mask),
  corrected to total by **adding the Sérsic model's beyond-aperture light** (the GAMA/Taylor+2011
  `fluxscale` approach; additive variant). Methodologically: **Taylor et al. 2011** (GAMA aperture→total
  for M⋆), **Sonnenfeld et al. 2013** (SL2S lens-deflector Sérsic photometry, mask arc + interlopers),
  **Graham & Driver 2005** (Sérsic total formalism).
- **Report FIVE estimators (empirical → model), M⋆ for each, at 2 R_e:**
  | estimator | log M⋆ | what it adds |
  |-----------|--------|--------------|
  | raw (empirical aperture, masked) | **11.22** | nothing — lower bound |
  | raw + apcorr (most-empirical total) | **11.35** | model **wings** only (masked interior NOT filled) |
  | filled (within aperture) | **11.36** | model **fill** of masked pixels |
  | **total (aperture-corrected) — headline** | **11.47** | fill + wings |
  | Sérsic-total (pure model) | **11.41** | everything modeled |
  The **total converges across 1/2/2.5 R_e** (11.45–11.49) → the aperture correction is internally
  consistent. raw converges 2↔2.5 R_e (11.22).
- **Sérsic-only (full-model) systematic budget = ±0.19 dex** (`scripts/sersic_total_systematic.py`;
  elliptical multi-start fit, 2026-06-11): stat ±0.12 ⊕ sys ±0.15, with **arc-mask choice ±0.125
  dominant**, model-form ±0.057, Sérsic-fit-n ±0.027, flux-floor ±0.023, apcorr↔pure-model ±0.010.
  (Was ±0.17 with the circular fit; the mask term grew because the elliptical model widens the
  expert↔global Sérsic-total gap.)
- **Apcorr chain validated vs established codes (2026-06-11, `scripts/validate_apcorr_established.py`;
  petrofit/statmorph absent → photutils + scipy incomplete-gamma references):** (1) b_n Ciotti&Bertin99
  vs exact `gammaincinv(2n,½)` ≤0.05%; (2) the Sérsic total-flux formula (Graham&Driver05 eq.4) vs a
  numeric render-to-20 R_e ≤0.03%; (3) the aperture correction — rendered-model enclosed light vs the
  **analytic Sérsic incomplete-gamma law** `γ(2n, b_n·2^{1/n})/Γ(2n)` (Graham&Driver05) — Δ≤0.0007;
  (4) the empirical hard-edge aperture sum (`in_aperture`) vs photutils `EllipticalAperture(method='exact')`
  ≤0.09%; (5) the finite-cutout truncation of `sum(model_full)` (the apcorr denominator) vs the
  to-∞ total ≤0.19% (≤0.002 mag). **No production bug; M⋆ unchanged.**
- **Single-Sérsic ellipticity fix (2026-06-11).** The deflector light model is a single Sérsic per band
  fit on the best mask. The fitter originally had a single low-ellipticity start that railed into a
  spurious **circular** χ² minimum (b/a=1 for all bands), under-reporting the shape and biasing the
  `(1−ellip)` total flux. Fixed with a **multi-start** (flux-weighted-moment seed à la SExtractor/GALFIT
  + a PA grid + circular fallback, lowest-residual). It recovers the true shape; **M⋆ headline (total)
  is unchanged at 11.47** (empirical-aperture-dominated). Validated by injection-recovery against
  astropy's peer-reviewed `Sersic2D` (`scripts/validate_sersic_fitter_synthetic.py`): clean recovery for
  the deflector regime (n≈1.2–1.6, b/a≤0.85), with a realistic **b/a uncertainty ~±0.06** (the formal
  bootstrap is too tight). GALFIT/imfit/petrofit/statmorph are not installed (conda-env policy).
- **The M⋆ posterior has a longer LOW tail** (e.g. total 11.47 +0.10/−0.13). This is **NOT** the
  total-Sérsic estimator — it is present in all five estimators including the purely empirical `raw`.
  It is the **age–dust–M/L "outshining" degeneracy** in the 4-band Bagpipes fit: the low-M⋆ tail
  correlates with young age (r=+0.90), high dust (r=−0.58), and high sSFR (r=−0.83) — young+dusty
  low-M/L solutions reproduce the same broadband points. We keep the exponential-SFH prior and report
  the asymmetric posterior as-is (the tail is real model uncertainty, not an estimator artifact).
- **Per-band single-Sérsic parameters (appendix table; `scripts/sersic_parameter_table.py`):**

| Band | $\lambda_{\rm piv}$ (Å) | $r_e$ (″) | $r_e$ (kpc) | $n$ | $b/a$ | PA (°) | $\mu_e$ (mag/″²) | $m_{\rm AB}^{\rm tot}$ |
|---|---|---|---|---|---|---|---|---|
| F200LP | 4972 | 2.002±0.122 | 14.09±0.86 | 1.24±0.15 | ≳0.95† | —† | 25.23 | 20.92±0.07 |
| F140W | 13923 | 1.765±0.015 | 12.42±0.11 | 1.44±0.15 | 0.86±0.06 | 4±20 | 22.36 | 18.42±0.01 |
| F150W2 | 16720 | 1.643±0.005 | 11.57±0.03 | 1.22±0.15 | 0.80±0.06 | 11±20 | 21.91 | 18.29±0.00 |
| F322W2 | 32470 | 1.757±0.010 | 12.37±0.07 | 1.59±0.15 | 0.85±0.06 | 6±20 | 21.79 | 17.82±0.01 |

  † F200LP: the deflector is faint in this band, so the single-Sérsic ellipticity/PA are unconstrained
  (consistent with circular); the IR bands set the shape. Errors are bootstrap, with $b/a$/PA/$n$ floored
  by the injection-recovery scatter; the mask-choice systematic (§2.1.1b) dominates and is separate.
  The deflector is **mildly elliptical (b/a≈0.80–0.86, PA≈4–11° E of N)** with low Sérsic index
  (n≈1.2–1.6) — consistent with the b/a=0.75 adopted for the aperture from the independent light fit.
- **Aperture provenance / why this supersedes 11.36:** the old per-band apertures were hand-drawn in the
  GUI and **inconsistent** — F200LP a=1.19″ (~0.5 R_e) captured only **17%** of the light vs IR ~46% —
  so F200LP was under-measured → SED artificially red → M⋆ biased **HIGH** (11.36). The matched aperture
  measures the same region in every band → bluer SED → lower M/L → the empirical raw drops to 11.22.
- **Flux uncertainties:** the *propagated statistical* errors are **0.03–1.6%** (F322W2 0.03%, F150W2
  0.05%, F140W 0.4%, F200LP 1.6%) — 6–362× below the adopted **10% floor**. The 10% is a deliberate
  **systematic floor** (zeropoint, aperture, drizzle correlated noise ~1.3–2×, SED-model mismatch), not
  the measurement error — report the true errors and state the floor (~5% is better-motivated; precedent
  in GAMA/Taylor+2011, whose error budgets are calibration-dominated, not photon-noise). **Floor sensitivity (`results/aperture_floor5_check.npz`):** 5%↔10% shifts every estimator's central by ≤+0.04 dex (headline total 11.47→11.50) with tighter stat at 5% — well within the ±0.12 systematic, so the headline is robust to the floor; 10% kept for the headline.
- **Superseded values (cross-checks):** mismatched-aperture 11.36 (F200LP-biased high); 2-component-mask
  11.16 +0.31 (over-masked → biased low); expert-aperture 11.33.
- **"Typical elliptical at z~0.7":** M⋆ ≈ 2.3 × 10¹¹ M☉ (raw headline) to ≈ 3.0 × 10¹¹ (fill-in upper) at z = 0.6756, passive/quiescent SFH
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
- **σ_e definition:** **Gültekin et al. (2009)** eq. 1 (the luminosity-weighted second-moment integral —
  verified 2026-06-12); co-add-then-fit method **Cappellari et al. (2006, MNRAS 366, 1126, SAURON IV §2.3)**
  [NB: that paper's eq. 1 is the *aperture correction*, NOT σ_e — do not cite "Cappellari 2006 eq. 1" for
  σ_e]; the discrete spaxel/bin sum (§7) is **Cappellari et al. (2013, ATLAS3D XV, MNRAS 432, 1709) eq. 29**
  (= the discretised Gültekin 2009 eq. 1; eq. 29 also has an inclination/projection term that our
  observed-plane sum omits) [verified 2026-06-12]; **Kormendy & Ho (2013, ARA&A 51, 511)** eq. 3 (scatter 0.29 dex);
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

1. **R_e method family (was "non-converged"; largely ADDRESSED 2026-06-11).** Headline R_e = 2.097″
   (best-mask raw CoG, r_max=6″) is now **photutils-validated to ±0.002″** (`validate_Re_photutils.py`:
   independent `RadialProfile` integration), and the raw CoG sits at the **top** of a tight best-mask
   method family — sky-sub CoG 1.922″, Sérsic r_eff 1.897″ → method systematic **±0.100″** carried
   explicitly. The raw-CoG-is-top-of-family caveat remains (no sky pedestal subtraction), but it is now
   bounded and documented, not a flagged unknown. Knock-on: R_e sets the σ_e aperture and (weakly) M★;
   both effects are inside their systematic budgets.
2. **R_e-source budget term (RESOLVED 2026-06-11).** Now uses the **best-mask CoG light family only**
   {1.912/2.097/2.281″} → ±5.13 — the referee-preferred light-family treatment. CaHK+G 2.90″ (a
   different I-map definition) and the full grid (±9.98) are reported as cross-checks, NOT folded, so
   the budget no longer conflates aperture-uncertainty with σ(R)-physics.
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

> **Superseded by "best mask throughout" (2026-06-11).** The numbers in THIS appendix (σ_e=269.62,
> total ±13.27, R_e-source ±6.13, R_e=2.305″) are the **M12 (2026-06-08) state**. They were the
> headline until the best-mask cascade adopted R_e=2.097″ (best-mask CoG, photutils-validated),
> which moved σ_e → **267.31 ± 12.79** and re-derived R_e-source → **±5.13** (best-mask CoG family;
> `run_sigma_e_Re_grid.py`). The current state is §2.1.2 / §2.4.3 / §3.1 above. This trail is kept
> verbatim as M12 provenance; do not quote its numbers as current.

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
  cross-check (±4.58; overlaps the F200-mask term, not double-counted) and resolved the Hδ open item
  (**keep unmasked** — well-fit, masking it discards LOSVD information). Central σ_e 269.62 unchanged
  across M10→M12; only the error grew.

**Net effect of the trail:** central σ_e settled at 269.62 km/s by M10 and has not moved since; the
error bar evolved ±17.78 (pre-M11) → ±11.77 (M11) → ±13.27 (M12). The mask set (He I + 35-entry
sky/CR + no-Balmer) and the 7-component budget that the final §2.4 procedure uses are the M12 result.
