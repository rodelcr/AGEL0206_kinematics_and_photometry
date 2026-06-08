# Methods and Systematics — AGEL0206 σ_e ApJL paper

**Last updated:** 2026-05-18 (narrative); **headline refreshed 2026-06-08.**

> **⚠ CURRENT HEADLINE (2026-06-08, post-M12) — authoritative source `results/PAPER_VALUES.json`:**
> **σ_e(<R_e) = 269.62 ± 13.27 km/s** (asym −13.45/+13.10) · **log M★/M☉ = 11.16 +0.32/−0.08** ·
> **R_e = 2.305″ = 16.23 kpc** · z_defl = 0.67564 · z_src = 1.302.
> The narrative BODY below predates the M8–M12 mask/systematic audits and prints earlier values
> (254.85 / 268.98, total ±18) in prose — those are HISTORICAL. For current numbers read this banner,
> `CLAUDE.md`, or `DRAFTING_FACTS_paper_2026-05-29.md`; all are emitted by `scripts/paper_values.py`.

This document is a single-pass narrative of the spectroscopic
(σ_e) and photometric (R_e, M★) pipelines used for the AGEL J020613−011417
deflector, with the systematic-error budget worked into the text rather
than relegated to a separate appendix. It is intended to:

- serve as a methods-section draft for the ApJL paper;
- give a new collaborator (or a future Claude session) a complete picture
  of how every headline number was produced and what was checked;
- complement `TESTS_AND_DIAGNOSTICS.md`, which is the long-form catalog
  of every individual test (A1–L3) with status and cache paths.

Where a test code (e.g. **D3**, **J5**) appears in brackets, it points
to the matching row of `TESTS_AND_DIAGNOSTICS.md` §0 for the full
provenance.

---

## 0. Headline numbers

| Quantity | Value | Source / notebook |
|---|---|---|
| **σ_e(<R_e)** — **paper headline**, wide arc-masked window, NEW `_mtwdo_` reduction, bad-pixel mask + Balmer kept + He I 3819 mask + M10 full sky audit + **M11 rigorous re-derivation of carried systematics** | **269.62 ± 11.77 km/s** (symmetric; asymmetric form **−11.98 / +11.57** preserves the stat-side bootstrap skew). M11 (2026-05-27) replaced the carried I-shape (±1.5), F200-mask (±3.8), and fit-window (±15) components with cube+mask-matched values (±2.27 / ±6.65 / ±3.82). The fit-window drop is the dominant change (the 3 windows now agree to ±4 km/s on the cleaned NEW cube). Reduction-pass component ±3.45 (M10). | `scripts/run_wide_sigma_e.py --cube new_clean_hei`, caches at `results/run_wide_sigma_e/new_clean_hei/`. M11 sweep caches at `results/{ishape_sweep,maskweight_sweep}_wR3800_5400_arcmask_M10/` + `results/nb09a_wavelength_sweep_M10/` |
| σ_e(<R_e) — old reduction cross-check, cleaned + He I + M10 sky masks | 262.72 ± 17.99 km/s (asymmetric −18.10 / +17.88) | `scripts/run_wide_sigma_e.py --cube headline_clean_hei`, caches at `results/run_wide_sigma_e/headline_clean_hei/`. The 6.90 km/s Δ to the cleaned + He-I + M10-sky-masked new cube is the source of the refined ±3.45 km/s reduction-pass systematic |
| σ_e(<R_e) — pre-M10 cross-check (cleaned + He I, NEW cube) | 271.87 −17.99/+17.74 km/s | Pre-M10 reference. +2.25 km/s shift to current headline = bias from un-masked OH/sky residuals folded into the M10 systematic |
| σ_e(<R_e) — pre-He-I cleaned cross-check (NEW cube) | 268.98 −18.19/+17.98 km/s | `--cube new_clean`. +2.89 km/s offset (un-masked He I 3819 contamination biased σ_e *down*) |
| σ_e(<R_e) — pre-He-I cleaned cross-check (OLD cube) | 260.44 −18.36/+17.97 km/s | `--cube headline_clean` |
| σ_e(<R_e) — legacy un-cleaned wide-window cross-check (NEW cube) | 265.76 ± 17.87 km/s | Pre-clean reference; `scripts/run_wide_sigma_e.py --cube new`, caches at `results/run_wide_sigma_e/new/` |
| σ_e(<R_e) — legacy un-cleaned wide-window cross-check (OLD cube) | 254.85 ± 17.87 km/s | Pre-clean reference; `scripts/run_wide_sigma_e.py --cube headline`, caches at `results/run_wide_sigma_e/headline/` |
| σ_e(<R_e) — narrow Ca H+K+G cross-check window (legacy old cube, no clean) | 267.95 ± 30.10 km/s | nb09d (`results/sigma_e_final_systematics_nb09d.npz`) |
| σ_e(<R_e/2) — gradient diagnostic | ~225 km/s (NARROW) | nb07c §6cum |
| **R_e** (paper headline) | **2.305″ = 16.23 kpc** | mean F140W + F200LP masked CoG; `scripts/final_sigma_e.py:curve_of_growth` |
| **log(M★/M☉)** (paper headline, aperture SED) | **11.33 +0.07 / −0.09** | Bagpipes, nb02 |
| log(M★/M☉) Sersic-total cross-check | 11.40 +0.11 / −0.15 | Bagpipes refit, nb08 |
| z_deflector | 0.67564 | nb04 multi-line Gaussian fit (NIST air rest λ) |
| z_source | 1.302 | AGEL DR2, independently cross-checked via Fe II UV multiplet in KCRM blue arm (I.2.2) |

The headline σ_e and σ_e cross-check are consistent at 0.4σ
(13.1 km/s difference, both budgets cover it).

---

# Part I — Spectroscopic analysis (σ_e)

## I.1 KCWI integral-field spectroscopy

**Cube actually used for σ_e (KCRM red arm):**
`Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`
(265 MB, shape 3317 × 100 × 100, λ = 5625–8941 Å, Δλ = 1.0 Å/pix,
0.300″ spaxels, ~30″×30″ FoV). Symlinked at the repo root from
`../velocity_dispersion_from_IFU/`.

**Complementary blue-arm reduction:**
`raw_KCWI/blue/DESJ0206_medium_combined.fits` (299 MB, shape
2595 × 120 × 120, λ = 3278–5872 Å, Δλ = 1.0 Å/pix, 0.300″ spaxels,
~36″×36″ FoV). Same nights as the red cube, KCRM Blue-Low grating.
**Not used for σ_e** (the deflector absorption features that anchor
σ_e are in the red), but used for the source-redshift cross-check
in I.2.1 below — at z_source = 1.302 the blue cube covers source-rest
~1424–2551 Å, which contains the Fe II UV resonance lines.

**Provenance correction — verified against the Keck observing logs (2026-05-26).**
The cube's FITS header (`DATE-OBS = 2025-11-17`, `PROGNAME = U002`, `PROGPI =
Jones`, `DATE-BEG = 2025-11-17T09:02:15`) describes only the FIRST input frame
of the most recent contributing night (`kr251117_00129`, the first DESJ0206
RED frame on Nov 17), not the full input set. The cube comment
*"Reduced and stacked by Kaustubh Rajesh Gupta at UTC 2026-01-07 00:39:38"*
records the actual stacking date and reducer. The science data is a
multi-night combine of **two confirmed AGEL collaboration nights** plus
a possible third:

| Night (UT) | Program | PI | DESJ0206 frames | Source |
|---|---|---|---|---|
| **2025 Aug 30** | K409 | TBD | 12 RED `kr250830_00090–00101` (60 min) + 4 BLUE `kb250830_00052–00055` (66 min) | local `raw_KCWI/provenance/K409_2025-08-30_UTC.html` (Drive id `1yoH_elZqEsFHPRc6f-XCB98l-7dGJewN`) |
| **2025 Nov 17** | U002 | Jones | 20 RED `kr251117_00129–00148` (90 min) + 5 BLUE `kb251117_00087–00091` (110 min) | Drive `NightLog_KCWI_2025-11-17.pdf` (id `1HVNoH2CSAc214b_PieBQ_UkeEasVl6q4`) |
| 2024 Dec 29 | TBD | TBD | 4 RED `kr241229_00092–00095` per local `dec29_2024_rl8000.list` — **NOT independently confirmed as DESJ0206 pointings**; no machine-readable observing log exists in Drive (only image-only PDF `16zWO8yaCIvRtoCGenNTndvB7dYdrCkMH`) | (raw frames on Yuguang Chen's machine; alignment QA PNGs in Drive id `1ELcTQiXV9m04Y8sA9st3OlpsNk7LE0D7`) |

Total confirmed on-source integration: **≈ 170 min RED + 176 min BLUE** over
Aug 30 + Nov 17. RL grating central wavelengths differ between nights
(Aug 30: 7400 Å; Nov 17: 7150 Å) — both rebinned onto a common
1.0 Å/pix grid spanning 5625–8941 Å by the kcwiRedux pipeline.

The user (Cordova) is on the observer list for **both** Aug 30 (K409 crew:
Alcorn, Cordova, Glazebrook, Gonzalez-Lopezlira, Jones, Kacprzak, Tran,
Vasan G C) and Nov 17 (U002 crew: Alcorn, Barone, Chen, Cordova,
Glazebrook, Gupta=Kaustubh, Jones, Kacprzak, Rhoades, Tran, Vasan G C).
A previously-recorded "Sept 29 K409 contribution" is **incorrect** — the
local Sept 29 log shows zero DESJ0206 entries. **For the paper, cite the
two confirmed nights** (Aug 30 K409 + Nov 17 U002) and flag the Dec 29
frames as a possible third-night contribution pending header verification
on the raw FITS files. See `reference_kcwi_data_properties.md` and **A1**.

### Instrumental line-spread function (**B1**)

`DISPSCAL = 0.2940` (instrumental σ in pixels) × `CD3_3 = 1.000 Å/pix`
→ **FWHM_inst = 0.692 Å** (constant in observed Å). Across the fit
window the instrumental dispersion is σ_v,inst ≈ 12.6 km/s observed
(R ≈ 10000), or 7.5 km/s rest at z_def = 0.67564. The deflector velocity
dispersion (~270 km/s) is therefore resolved by a factor of ~20.

The instrumental LSF is passed to ppxf as a per-pixel rest-frame FWHM
dict (`{"lam": lam_gal_rest, "fwhm": fwhm_gal_rest}` at
`scripts/bootstrap_ppxf.py:_prep_spectrum_for_ppxf`), the canonical
Cappellari (2023) interface. `sps_lib` interpolates this onto each
template's native wavelength grid and pre-broadens via
`util.varsmooth`. We verified (**B2**) that scaling DISPSCAL across the
range {×0.5, ×0.75, ×1.0, ×1.25, ×1.5, ×2.0} shifts σ by at most 0.83
km/s — most templates are intrinsically broader than the KCWI LSF, so
ppxf's `fwhm_diff² = (FWHM_gal² − FWHM_tem²).clip(0)` clamp at
`sps_util.py:169` means no convolution is applied for FSPS or EMILES, and
only sub-km/s deconvolution for XSL. The LSF subtraction is rock-solid
relative to the headline precision.

### Seeing

Guide-camera FWHM = 1.27″ (`GUIDFWHM` header from the Nov 17 stacking
pass). Used as a conservative seeing estimate; the K409 Aug 30 log gives
airmass 1.07–1.08 but the per-night DIMM seeing has not yet been
extracted from the observing-log header (open item flagged in
`reference_kcwi_data_properties.md`). The R<R_e/8 aperture (0.288″) was
dropped because it sits inside the seeing FWHM/2 = 0.64″ (**D5**).

## I.2 Redshifts

### I.2.1 Deflector systemic redshift (**A2**)

The deflector redshift was redetermined by independent Gaussian fitting
of seven absorption / emission lines in the integrated spectrum,
anchored on NIST **air** rest wavelengths (Ca K 3933.66, Ca H 3968.47,
G-band 4304.40, Hδ, Hγ, Hβ, [OII] doublet with shared parameters). The
result, z_systemic = **0.67564**, supersedes the AGEL DR2 catalog value
0.67511 (which is offset by ~95 km/s — used only in the older notebooks
03/03b). Notebook 04; module `scripts/redshift_verify.py` provides the
`ABSORPTION_LINES`, `EMISSION_LINES`, and fitting routines.

All wavelength bookkeeping downstream uses NIST air rest wavelengths;
see I.5 for how the air↔vacuum convention is handled between the galaxy
spectrum and the SPS template libraries.

### I.2.2 Source redshift cross-check via Fe II UV absorption (blue arm)

The AGEL DR2 source redshift z_source = 1.302 was independently
cross-checked using the **KCRM blue-arm cube**
(`raw_KCWI/blue/DESJ0206_medium_combined.fits`, λ_obs = 3278–5872 Å, see
I.1). At z = 1.302 this observed-frame range covers source-rest
~1424–2551 Å, containing the canonical Fe II UV resonance multiplets
that are the standard low-ionization absorption-line diagnostic in
intermediate-z galaxy ISM spectra. The Fe II UV3 / UV2 lines fall
cleanly inside the blue cube; the strong UV1 doublet (Fe II λ2587, λ2600)
lies marginally past the blue-arm red edge at 5872 Å (would appear at
obs ~5957 / 5988 Å):

| Source-rest line | Source-rest λ (Å, air) | Observed λ at z=1.302 (Å) | In blue cube? |
|---|---|---|---|
| Fe II UV3 | 2344.21 | 5397.4 | ✓ |
| Fe II UV2 | 2374.46 | 5466.9 | ✓ |
| Fe II UV2 | 2382.76 | 5486.0 | ✓ |
| Fe II UV1 | 2586.65 | 5955.5 | (just past red edge) |
| Fe II UV1 | 2600.17 | 5986.6 | (just past red edge) |

These features, fit jointly on the arc spectrum extracted from the blue
cube, recover a source redshift consistent with the AGEL DR2 value
z_source = 1.302. **[TBD: insert measured z_source ± uncertainty,
notebook/script reference, and result-cache path once the analysis is
committed — the work was done on the blue arm but the artifact is not
yet in the repo as of 2026-05-19.]**

The Fe II UV multiplet identification is the cleanest available
verification for z_source in this dataset: the source is a star-forming
spiral (so the resonance Fe II absorption is expected to be present);
all five lines have published air-rest wavelengths to ≪ 0.01 Å (NIST);
and joint fitting with shared (z, σ) collapses the per-line statistical
uncertainty by √N. The Mg II λ2796, λ2803 doublet from the same source
(red cube, def-rest 3835–3855 Å) is independently the strongest
detection used for the source-emission spectral mask (I.4), giving a
second handle on z_source via the same physical absorber.

## I.3 Aperture and I-weight construction

The headline σ_e is measured on a single I-weighted, 1-D aperture
spectrum extracted at R<R_e from the IFU cube. This matches the
literature standard (Cappellari et al. 2006, eq. 1, as adopted by
Kormendy & Ho 2013, Greene et al. 2020 ARA&A, SAURON/ATLAS3D, MaNGA).

### Center and aperture radius

- **Center.** F140W and F200LP cutouts are centroided with
  `photutils.centroid_2dg` on the deflector core (each with its own
  `_mask.fits`), averaged in world coordinates, then propagated to IFU
  sub-pixel via the KCWI WCS (`scripts/final_sigma_e.py:find_center`).
  Adopted: **RA = 31.55613°, Dec = −1.23819°** (ICRS). The F140W vs
  F200LP centroid offset is 0.36″ — driven by the F200LP arc dragging
  that filter's centroid; F140W is the cleaner bulge centroid, and the
  HST-mean is the robust compromise. The 5-center sweep (**F1**) gives
  σ_e spread = 3.7 km/s, carried as the ±4 km/s "centering" component
  of the systematic budget.

- **Aperture radius.** R_e = 2.305″ = 16.23 kpc, the mean of the F140W
  and F200LP masked curves-of-growth (see II.3). Three apertures were
  considered: R<R_e/8, R<R_e/2, R<R_e (**D5**). R<R_e/8 was dropped
  because at 0.288″ it sits inside seeing FWHM/2; the headline is at
  R<R_e (KH13 / Greene+20 σ_e definition).

### Spatial arc mask (**E1**)

The lensed source at z=1.302 contributes flux at several IFU spaxels
along the arc that pass into the R<R_e aperture as low-velocity
continuum. We mask them using the **F200LP segmentation mask**
(`AGEL020613-011417A_F200LP_WFC3_cutout_L3_mask.fits`, 500 × 500 at
HST 0.05″/pix). This mask is arc-only by construction; the F140W
`_mask.fits` is broader and also covers the deflector core, so it is
*not* an arc mask and never used as one (see `feedback_masking_strategy.md`,
**A4**/I2 in TESTS_AND_DIAGNOSTICS).

The HST mask is reprojected to the IFU 0.300″ grid via
`scipy.ndimage.map_coordinates(order=0)` (nearest-neighbor). Inside
R<R_e: 38 of 184 spaxels are flagged (~21% by count, ~27% by I-weight).

### I-weight construction

For each spaxel inside R<R_e (and unflagged by the spatial mask), the
flux weight is the per-spaxel collapse of the IFU 6500–7500 Å white-
light band (the same window used in the narrow-window ppxf fit, so the
weighting and the fit have matched footprints). The weighted 1-D
spectrum is summed (not averaged) so the wild bootstrap can preserve
the per-spaxel noise contribution.

The choice of I-weight map is one of the systematic axes (**D3**, **J6**).
We bracket it with the 10-shape sweep
(`scripts/run_isource_shape_sweep.py`): IFU_band (headline), IFU
white-light, F140W and F200LP × {raw, arcmasked, 1D-CoG smooth,
2D-Sersic smooth}. After the 2026-05-11 Sersic2D bound-fix (J5; see
I.5 below) the I-shape spread is **±1.5 km/s** at the wide arc-masked
window (was ±13 km/s at narrow before the bound-fix removed the
F200LP n=0.30 pathology).

## I.4 Wavelength fit window

This is the methodological pivot of the paper: there are two windows
in play, and the headline moved to the wider one during May 2026.

### Narrow Ca H+K + G window (cross-check)

`w6500_7500` = obs 6500–7500 Å = **rest 3879–4476 Å** at z = 0.67564.
Anchored on Ca H+K and the G-band. This is the window used by the older
nb07c §6cum pipeline (now retained as the **method cross-check**).
927 good pixels. σ_e(<R_e) = 267.95 ± 30.10 km/s. <!-- pv-skip: narrow-window method cross-check, not the headline -->

### Wide arc-masked window (current paper headline)

`wR3800_5400_arcmask` = **rest 3800–5336 Å**, with explicit masks at four
deflector-rest bands that are contaminated by emission from the z=1.302
source (factor (1 + z_s)/(1 + z_l) = 1.374 maps each source-rest line into
the deflector rest frame). 2161 good pixels — a 2.6× increase in spectral
information vs the narrow window. σ_e(<R_e) = 254.85 ± 17.87 km/s.

**Spectral-mask catalog (J3, codified at
`scripts/run_window_sweep.py:ARC_MASKS_REST`): three source-emission
masks plus one telluric.**

| Identification | Source-rest λ (Å) | Deflector-rest λ (Å) | Obs λ (Å) | Masked band, def-rest (Å) |
|---|---|---|---|---|
| Mg II λ2796, λ2803 (z=1.302) | 2796, 2803 | 3842, 3852 | 6438, 6455 | 3835–3855 |
| **O₂ A-band leading-edge telluric** | — | — | **7593–7626** | 4525–4545 |
| [O II] λ3727, λ3729 (z=1.302) | 3727, 3729 | 5121, 5124 | 8581, 8585 | 5115–5135 |
| [Ne III] λ3869 + Mg b cluster (z=1.302) | 3869 | 5314 | 8903 | 5260–5340 |

Total: 212 of 2374 pixels (8.9% of the fit window) excluded. The catalog
was identified empirically (>4σ excesses above local continuum in the
wild-bootstrap posterior, **J2**). Mg II, [O II], and [Ne III] are
physically identified at z = 1.302. The fourth band, originally tagged
"source rest ~3300 Å unidentified," was confirmed on 2026-05-18 to be a
**telluric residual at the leading edge of the atmospheric O₂ A-band**
(obs 7593–7626 Å): the spike sits at the same observed wavelength in
the deflector aperture, the CLAUDE.md sky box, and a 4–8″ off-deflector
annulus, with deflector/off-source amplitude ratio 1.11× despite a
17.5× continuum-brightness ratio — the unambiguous signature of an
additive wavelength-locked systematic, not a localized source line. See
`NOTES_4534A_spike_investigation_2026-05-18.md` and
`results/figures/spike_4534A_diagnostic.png`. The mask was placed for
the right reason (kills the spike) but for the wrong attribution; the
mask itself is unchanged.

This catalog is the first-of-its-kind methodology for AGEL deflector
kinematics. It is reusable for any AGEL target where z_source is known.

**Why the wide window is the headline (not just incrementally better).**
Three independent improvements stack at the wide arc-masked window
(documented in nb09d and figured in
`results/figures/nb09d_per_sps_both_windows.png`):

1. **4× tighter statistical error.** 2161 good pixels vs 927; the
   wild-bootstrap stat 1σ shrinks from ±23.9 to ±6.1 km/s (**J6**).
2. **SPS-template spread collapses 26 → 4 km/s.** At the narrow window
   the three SPS libraries disagree by 26 km/s with strong ordering
   (FSPS=253.7 < EMILES=267.9 < XSL=279.6; see also
   `reference_sps_systematic.md`). At the wide arc-masked window they
   agree to 4.2 km/s (FSPS=253.9, EMILES=253.3, XSL=257.5) (**J4**).
   Mg b + Fe5270 + Fe5335 + the broader feature set anchor the templates
   well enough that the FSPS/EMILES/XSL ordering vanishes.
3. **F200 spatial-mask sensitivity drops 4×.** From ±7.5 km/s (narrow)
   to ±3.8 km/s (wide) (**J7**) — the spectral arc-emission masks
   absorb most of what the F200 spatial mask used to absorb; the
   spatial mask becomes a sub-budget perturbation.

**Provenance of the wavelength range.** Not from a specific literature
prescription. The 3800 Å blue edge is set by the XSL template floor at
def-rest ~3675 Å plus a 5% velocity padding; the 5336 Å red edge is set
by KCWI red-channel coverage at z = 0.67564. The closest literature
precedent for going as blue as 3800 Å is LEGA-C (Bezanson et al. 2018;
van der Wel et al. 2021) at z ~ 0.7 using ~3500–5500 rest.
SAURON/ATLAS3D's canonical kinematic window is narrower (4800–5380; Mg b
+ Fe5270 + Fe5335 only). The Lick-tradition "rest 4000–5400 (no Ca H+K)"
range is captured by the `wR4000_5400_arcmask` sensitivity check (see
I.6 below), not the headline. **A drafted statement claiming "Cappellari
2017 recommends extending blueward" was a hallucination — retracted.**

## I.5 ppxf setup

We use `ppxf` v9.x (Cappellari & Emsellem 2004; Cappellari 2017; 2023)
with three SPS template libraries pooled at the posterior level.

### SPS libraries and the vacuum/air fix (**B3, C1–C4**)

All three libraries shipped with ppxf are used and pooled:

| Library | File | Native wavelength frame | V_sys diagnostic |
|---|---|---|---|
| FSPS | `spectra_fsps_9.0.npz` | **vacuum** | Ca K minimum at 3934.86 Å (+1.20 Å vs air rest) |
| EMILES | `spectra_emiles_9.0.npz` | **air** | Ca K minimum at 3933.82 Å (+0.16 Å) |
| XSL | `spectra_xsl_9.0.npz` | **air** | end-to-end V_sys closure at +9 km/s with air galaxy |

The original pipeline (nb07c era) converted the galaxy to air for *all
three* libraries. Air-galaxy ↔ vacuum-FSPS produced a spurious −90 km/s
FSPS V_sys offset (matches the air↔vacuum differential at 6500–7500 Å,
−83 km/s per Ciddor 1996). The fix
(`scripts/bootstrap_ppxf.py:SPS_NATIVE_FRAME` + `frame_galaxy='auto'`):
keep the galaxy in vacuum for FSPS, apply the scalar-median
`util.vac_to_air` ratio for EMILES and XSL. After the fix the V_sys
split-track collapses from ~110 km/s to ~15 km/s and the σ shift on the
3-SPS pool is ≤ 2 km/s (**C2**, **C3**). The full derivation is in
`NOTES_nb09_frame_fix_and_final_sigma_e_2026-04-28.md`, including the
canonical Cappellari pattern citation (`ppxf_example_kinematics_sdss.py`
:127-134) and Ciddor (1996) for the conversion accuracy.

We follow Cappellari's pattern of applying `vac_to_air` at observed
wavelengths (rather than rest); the alternative would shift V_sys by
1.83 km/s at our z (**B5**), sub-budget.

### Per-SPS V_sys subtraction before pooling (**C4**)

The σ posteriors from each SPS are concatenated raw (σ is invariant to
V offsets), but V_sys medians are subtracted per-SPS before pooling so
that residual template-systematic offsets (now FSPS −19, EMILES −4,
XSL −1 km/s post-frame-fix) don't inflate the V posterior. The pooled σ
posterior is the 16/50/84 percentile of the concatenated all-3-SPS
sample (typically ≥10⁴ samples after 30 polynomial degrees × N=500
bootstrap × 3 SPS).

### Polynomial degree (**G3**)

Multiplicative polynomial degrees 15–29 (15 values) are stacked into
the per-SPS posterior. The σ–vs–degree diagnostic
(`results/figures/nb09_sigma_vs_degree.png`) is flat within the
bootstrap envelope across this range — the polynomial is saturated and
no monotonic trend (which would indicate the polynomial absorbing LOSVD
signal) is present. The lower bound (deg=15) is set by the requirement
to flatten the global continuum across the wider arc-masked window;
deg<15 leaves systematic residuals.

### Wild bootstrap (**G1, G2**)

For each (SPS × degree) we run a hybrid Rademacher × local-residual wild
bootstrap (`scripts/bootstrap_ppxf.py:run_bootstrap_single_degree`). Per
iteration: residual = (galaxy − bestfit) × Rademacher sign, locally
scaled by σ from a 75-pixel rolling window, re-fit, record V and σ.
Production statistics is N=500 per (SPS × degree); N=50 smoke runs are
used for development. N=50 vs N=500 agree on σ within 1 km/s (**G4**).
Parallelization uses joblib with BLAS pinned to 1 thread per worker
(`scripts/bootstrap_ppxf_parallel.py`) — without the pin, BLAS thread
oversubscription on Mac multicore caused a >4 h pathology on one fit.

## I.6 σ_e measurement

The pipeline applied per window:

```
1. Extract the I-weighted aperture spectrum at R<R_e (I.3).
2. For each SPS:
     a. Set the galaxy wavelength frame to match the SPS native frame (I.5).
     b. For each polynomial degree d ∈ {15, …, 29}:
          Run wild bootstrap × N=500.
3. Subtract per-SPS V_sys median; concatenate σ samples across SPS × degrees.
4. Headline σ_e = median of the pooled posterior; stat 1σ = pooled (84−16)/2.
```

The headline σ_e values are reported in §0. Full per-SPS,
per-aperture, per-polynomial-degree posteriors are cached in
`results/nb09a_wavelength_sweep/{window}_{sps}_T*_N500.npz` for the
wide window, and `results/final_sigma_e_paper/` for the narrow window.

### Three-window cross-check (**J8**)

A third window beyond the two headlines, `wR4000_5400_arcmask` (rest
4000–5400, no Ca H+K), was run at full N=500 as an orthogonal
feature-set check (Hβ + Mg b only, the Lick-tradition window). Its σ
agrees with both headline windows within the budget:

| Window | rest pixels | σ_e (3-SPS pool, N=500) | role |
|---|---|---|---|
| w6500_7500 | 3879–4476 | 269.7 km/s | NARROW cross-check |
| wR3800_5400_arcmask | 3800–5336 | **254.8 km/s** | **HEADLINE** |
| wR4000_5400_arcmask | 4000–5400 | ~265 km/s | orthogonal Hβ+Mg b check |

The 15 km/s spread among these three windows defines the **fit-window
systematic** carried in both budgets (now the dominant residual).

## I.7 Independent cross-check estimators

Three architecturally-independent σ_e estimators were built and compared
at the narrow window (none share intermediate state with the cumulative
I-weighted aperture fit). They all agree within budget; see
`reference_cumulative_vs_annular_sigma_e.md` for the design rationale,
`reference_sigma_e_estimator_independence.md` for the independence
guarantee.

| Estimator | Method | σ_e(<R_e) NARROW | vs §6cum |
|---|---|---|---|
| **§6cum** (paper cross-check) | Cumulative I-weighted aperture ppxf | 267.32 ± 24 km/s | — |
| **§7** (annular Gültekin) | Discrete annular sum of F·(V²+σ²) (**H1**) | 256.17 −13.0/+12.7 km/s | −11 (<1σ) |
| **§7b** (flat-σ extrapolation) | §7 + outer-annulus σ extrapolation (**H2**) | 274.37 −16.2/+17.4 km/s | +7 (<1σ) |
| **nb07e** (arc subtraction) | §6cum with α·arc subtracted as additive template (**E7**) | matches §6cum within 0.1 km/s | — |
| **nb07f** (arc-as-sky) | §6cum no-mask + arc as ppxf `sky=` template (**E8**) | 274.64 (recovery=145%) | confirms dilution mechanism |

§7 (discrete annular sum) is the σ_e definition often associated with
Gültekin et al. (2009); we use it for the radial σ(R) profile and the
per-annulus systematics inspection. The §6cum cumulative I-weighted
aperture is the literature-standard σ_e definition for M•–σ work
(Cappellari et al. 2006, eq. 1, also used by KH13 / Greene+20 /
SAURON / ATLAS3D / MaNGA), and is the cross-check the paper headline is
calibrated against. The arc-as-sky test (E8) is the rigorous test of
the no-mask Δ mechanism: with the arc spectrum as a free-amplitude
additive template, σ recovers to 274.6 km/s, +6.8 km/s above the masked
headline — the **145% recovery fraction** confirms that the no-mask
bias is entirely continuum dilution and there is no kinematic-blend
mechanism contributing in the opposite direction. The 7 km/s overshoot
is consistent with a small amount of real deflector light leaking into
the outer-arc-mask spaxels that get oversubtracted.

## I.8 Resolved kinematics (supplementary)

The wide arc-masked window also enables per-spaxel and PowerBin
kinematic maps that failed at the narrow window (Test 19 in nb05x).
Documented in `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` and
`project_nb11_resolved_kinematics.md`.

- **Per-spaxel at S/N ≥ 5 (K1):** 17 central spaxels inside R<1.5 R_e
  with σ ∈ [144, 251] km/s, zero outliers above 400 km/s, median
  σ = 201 km/s. σ(R) drops with radius — consistent with the
  elliptical-bulge expectation (see `feedback_deflector_morphology.md`).
- **PowerBin at target S/N=15 (K2):** 7 bins inside R<1.5 R_e, median
  σ = 294 km/s. 5 of 7 inner bins are clean; the 2 outer bins still hit
  σ > 400 km/s (KCWI per-spaxel S/N limit at the aperture edge, NOT a
  binning failure).
- **K3:** at the narrow window neither approach worked — outer bins
  σ > 800, no per-spaxel S/N above floor. The wide arc-masked window
  delivered a 4× increase in per-pixel S/N that unlocked the maps.

These are supplementary figures, not the headline.

## I.9 Spectroscopic systematic budget — both windows

The headline σ_e is the median of the pooled posterior on the NEW
`_mtwdo_` reduction; the error is the quadrature sum of **seven**
independent components (one more than the original budget — the
"reduction-pass" component added 2026-05-26, refined again 2026-05-27
when the He I 3819 source-emission mask was added):

| Component | Wide arc-masked, NEW `_mtwdo_` + clean + He I + M10 + M11 (headline) | Source |
|---|---|---|
| Statistical (pooled 1σ) | **±4.6** (asym −5.16/+4.13) | wild bootstrap N=500 × 3 SPS × 15 degrees on cleaned + He-I + M10-masked NEW cube |
| **I-shape (10-shape spread, M11)** | **±2.27** | `results/ishape_sweep_wR3800_5400_arcmask_M10/` — peak-to-peak/2 of 10 I-shape posteriors (range 266.83–271.37 km/s) |
| **F200 spatial mask (peak-to-peak/2, M11)** | **±6.65** | `results/maskweight_sweep_wR3800_5400_arcmask_M10/` — (w00=269.69, w50=261.95, w100=256.39)/2 = 6.65 |
| Frame (vac/air per SPS, structural) | ±5.0 | `audit_ppxf_methodology.py` audit 1 |
| Centering (5-center sweep ±0.4″, carried) | ±4.0 | `NOTES_centering_investigation_2026-04-27.md` |
| **Fit-window (3-window spread, M11)** | **±3.82** | `results/nb09a_wavelength_sweep_M10/` — wR3800_5400_arcmask 269.66, wR4000_5400_arcmask 268.58, w6500_7500 276.22 → peak-to-peak/2 = 3.82. **Dramatic drop from carried ±15** because the cleaned NEW cube has wide/narrow agreement (vs ±15 spread on OLD cube) |
| **Reduction-pass (M10)** | **±3.45** | half-Δ between cleaned + He-I + M10-sky-masked new and old reductions = (269.62−262.72)/2 = 3.45 km/s |
| **Quadrature total (symmetric)** | **±11.77** | sqrt(2.27² + 6.65² + 5.0² + 4.0² + 3.82² + 3.45²) |
| Asymmetric total (lower / upper) | **−11.98 / +11.57** | preserves stat-side skew |

**2026-05-27 M11 rigorous re-derivation result:** total sys dropped from ±17.78 → ±11.77 km/s (down 6.0). The fit-window component dropping from ±15 → ±3.82 is the dominant change. The carried ±15 was from pre-clean OLD-cube nb09 §9; the cleaned NEW cube + M10 masks bring the three windows into much closer agreement.

OLD cube cross-check at the same M10 mask configuration (carried components still apply for it, since it's not the headline): σ_e = 262.72 ± 17.99 km/s. Narrow Ca H+K+G cross-check on OLD cube (legacy): σ_e = 267.95 ± 30.10 km/s.

The fit-window systematic (±15 km/s) is still the **dominant residual**.
The other six components quadrature-sum to ~10.5 km/s. Tightening the
window systematic would require additional N=500 windows at different
feature subsets (Hβ-only, Mg b-only) — not cheap, and only marginally
useful given the dominant role of the 13 km/s wide-vs-narrow shift.

**On using "old" component values for the new headline**: the I-shape,
F200-mask, and fit-window systematics were measured on the OLD cube.
We assume they scale similarly on the NEW cube (component-magnitude
phenomena from the same input frames + same arc + same I-weight maps,
just different flat-fielding) and carry the OLD values unchanged. Strict
re-derivation on the new cube would cost ~50 min N=500 × 3 SPS × 10
shapes + 3 mask weights — possible but expected to give the same
~1-km/s-level numbers since the underlying frames and aperture are
unchanged. Flagged for re-derivation if needed.

**Note on independence.** The budget components are not strictly
statistically independent: I-shape is bounded by mask choice; frame is
bounded by SPS choice; centering is partially absorbed by the
multiplicative polynomial. Quadrature combination is therefore
*conservative* — the true total error could be ~10% smaller. We
report the conservative number.

### Cross-references for individual budget items

- Statistical: bootstrap convergence at N=500 (**G4**); polynomial
  saturation (**G3**); per-SPS V_sys subtraction before pooling
  (**C4**).
- I-shape: 10-shape sweep (**D3**, **D4**, **J6**); Sersic2D bound-fix
  prevents the n=0.30 pathology (**J5**, applied 2026-05-11 at
  `scripts/run_isource_shape_sweep.py:181-207`).
- F200 mask: hard/no-mask Δ (**E2**, **E3**); 5-point mask-weight
  sweep with quadratic c = +7.3 (threshold-dominated bias, **E5**); arc
  dilation rejected as inflating σ via noise (**E6**, negative finding);
  arc-subtraction sibling matches §6cum (**E7**); arc-as-sky 145%
  recovery confirms dilution mechanism (**E8**).
- Frame: per-SPS frame identification + V_sys closure (**B3, C1, C2**);
  σ shift sub-budget (**C3**).
- Centering: HST-mean centroid_2dg (**A4, F2**); 5-center sweep (**F1**).
- Fit-window: nb09 §9 three-window spread (**J8**); window sweep
  underpinning (**J1**); source-emission catalog (**J2, J3**).

### What is *not* in the budget (intentionally)

- **R_e source.** Four R_e definitions (mean, F140W-only, F200LP-only,
  Ca H+K+G depth-map) give σ_e at the narrow window with 16.9 km/s
  spread (**D7**) — at the mask-budget level, but sub-budget vs the
  total ±30. σ_e rises monotonically with R_e (more outer-bulge
  spaxels in the I-weighted aperture). This is treated as a method
  sensitivity rather than a budget addition.
- **Reduction-pass (M1)** — **PROMOTED to the formal budget (2026-05-26).**
  The new `_mtwdo_` reduction (with hybrid master-twilight + dome flats)
  is now the headline σ_e source; the OLD reduction is retained as a
  second-reduction cross-check. The half-Δ between the two reductions
  (= 5.46 km/s) is carried as the **"reduction-pass" component of the
  systematic budget** (see I.9 table above). Flag: only two reductions
  available — revisit and refine if a 3rd reduction lands. Audit:
  `NOTES_nb09e_audit_2026-05-20.md` + memory
  `project_nb09e_reduction_systematic.md`; cross-test for localized
  mechanism: rejected (the +10.9 km/s shift is not from a single
  wavelength boundary — see Addendum 2026-05-26 (b) in the project memory).
- **Local-MAD bad-pixel cleaning (M5, 2026-05-26 sensitivity, NOT in budget).**
  ppxf `clean=True` rejects 0 pixels on this dataset because the noise
  array is overestimated relative to the actual residual scatter. Using
  local-MAD-based outlier detection (rolling 75-pixel window, |residual|/
  local_MAD > 3σ) on the canonical residuals identifies **52 cosmic-ray-
  like spikes** that the noise-scaled clean=True missed. The biggest is a
  6-pixel cluster at def-rest 5232–5236 Å (obs 8768–8774 Å) with a
  26σ amplitude in local-MAD units — a clear unresolved CR residual.
  Clipping these 52 pixels and re-running at N=100 shifts σ_e UP by:
    - NEW cube: 265.76 → **267.18 km/s** (Δ = +1.42)
    - OLD cube: 254.85 → **256.01 km/s** (Δ = +1.16)

  **Both cubes shift by the same amount → the 52 bad pixels are
  intrinsic to the data (CR residuals shared by both reductions),
  NOT reduction-specific** (M6, 2026-05-26 replication). The +10.9 km/s
  reduction-pass gap is preserved (10.91 un-clipped → 11.17 clipped).
  An alternative noise-scaled approach (scale noise array by 0.3 and
  re-run ppxf clean=True) finds only 19 pixels per fit but gives a
  larger σ_e shift (+5.28 km/s to 271.04) — different pixel-selection
  algorithm; the local-MAD result is more defensible because the same
  flagged pixels replicate the cleaning shift on the OLD cube.

  We carry this as a stated sensitivity rather than folding into the
  formal budget because (a) the shift is well below stat 1σ (±6.1),
  (b) cleaning preserves the reduction-pass gap, (c) the choice of
  threshold (3σ local-MAD) is not strongly motivated vs alternatives.
  Caches: `results/run_wide_sigma_e/new_local_mad_clip_N100/`,
  `results/run_wide_sigma_e/old_local_mad_clip_N100/`, and
  `results/run_wide_sigma_e/new_noise_scaled_clean_N100/`. Scripts:
  `/tmp/test_local_mad_clip.py`, `/tmp/test_local_mad_clip_OLD_cube.py`,
  `/tmp/test_noise_scaled_clean.py`.
- **R<R_e/2 vs R<R_e aperture choice.** Aperture is a methodological
  choice, not a systematic. R<R_e is the KH13/Greene+20 convention and
  the paper headline; R<R_e/2 is reported as a gradient diagnostic.
- **σ–vs–degree trend.** Flat across deg 15–29 (**G3**); no
  contribution to budget. If a future referee asks, the trend is shown
  inline in nb09 §6.5.

---

# Part II — Photometric analysis (R_e, M★)

## II.1 Imaging data

Four bands enter the SED fit and the morphological measurements; full
provenance in `reference_hst_jwst_data_properties.md`.

| Telescope | Instrument | Filter | Program | PI | Date (UT) | Drizzled scale | Exptime | ZP_AB |
|---|---|---|---|---|---|---|---|---|
| HST | WFC3/UVIS | F200LP | 16773 | Glazebrook | 2022-07-14 | 0.05000″/pix | 600.0 s | **27.344** |
| HST | WFC3/IR | F140W | 16773 | Glazebrook | 2022-07-14 | 0.08000″/pix | 597.7 s | **26.446** |
| JWST | NIRCam SW (Mod B) | F150W2 | 05594 (Cluster SLICE) | Mahler | 2024-08-27 | 0.03075″/pix | 1836 s | 28.033 |
| JWST | NIRCam LW (Mod B) | F322W2 | 05594 | Mahler | 2024-08-27 | 0.06301″/pix | 1836 s | 26.475 |

AB zeropoints are derived per image from the FITS headers (the standard
HST formulation `ZP_AB = −2.5·log10(PHOTFLAM) − 21.10 −
5·log10(PHOTPLAM) + 18.6921`; JWST `ZP_AB = −6.10 − 2.5·log10(PIXAR_SR)`
on the MJy/sr i2d products). **Zeropoints are read from headers; never
hardcoded.**

Plate scale at z_def = 0.67564 (FlatLambdaCDM, H₀=70, Ω_m=0.3):
1″ ≈ 7.04 kpc. At the four drizzled scales above this is 352, 564, 217,
and 444 pc/pix respectively.

**Drizzle pixel correlation.** Noise in drizzled images is correlated
between adjacent pixels — empirical noise is ~1.3–2× higher than the
idealized CCD-equation prediction. The STScI ETC does **not** account
for this. Relevant for aperture photometry error bars; we propagate via
the per-pixel noise map produced by the drizzle pipeline rather than
recomputing from photon statistics.

## II.2 Centering

Identical to I.3 — `photutils.centroid_2dg` on F140W and F200LP cutouts,
averaged in world coordinates. RA = 31.55613°, Dec = −1.23819° (ICRS).
The F140W/F200LP centroid offset is 0.36″ (driven by the F200LP arc
contamination); the HST-mean is the robust adopted center.

## II.3 Effective radius (**A3**)

R_e = **2.305″ = 16.23 kpc**, the mean of the F140W and F200LP
masked curves-of-growth (`scripts/final_sigma_e.py:curve_of_growth`).
The masked CoG zeros mask-flagged HST pixels before the radial flux
integral, so the lensed arc and any companion contributions do not bias
the bright-end half-light point.

| Band | Masked CoG R_e (″) |
|---|---|
| F140W | 2.168 |
| F200LP | 2.441 |
| **Mean (headline)** | **2.305** |

**Production vs measurement script** — important not to confuse them.
The production R_e (2.305″) comes from `final_sigma_e.py:curve_of_growth`
using HST-mean `centroid_2dg` centers. A separate older script
`scripts/measure_Re.py` (`results/Re_measurements.npz`) uses the cutout
image-geometric center and produces 10 variants under different masking
strategies; its "proper mask" rows give a different mean (2.633″).
Cite the production value.

The earlier IFU white-light value (R_e = 2.61″, used in nb05/06 era) is
superseded by the HST-based R_e = 2.305″, which has the F140W + F200LP
arc masks built in.

**R_e source as a σ_e systematic (D7).** Four R_e definitions (mean,
F140W-only, F200LP-only, Ca H+K+G depth-map at R_e = 2.866″) shift the
narrow-window σ_e by up to 14 km/s. σ_e rises monotonically with R_e —
larger aperture catches more outer-bulge spaxels. This sensitivity is
documented in nb09 §7b and figured at
`results/figures/nb09_re_source_systematic.png`; it is held outside the
formal budget (see I.9).

## II.4 Aperture photometry

Performed interactively with `scripts/photometry_masking_HST.py` (HST)
and `scripts/photometry_masking_JWST.py` (JWST) — interactive
matplotlib GUI tools (click-to-mask + sliders); not headless. The HST
script uses `PHOTFLAM`/`PHOTPLAM`; the JWST script uses `PIXAR_SR`.

The deflector aperture excludes the arc using the F200LP `_mask.fits`
(reused from the σ_e pipeline); for the JWST bands the arc geometry is
the same and re-projected from the HST WCS. The four extracted AB
magnitudes are the input to Bagpipes (II.5):

| Filter | λ_pivot (Å) | AB observed |
|---|---|---|
| F200LP | 4972 | 22.6126 |
| F140W | 13 923 | 19.1335 |
| F150W2 | 16 720 | 18.9425 |
| F322W2 | 32 470 | 18.6042 |

Error bars are set to **10% fractional** on each band as a conservative
systematic budget (covers both the formal photon noise and the
drizzle-correlation underestimate noted in II.1).

## II.5 Bagpipes SED fitting (**A5**)

We fit the four-band SED with `Bagpipes` (Carnall et al. 2018, 2019).
Configuration:

```python
fit_instructions = {
    "redshift":   (0.674, 0.676),       # narrow spec prior, brackets z = 0.67564
    "exponential": {
        "age":          (0.1, 15.),     # Gyr
        "tau":          (0.3, 10.),     # e-folding time, Gyr
        "massformed":   (1., 15.),      # log10(M_formed / M☉)
        "metallicity":  (0., 2.5),      # Z / Z☉
    },
    "dust": {"type": "Calzetti", "Av": (0., 2.)},
}
```

- **SPS:** BPASS + CLOUDY (Bagpipes default), Calzetti dust law, free
  metallicity.
- **Sampler:** Nautilus (n_live = 400; pure-Python fallback to
  MultiNest if needed). 500-sample posterior.
- **Cache:** `pipes/posterior/AGEL0206/0206_real.h5` (delete to refit).

The fit was executed end-to-end in
`notebooks/02_streamlined_Bagpipes_SED.ipynb` (and the older exploratory
`02_Bagpipes_SED_fitting.ipynb`). Saved posteriors in
`results/bagpipes_sed_results.npz`.

### Posteriors

| Parameter | p16 | **p50** | p84 |
|---|---|---|---|
| **log(M★/M☉)** (headline) | **11.24** | **11.33** | **11.41** |
| mass-weighted age (Gyr) | 3.69 | 5.10 | 6.11 |
| exponential SFH age (Gyr) | 4.46 | 5.98 | 6.97 |
| exponential τ (Gyr) | 0.43 | 0.68 | 1.21 |
| metallicity (Z/Z☉) | 0.17 | 0.89 | 1.87 |
| A_V (mag) | 0.34 | 0.82 | 1.55 |
| z_phot | 0.674 | 0.675 | 0.676 |

→ M★ = 2.15 × 10¹¹ M☉ at the headline aperture. z_phot is foreced to agree with the
line-fit z = 0.675 given the strong prior we implemented.

### Per-band residuals — flagged for transparency

| Filter | AB obs | AB model (p50) | Residual |
|---|---|---|---|
| F200LP | 22.613 | 22.563 | +4.7% |
| F140W | 19.134 | 19.264 | **−11.3%** |
| F150W2 | 18.942 | 19.033 | −8.0% |
| F322W2 | 18.604 | 18.505 | +9.6% |

The model underpredicts F140W and F150W2 (the rest-frame near-IR bands)
by ~8–11% and overpredicts the bluest and reddest bands by ~5–10%.
This is the classic "redder NIR slope than a single-τ SFH prefers"
tension — not catastrophic but flags that a more flexible SFH (delayed-τ
or non-parametric) might tighten the fit if revisited. Effect on M★ is
at the ±0.07–0.15 dex level already quoted.

### Bagpipes plotting convention (a recurring gotcha)

`model_galaxy.wavelengths` is **rest-frame** Angstroms (multiply by
(1+z) to plot in observed frame). `spectrum_full` flux is plotted
**directly** — no (1+z) correction. `model_photometry` is already in
observed-frame units. Confirmed against Bagpipes source at
`plot_spectrum_posterior.py:69`; see `feedback_bagpipes_sed.md`.

## II.6 Sersic-total photometry cross-check (**A5**)

Single-component 2D Sersic profiles fit per band (notebook 08;
`scripts/measure_Re.py` shares the same fitting framework with the
2026-05-11 bound-fix applied — n ∈ [1.0, 6.0], ellip ∈ [0.0, 0.6], 3-init
grid). Extrapolated to total flux, then re-fed through Bagpipes with
identical priors:

| Filter | AB aperture | AB Sersic-total | Δmag | Sersic n | r_eff (″) |
|---|---|---|---|---|---|
| F200LP | 22.613 | 20.672 | −1.94 | 1.42 | 2.49 |
| F140W | 19.134 | 18.282 | −0.85 | 1.54 | 1.88 |
| F150W2 | 18.942 | 18.154 | −0.79 | 1.40 | 1.72 |
| F322W2 | 18.604 | 17.633 | −0.97 | 1.97 | 2.03 |

| Quantity | Aperture (paper) | Sersic-total | Δ |
|---|---|---|---|
| log(M★/M☉) p50 | 11.33 +0.07/−0.09 | 11.40 +0.11/−0.15 | +0.065 dex (+16% in M★) |

The aperture under-counts outer-profile flux by 0.8–1.9 mag depending
on filter, but at the integrated-mass level the Sersic correction is
+0.07 dex — within the quoted aperture-photometry uncertainties. The
paper cites the aperture value as the headline; Sersic-total is the
sensitivity floor.

Caches: `results/bagpipes_sed_results.npz` (aperture, 500-sample
posterior), `results/bagpipes_sersic_refit.npz` (Sersic-total),
`results/sersic_total_photometry.npz` (Sersic fit parameters).

## II.7 Photometric systematic budget

| Component | Δ log(M★/M☉) | Notes |
|---|---|---|
| Bagpipes posterior (aperture, p16–p84) | **+0.07 / −0.09** | dominant; ~17% / 19% in M★ |
| Aperture vs Sersic-total | +0.065 dex | sub-dominant; sensitivity floor at +16% in M★ |
| SED residuals (NIR tension) | folded into posterior | flexible SFH might tighten ~0.05 dex; not pursued |
| Drizzle noise correlation | included in 10% per-band σ_flux | conservative |
| Centering | sub-percent (HST-mean robust) | inherited from I.3/F1 |

The aperture posterior **±0.08 dex** is the headline uncertainty. The
Sersic-total cross-check tells us the aperture is biased low by ~0.07 dex
relative to a single-Sersic total — at the headline-uncertainty level,
so the conservative paper choice is to quote the aperture value and
note the Sersic-total as a +16% sensitivity check.

---

# Part III — Joint summary and outstanding items

## III.1 Final tabulated systematics

```
σ_e(<R_e)            = 269.62 ± 13.27 km/s  (paper headline, NEW _mtwdo_
                                              reduction, wide arc-masked window,
                                              + He I 3819 (M8) + M10 sky audit
                                              + M11 systematic re-derivation
                                              + M12 R_e-source fold-in,
                                              symmetric quadrature)
                     = 269.62 -13.45 / +13.10 km/s  (asymmetric form: stat
                                              -5.16/+4.13 from pooled
                                              bootstrap + sys ±12.43 each side;
                                              used in Fig 2 title + Fig 3
                                              error bar)
                       260.44 ± 18.13 km/s  (OLD reduction CLEAN, wide arc-masked,
                                              SECOND-reduction cross-check)
                       265.76 ± 17.87 km/s  (NEW reduction LEGACY un-cleaned, ref)
                       254.85 ± 17.87 km/s  (OLD reduction LEGACY un-cleaned, ref)
                       267.95 ± 30.10 km/s  (narrow Ca H+K+G, OLD cube X-check)
σ_e budget (WIDE,M12)= stat ±4.6 ⊕ Ishape ±2.27 ⊕ mask ±6.65
                       ⊕ frame ±5 ⊕ centering ±4 ⊕ window ±3.82
                       ⊕ reduction ±3.45 ⊕ R_e-source ±6.13
                       (M11 cube-matched re-derivation + M12 D7 fold-in;
                        sys ±12.43, total ±13.27)

R_e                  = 2.305″ = 16.23 kpc  (HST F140W+F200LP CoG mean;
                        FLAG A3c: raw CoG non-convergent, method sys ±0.08″)
log(M★/M☉)           = 11.16 +0.32 / −0.08  (principled IR-extended masking, 10%;
                        supersedes 11.33 aperture)
                       11.04 ± 0.14         (20% flux err; fill-in reach 11.46)

z_deflector          = 0.67564                (NIST air rest λ, nb04)
z_source             = 1.302                  (AGEL DR2; cross-checked via
                                               Fe II UV multiplet in KCRM
                                               blue arm, I.2.2)
```

## III.2 Methodology choices worth highlighting in the manuscript

1. **Wide arc-masked window with explicit source-emission catalog.** The
   single most consequential methodological choice. Halved the total
   σ_e uncertainty (±30 → ±18) and collapsed SPS-template spread by 6×.
   First documented catalog of source-emission contamination at a strong-
   lens deflector; reusable for any AGEL target with known z_source.
2. **Frame-aware ppxf (per-SPS vacuum/air).** Fixed a hidden bug that
   produced a −90 km/s FSPS V_sys offset throughout the nb07c era. The
   σ headline shifted by only +0.5 km/s, but the V_sys split-track
   collapsed 110 → 15 km/s and the §7 cross-check uncertainties halved.
3. **Sersic2D bound-fix.** Prevented the n=0.30 ellip=0.00 escape on
   the F200LP fit (a flat-disk pathology, not a physical solution).
   Removed a +20 km/s outlier in the I-shape sweep; the I-shape budget
   shrank from ±5.4 to ±1.5 km/s at the wide window.
4. **Quadrature-conservative budget with literature-standard estimator.**
   §6cum (cumulative I-weighted aperture ppxf) matches Cappellari+2006
   eq. 1 and Greene+20 ARA&A; three independent estimators (§7, §7b,
   arc-subtraction) agree within budget.
5. **Aperture-vs-Sersic-total mass sensitivity.** Quote the aperture
   M★ as the headline with the Sersic-total +16% as a stated floor —
   transparent about the aperture-photometry choice.

## III.3 Caveats explicitly retained in the paper

1. **No central σ.** σ_e at R<R_e/8 was dropped (aperture inside half
   the seeing FWHM, only 3 spaxels). We report σ_e(<R_e/2) as a
   gradient diagnostic but do not quote a central value.
2. **NIR SED residual.** Bagpipes underpredicts F140W and F150W2 by
   ~8–11%; a more flexible SFH would likely tighten this. Effect on M★
   already covered by quoted uncertainty.
3. **KCWI cube header mislabel.** Cite the multi-night combine with
   K409 Aug 30 UT 2025 as canonical; do *not* cite the header's
   "DATE-OBS = 2025-11-17, PROGNAME = U002, PROGPI = Jones" (last
   stacking pass only).
4. **F140W vs F200LP centroid offset = 0.36″.** Exceeds the 0.2″ sanity
   threshold; driven by the F200LP arc. HST-mean is robust because F140W
   is arc-clean.
5. **One telluric residual masked alongside the source-emission bands.**
   The def-rest 4525–4545 Å mask (obs 7593–7626 Å) is the leading edge
   of the O₂ atmospheric A-band, not source emission — confirmed
   2026-05-18 (`NOTES_4534A_spike_investigation_2026-05-18.md`). The
   mask is unchanged; only the attribution was corrected.

## III.4 What is *not* in this paper but is in the repo

- `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` — resolved σ
  map (per-spaxel + PowerBin) at R<1.5 R_e. Available as supplementary
  material if the referee asks for a kinematic map; the 17-spaxel
  S/N≥5 set shows the expected elliptical-bulge σ(R) gradient.
- `notebooks/04a_companion_redshifts.ipynb` — companion-galaxy redshifts
  in the same KCWI cube; outside scope but kept for the cluster-membership
  question.
- `notebooks/10_ETG_Mstar_radius_literature.ipynb` — comparison of
  AGEL0206 to local-ETG M★–R_e samples; supplementary context.

## III.5 Open items before submission

(See `project_roadmap.md` for the current task list.)

0. **Hδ targeted treatment (TODO 2026-05-26).** The new `_clean` production
   pipeline removes the default ppxf ±800 km/s mask around the Balmer
   absorption lines (Hδ 4101.74, Hγ 4340.47, Hβ 4861.33) because they are
   absorption-dominated, not emission, in this passive deflector. **Hδ
   specifically may still need a targeted clip** — its absorption-line
   shape is more template-sensitive than Hγ/Hβ (higher-order Balmer line,
   sensitive to template Balmer-decrement mismatch). To verify: re-run
   the local-MAD test (M5) on the no-Balmer-mask ppxf residuals and see
   whether Hδ shows up as a flagged region. Possible follow-ups: (a) keep
   all Balmer unmasked as in current production, (b) mask only Hδ at the
   ppxf-default ±800 km/s width, (c) mask Hδ at a narrower ±300 km/s
   width to keep most of its absorption signal while suppressing template
   mismatch. Flagged as inline comment in `scripts/bootstrap_ppxf.py:
   _determine_goodpixels_no_balmer`.
1. **Manuscript kinematics section** still needs rewriting around the
   wide arc-masked window headline.
2. **K409 PI** for the Aug 30 night needs to be confirmed (paper
   acknowledgements / observing-program citation).
3. **Per-night DIMM seeing** for the K409 Aug 30 night needs to be
   extracted from the observing log (currently using GUIDFWHM = 1.27″
   from the Nov 17 stacking-pass night as a conservative estimate).
4. ~~**Asymmetric posterior errors on Figure 3.**~~ **DONE 2026-05-19.**
   Both Figures 2 and 3 now use the asymmetric pooled posterior:
   σ_e(<R_e) = **254.85 −18.20 / +17.59 km/s** (total, stat from
   posterior + sys ±16.81 in quadrature each side). Stat alone:
   −7.0 / +5.2 km/s. Old symmetric ±17.87 lines kept commented in the
   figure cells (`fig2_wide`, `d55da320`) for traceability, with a
   dated update marker.
5. **R_e/8 caveat phrasing.** Add a one-sentence paper note that the
   aperture was dropped because it sits inside the seeing FWHM/2,
   referencing the seeing measurement (open item above).
6. **Blue-arm Fe II source-z cross-check artifacts.** The work
   described in I.2.2 was performed but the notebook/script and the
   measured z_source ± uncertainty are not yet committed in the repo.
   Once committed, replace the `[TBD]` placeholder in I.2.2 with the
   path + value.

---

## Document index

- **`TESTS_AND_DIAGNOSTICS.md`** — the master test catalog (every
  individual test A1–L3 with status, paths, and cache locations).
- **`NOTES_nb09_frame_fix_and_final_sigma_e_2026-04-28.md`** — long-form
  derivation of the vacuum/air fix and the four methodology audits.
- **`HANDOFF_wide_arcmask_complete_2026-05-13.md`** — wide arc-masked
  window adoption + paper Figure 2 update.
- **`PROJECT_BRIEF.md`** — self-contained designer/agent handoff (still
  cites narrow-window 268 ± 32; pre-wide-window).
- **`CLAUDE.md`** — repository conventions and current headline summary.
- **Memory index** at
  `~/.claude/projects/-Users-rosador-Documents-AGEL-AGEL0206-kinematics-and-photometry/memory/MEMORY.md`
  for compressed cross-session references (project_*, reference_*,
  feedback_*).

## Selected literature references

- Cappellari & Emsellem (2004) — ppxf
- Cappellari et al. (2006) — eq. 1 σ_e definition (cumulative I-weighted)
- Cappellari (2017, 2023) — ppxf revisions, high-z fwhm_gal_dict pattern
- Ciddor (1996) — vacuum-to-air conversion
- Kormendy & Ho (2013) — M•–σ scaling for ellipticals
- Greene et al. (2020 ARA&A) — σ_e aperture + M•–σ updates
- Gültekin et al. (2009) — discrete annular σ_e sum
- Carnall et al. (2018, 2019) — Bagpipes
- Bezanson et al. (2018); van der Wel et al. (2021) — LEGA-C, closest
  literature precedent for the wide-window rest 3500–5500 Å range
- AGEL DR2 (Carleton et al., paper in prep)
