# Tests & Diagnostics — AGEL0206 σ_e ApJL paper

**Last updated:** 2026-04-29
**Headline:** σ_e(<R_e) = **268 ± 32 km/s** (stat ±24 ⊕ I-shape ±13 ⊕ mask ±16 ⊕ frame ±5 ⊕ centering ±4)
**Method:** §6cum cumulative I-weighted ppxf at R<R_e=2.305" inside the F200LP-masked, frame-aware, SPS-pooled (FSPS+EMILES+XSL) bootstrap × polynomial-degree posterior

This document is the canonical index of every test, sweep, audit, and
sensitivity check run for the σ_e measurement. Each row points to the
script/notebook that runs it and the result file or note it produced.

---

## 0. Quick-reference test matrix

| # | Test | Where | Result | Status |
|---|---|---|---|---|
| **A. Foundational** | | | | |
| A1 | KCWI cube provenance | `NOTES_methodology_2026-04-27.md` | Multi-night Aug/Sept/Dec 2025; header mislabel only | ✅ |
| A2 | Systemic redshift via line fitting | `notebooks/04_redshift_verification.ipynb`, `scripts/redshift_verify.py` | z = 0.67564 (vs DR2 0.67511) | ✅ |
| A3 | R_e from masked curve-of-growth (F140W+F200LP) | `scripts/measure_Re.py`, `scripts/final_sigma_e.py:curve_of_growth` | R_e = 2.305" = 16.23 kpc | ✅ |
| A4 | HST-mean center via centroid_2dg | `scripts/final_sigma_e.py:find_center` | F140W & F200LP agree to 0.36" | ✅ |
| A5 | Bagpipes SED fitting (M★) | `notebooks/02_streamlined_Bagpipes_SED.ipynb`, `08_sersic_total_photometry.ipynb` | log M★ = 11.33 +0.07/−0.09 | ✅ |
| **B. Pipeline correctness audits** | | | | |
| B1 | Instrumental LSF audit (DISPSCAL=0.294) | `scripts/ifu_spectral_resolution.py` | FWHM=0.692 Å; σ_v_inst≈12-14 km/s | ✅ |
| B2 | σ_inst sensitivity (×0.5 to ×2.0 LSF) | `scripts/sigma_inst_sensitivity.py`, audit 3 | max \|Δσ\| = 0.83 km/s | ✅ |
| B3 | SPS template wavelength frames | `scripts/audit_ppxf_methodology.py`, NOTES Test 2 | FSPS=vacuum, EMILES=air, XSL=air (Ca K minimum + V_sys) | ✅ |
| B4 | V_sys air vs vacuum × 5 polynomial degrees | `scripts/audit_ppxf_methodology.py` audit 1 | ΔV consistent within ±2 km/s across degrees | ✅ |
| B5 | z × air-vac differential (obs vs rest application) | `scripts/audit_ppxf_methodology.py` audit 2 | 1.83 km/s differential — sub-budget | ✅ |
| B6 | fwhm_gal_dict frame consistency check | `scripts/audit_ppxf_methodology.py` audit 4 | dict in REST frame, matches Cappellari pattern | ✅ |
| **C. SPS templates & frame-fix** | | | | |
| C1 | Frame-aware ppxf: per-SPS native frame | `scripts/bootstrap_ppxf.py:SPS_NATIVE_FRAME` + `frame_galaxy='auto'` | V_sys split collapsed from ~110 to ~15 km/s | ✅ |
| C2 | End-to-end V_sys closure (frame swap) | NOTES Test 3, addendum 2026-04-29 | Frame fix matches diagnosis | ✅ |
| C3 | σ shift from frame fix | NOTES addendum, `error_budget()` | ≤5 km/s across SPS → carried as ±5 km/s budget | ✅ |
| C4 | 3-SPS pooled posterior | `scripts/final_sigma_e.py:pool_sps` | Per-SPS V_sys subtracted before pooling | ✅ |
| **D. Aperture & I(r)** | | | | |
| D1 | §6cum cumulative I-weighted ppxf | `scripts/final_sigma_e.py:run_aperture_sps` | σ_e(<R_e) = 267.82 km/s headline | ✅ |
| D2 | §6cum vs §7 (annular) cross-check | `notebooks/07c`, `~/.claude/.../reference_cumulative_vs_annular_sigma_e.md` | §7=255, §7b=271 — all <1σ from §6cum | ✅ |
| D3 | I(r) shape sweep (10 shapes, fixed mask) | `scripts/run_isource_shape_sweep.py`, `scripts/analyze_isource_shape_sweep.py` | Range = 12.9 km/s (excluding F200LP_Sersic2D outlier) → ±13 budget | ✅ |
| D4 | I(r) shape sweep refresh post-frame-fix at N=500 | NOTES addendum (d), 2026-04-29 | Frame fix has <0.5 km/s impact; ±13 budget validated | ✅ |
| D5 | Aperture choice: R<R_e/8, R<R_e/2, R<R_e | `scripts/final_sigma_e.py:APERTURE_FRACS` | R<R_e/8 dropped (inside seeing FWHM/2 = 0.64") | ✅ |
| D6 | Equal-N vs equal-width annular binning | `notebooks/07b` (5-bin), `notebooks/07c` (equal-N) | Equal-N inside R_safe=3R_e/4 chosen | ✅ |
| **E. Mask treatment** | | | | |
| E1 | F200LP arc mask reprojected to IFU grid | `scripts/final_sigma_e.py:load_setup` | 0.7% of all spaxels flagged (~38 inside R<R_e) | ✅ |
| E2 | Hard mask (w=0.0) headline | `scripts/final_sigma_e.py`, track A | σ_e = 267.82 ± 24 km/s | ✅ |
| E3 | No-mask sensitivity (w=1.0) | track B | σ_e = 252.82, Δ = −15.0 km/s | ✅ |
| E4 | Soft-mask single-point (w=0.5) | `scripts/soft_mask_track.py`, track C, NOTES addendum (b) | σ_e = 258.54, Δ = −9.29; super-linear | ✅ |
| E5 | 5-point mask-weight sweep w∈{0,0.25,0.5,0.75,1} | `scripts/soft_mask_track.py --weight ...`, `scripts/analyze_mask_weight_sweep.py`, NOTES addendum (c) | Per-step drops 5.0→4.3→3.5→2.2; quadratic c=+7.3 (concave-up); threshold-dominated | ✅ |
| E6 | Mask dilation (over-masking diagnostic) | `notebooks/07d_sigma_e_forceful_mask.ipynb`, `~/.claude/.../project_nb07d_overmasking_finding.md` | Dilation inflates σ via noise; do NOT dilate (negative finding) | ✅ |
| E7 | Arc-spectrum subtraction sibling | `notebooks/07e_sigma_e_arc_subtract.ipynb` | Matches §6cum within 0.1 km/s — residual arc dilution sub-dominant | ✅ |
| **F. Centering** | | | | |
| F1 | 5-center sweep (±0.4" perturbations) | `NOTES_centering_investigation_2026-04-27.md`, `scripts/regen_s6cum_nomask_diagnostic.py` | σ_e spread = 3.7 km/s → ±4 km/s budget | ✅ |
| F2 | HST F140W vs F200LP centroid offset | `scripts/final_sigma_e.py:find_center` | 0.36" — driven by F200LP arc; HST-mean robust | ✅ |
| **G. Bootstrap & polynomial degree** | | | | |
| G1 | Wild bootstrap (Rademacher × local residual) | `scripts/bootstrap_ppxf.py:run_bootstrap_single_degree` | 75-pix rolling window for σ_loc | ✅ |
| G2 | Parallel bootstrap (joblib, BLAS=1 per worker) | `scripts/bootstrap_ppxf_parallel.py` | ~6 min for 6 fits at N=500 | ✅ |
| G3 | Polynomial degree sweep (deg 15-29) | `notebooks/03c`, `scripts/build_nb09.py:§6.5` | σ flat within bootstrap envelope; polynomial saturated | ✅ |
| G4 | N=50 smoke vs N=500 production agreement | `scripts/final_sigma_e.py --n_bootstrap` | Within 1 km/s | ✅ |
| H | **Cross-checks (independent estimators)** | | | |
| H1 | §7 discrete Gültekin annular sum | `notebooks/07c §7` | 254.99 −24.2/+28.4 km/s (arc-filtered) | ✅ |
| H2 | §7b flat-σ extrapolation into outer annulus | `notebooks/07c §7b` | 271 −33/+35 km/s | ✅ |
| H3 | F200 mask sensitivity at N=500 | `results/annular_bootstrap_07c_nomask/` | Δ_mask = −16.4 km/s → ±16 budget | ✅ |
| **I. Earlier sigma-discrepancy work (Tests 1-21)** | | | | |
| I1 | NB01 vs NB05 σ reproduction | `notebooks/05x` Tests 1-4 | Confirmed σ depends on aperture | ✅ |
| I2 | Source mask + PSF-aware masking | `notebooks/05x` Tests 7-10 | F140W mask flags deflector core; F200LP mask is the right arc mask | ✅ |
| I3 | Threshold sweep + contamination weighting | `notebooks/05x` Tests 11-13 | Pre-Gültekin masking strategy work | ✅ |
| I4 | V_rms profile + bootstrap | `notebooks/05x` Tests 14-16 | Pre-Gültekin radial work | ✅ |
| I5 | Redshift sensitivity (z=0.67511 vs 0.67564) | `notebooks/05x` Test 17 | <1 km/s effect on σ | ✅ |
| I6 | Voronoi binning attempt | `notebooks/05x` Test 18 | Failed — abandoned | ✅ |
| I7 | PowerBin spatial binning | `notebooks/05x` Test 19, `scripts/run_powerbin_test19.py` | Rejected as binning scheme — outer-bin σ>800 km/s | ✅ |
| I8 | R_e four-way comparison | `notebooks/05x` Test 20, `scripts/measure_Re.py` | F140W+F200LP masked CoG mean = 2.305" headline | ✅ |
| I9 | Per-spaxel vs PowerBin rotation map | `notebooks/05x` Test 21 | No rotation on Re-scale; σ-dominated | ✅ |

---

## 1. Headline numbers

### σ_e at three apertures (cumulative I-weighted ppxf, frame-aware, masked)

| Aperture | σ_e [km/s] | Notes |
|---|---|---|
| R<R_e (=2.305") | **268 ± 32** (stat ±24) | **paper headline** — KH13/Greene+20 σ_e |
| R<R_e/2 (=1.152") | 224 ± 18 | gradient diagnostic |
| R<R_e/8 (=0.288") | dropped | inside seeing FWHM/2 = 0.64" |

### Error-budget composition (R<R_e)

| Component | Value | Source |
|---|---|---|
| Statistical (bootstrap pooled 1σ) | ±24 | nb09, this paper |
| I-shape spread (10 shapes excluding outlier) | ±13 | nb07c sweep, refreshed N=500 post-frame-fix |
| Mask on/off Δ | ±16 | nb07c, nb09 §5 |
| Frame fix (max σ shift across SPS) | ±5 | NOTES addendum + audit 1 |
| Centering (5-center sweep, ±0.4") | ±4 | `NOTES_centering_investigation_2026-04-27.md` |
| **Quadrature total** | **±32** | `scripts/final_sigma_e.py:error_budget` |

### Sensitivity values (for paper text)

- Soft-mask (w=0.5) at R<R_e: 258.5 km/s — Δ from headline = −9.3 km/s
- No-mask (w=1.0) at R<R_e: 252.8 km/s — Δ from headline = −15.0 km/s
- σ_e(w) is essentially linear with weak super-linearity (quadratic c=+7.3)

---

## 2. Detailed test catalog

This expands each row of the matrix above. Subsections reference scripts,
notebooks, and result files — open the linked file for the detailed code
and the inline plots/output.

### A. Foundational data and geometry

**A1 KCWI cube provenance** —
The cube `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` is the final
multi-night reduction (Aug 29 K409 + Sept 29 K409 + Dec 29 + Nov 17 2025).
The header is mislabeled (DATE-OBS/PROGNAME/PROGPI describe only the last
stacking pass). Settled in `NOTES_methodology_2026-04-27.md`. **Cite the
multi-night provenance, not the header**.

**A2 Systemic redshift** —
Independent line-by-line Gaussian fitting of Ca H+K, [OII], G-band, etc.
gives z_systemic = 0.67564 (NIST air wavelengths). DR2 catalog z=0.67511
(used by older notebooks 03/03b) is offset by ~95 km/s from the line-fit
value; nb09 uses 0.67564 throughout.
Notebook: `notebooks/04_redshift_verification.ipynb`
Module: `scripts/redshift_verify.py` (provides ABSORPTION_LINES /
EMISSION_LINES dicts, `gaussian_fit`, `multi_line_fit`, `oii_doublet_fit`)

**A3 Effective radius R_e = 2.305"** —
Mean of F140W and F200LP masked curve-of-growth values (2.168 and 2.441").
The masked CoG zeros mask-flagged HST pixels before the radial integral so
the arc/companion contributions don't bias the bright-end half-light point.
Computed in `scripts/measure_Re.py` and re-derived in
`scripts/final_sigma_e.py:curve_of_growth` for the nb09 pipeline.
Cross-check: IFU white-light proper-masked CoG gave 2.61" earlier — the
HST-based 2.305" is now the headline.

**A4 HST-mean center** —
`photutils.centroid_2dg` on F140W and F200LP cutouts (with their own
`_mask.fits` as the mask) gives sub-pixel centroids in each HST WCS. We
average in world coordinates, then propagate to IFU sub-pixel:
`scripts/final_sigma_e.py:find_center`. F140W vs F200LP centroid offset =
0.36" — flags the F200LP arc-affected detection but the HST-mean is robust
(verified by the 5-center sweep, F1).

**A5 Stellar mass** —
Bagpipes SED fitting on HST + JWST aperture photometry → log M★/M☉ =
11.33 +0.07/−0.09 (notebook 02). Cross-checked at the Sersic-total level
in `notebooks/08_sersic_total_photometry.ipynb`.

### B. Pipeline correctness audits (`scripts/audit_ppxf_methodology.py`)

Four orthogonal correctness audits run end-to-end on the headline aperture
spectrum to verify ppxf usage matches Cappellari (2017, 2023). Saved in
`results/ppxf_methodology_audit.npz`. Documented in
`NOTES_nb09_frame_fix_and_final_sigma_e_2026-04-28.md` ADDENDUM 2026-04-29.

**B1 Instrumental LSF** —
KCWI header DISPSCAL = 0.294 (instrumental σ in pixels) × CD3_3 = 1.000
Å/pix → FWHM_inst = 2.355 × 0.294 = 0.692 Å (constant in observed Å). At
the rest-frame fit window (3879–4476 Å), σ_v_inst ≈ 12–14 km/s. Galaxy
σ ≈ 270 km/s is resolved by ~20×. Diagnosed in
`scripts/ifu_spectral_resolution.py`.

**B2 σ_inst sensitivity** —
Sweep DISPSCAL × {0.5, 0.75, 1.0, 1.25, 1.5, 2.0}. Max |Δσ| = 0.83 km/s
across all SPS. Reason: ppxf's `sps_lib` clips fwhm_diff² = (FWHM_gal² −
FWHM_tem²).clip(0) (sps_util.py:169). When templates are intrinsically
broader (FSPS, EMILES at our band), no convolution is applied — σ becomes
insensitive to DISPSCAL. LSF is rock-solid.
Script: `scripts/sigma_inst_sensitivity.py`

**B3 SPS template frames** —
Located the Ca K minimum in each shipped template:
- FSPS spectra_fsps_9.0.npz: Ca K at 3934.86 Å → +1.20 Å vs air rest = vacuum
- EMILES spectra_emiles_9.0.npz: 3933.82 Å → +0.16 Å = air
- XSL spectra_xsl_9.0.npz: V_sys closes at +9 km/s with air-conversion → air
End-to-end V_sys closure (NOTES Test 3) confirms.

**B4 V_sys air vs vacuum × 5 polynomial degrees** —
Run ppxf at degrees {15, 18, 21, 24, 27} for each SPS in both frames.
ΔV(air − vac) ≈ +83 km/s for FSPS, ≈ −82 km/s for EMILES & XSL,
consistent within ±2 km/s across degrees. Frame identification is robust.

**B5 z × air-vac differential** —
Apply `vac_to_air` at observed wavelengths (KCWI 6500–7500 Å) → v_offset
= −82.7 km/s. Apply at rest wavelengths (3879–4476 Å) → −84.5 km/s.
Differential = 1.83 km/s — sub-budget. Either convention is defensible;
we follow Cappellari's pattern (apply at obs).

**B6 fwhm_gal_dict frame check** —
Confirmed `_prep_spectrum_for_ppxf` builds the dict as
`{"lam": lam_gal_rest, "fwhm": fwhm_gal_rest}` in REST frame, matching
`ppxf_example_high_redshift.py:99-101`. FWHM_obs = 0.692 Å → FWHM_rest =
0.413 Å at z=0.67564. `sps_util.py:167` interpolates this onto template
`lam_temp` (also rest), then `varsmooth` pre-broadens.

### C. SPS templates and frame-fix

**C1 Frame-aware ppxf** —
`scripts/bootstrap_ppxf.py` exposes `frame_galaxy='auto'` which reads
`SPS_NATIVE_FRAME = {fsps: vacuum, emiles: air, xsl: air}`. The galaxy
spectrum is prepared in the matching frame (`util.vac_to_air()` applied
as scalar median when needed, matching `ppxf_example_kinematics_sdss.py`
:127-134). Diagnosed Apr 2026 — pre-fix the pipeline always converted to
air, producing a −90 km/s FSPS V_sys offset.

**C2 V_sys closure** —
Same galaxy spectrum, swap only the frame. Pre-fix: FSPS −90, EMILES +3,
XSL +9. Post-fix (frame-aware): FSPS −7, EMILES +3, XSL +9. Split-track
collapsed from ~110 → ~15 km/s. NOTES Test 3.

**C3 σ shift from frame fix** —
Side-effect on σ: ≤5 km/s per SPS, ≤2 km/s on the 3-SPS pool. Carried as
±5 km/s in the error budget (ironically dominated by the per-SPS shifts
that pool out, but conservative).

**C4 3-SPS pooled posterior** —
`scripts/final_sigma_e.py:pool_sps` concatenates per-SPS bootstrap σ
samples after subtracting per-SPS V_sys medians (so V offsets don't
inflate the V posterior; σ posteriors are concatenated raw — no V_sys
correction needed for σ). The headline σ_e(<R_e) is the 16/50/84
percentile of this combined-SPS pool.

### D. Aperture and I(r) framework

**D1 §6cum cumulative I-weighted ppxf** —
For each spaxel inside R<R_max, build I-weight from the IFU 6500–7500 Å
white-light band, drop arc-mask-flagged spaxels, normalize and sum to a
single 1-D spectrum, fit ppxf. This is the headline path — matches
KH13 / Greene+20 / SAURON / ATLAS3D / MaNGA σ_e definition (Cappellari+2006
eq. 1). Cross-references in
`~/.claude/.../reference_cumulative_vs_annular_sigma_e.md`.
Code: `scripts/final_sigma_e.py:run_aperture_sps`.

**D2 §6cum vs §7 (annular)** —
§7 discrete Gültekin annular sum at R<R_safe gives 254.99 km/s; §7b flat-σ
extrapolation gives 271 km/s; both within 1σ of §6cum's 267.82 km/s.
§6cum is the headline because (a) it matches the literature definition,
(b) it preserves LOSVD line-shape information that §7's moment-pooling
discards, (c) bright-center I-weighting auto-suppresses arc contamination.

**D3 I(r) shape sweep (10 shapes)** —
Hold the spaxel selection fixed (F200-raw mask, R<R_e), vary only the
I-weight map. Shapes:
1. IFU_band — 6500-7500 Å (headline)
2. IFU_wl — full IFU white-light
3. F140W_raw — HST F140W reprojected
4. F200LP_raw — HST F200LP reprojected
5. F140W_arcmasked — F140W with arc zeroed
6. F200LP_arcmasked — F200LP with own arc-mask zeroed
7. F140W_1Dcog — annular-mean F140W → interp1d to spaxels
8. F200LP_1Dcog — annular-mean F200LP → interp1d
9. F140W_Sersic2D — full 2D Sersic fit
10. F200LP_Sersic2D — full 2D Sersic fit (poor fit, n=0.30 — outlier)

Range excluding F200LP_Sersic2D outlier = 12.9 km/s → ±13 km/s budget.
Script: `scripts/run_isource_shape_sweep.py`
Analysis: `scripts/analyze_isource_shape_sweep.py`
Caches: `results/annular_bootstrap_07c_ishape/{shape}_{sps}.npz`
Figure: `results/figures/nb07c_ishape_sweep.png`

**D4 I-shape sweep refresh post-frame-fix** —
Original sweep at N=250 predated the frame fix. Re-ran at N=500 with
`frame_galaxy='auto'` on 2026-04-29. Per-SPS σ shifts <0.5 km/s; spread
unchanged at 12.9 km/s; ±13 budget validated. Pre-frame caches preserved
at `results/annular_bootstrap_07c_ishape_oldframe/`. NOTES addendum (d).

**D5 Aperture choice** —
Three apertures considered: R<R_e/8 (=0.288"), R<R_e/2 (=1.152"), R<R_e
(=2.305"). R<R_e/8 dropped because it sits inside half the seeing FWHM
(=0.64") — only 3 spaxels, and the inner σ would underestimate due to
seeing-broadened sampling. R<R_e is the paper headline; R<R_e/2 is
quoted as a gradient diagnostic.

**D6 Equal-N vs equal-width annular binning** —
For the §7 cross-check only. Equal-N inside R_safe = 3R_e/4 + 1 outer
flagged bin (`notebooks/07c`) gives balanced bootstrap variance per bin
and isolates arc contamination in the outer bin. Equal-width annuli
(`notebooks/07b`, 5 bins) shifted σ_e by ~10 km/s — at the SPS-systematic
level (±24). Equal-N is the §7 paper choice.

### E. Mask treatment

**E1 F200LP arc mask reprojection** —
F200LP `_mask.fits` (HST 0.04"/pix, arc-tuned) is reprojected to the IFU
0.30"/spax grid via `scipy.ndimage.map_coordinates(order=0)` (nearest-
neighbor). 69 of 10000 spaxels (0.69%) flagged. Inside R<R_e: 38 of 184
spaxels in the aperture, ~21% — but only ~27% of the I-weight (arc spaxels
are mostly off-center).

**E2 Hard-mask (w=0.0) headline** — track A.
**E3 No-mask (w=1.0) sensitivity** — track B. Δ_mask = −15.0 km/s.
**E4 Soft-mask single point (w=0.5)** — track C. Δ_soft = −9.3 km/s.
**E5 5-point mask-weight sweep (w∈{0,0.25,0.5,0.75,1})** —
Per-step σ drops at R<R_e: 5.0 → 4.3 → 3.5 → 2.2 km/s. Quadratic
c = +7.3 (concave-up, super-linear). Threshold-dominated bias: a few
heavily-arc-contaminated spaxels carry most of the bias. First 25% of
arc weight captures 34% of the total mask sensitivity.
Scripts: `scripts/soft_mask_track.py`, `scripts/analyze_mask_weight_sweep.py`
Figure: `results/figures/nb09_mask_weight_sweep.png`
NOTES addendums (b) and (c).

**E6 Mask dilation** —
Tested whether expanding the F200 mask further reduces σ. **Negative
finding**: dilation inflates σ via noise — the dropped spaxels are S/N
contributors, not bias contributors. Memory:
`~/.claude/.../project_nb07d_overmasking_finding.md`. Notebook: `nb07d`.

**E7 Arc-spectrum subtraction sibling** —
Build an arc-only spectrum from F200-masked spaxels, fit α × arc + β ×
deflector model, subtract. Result matches §6cum within 0.1 km/s →
*residual* arc dilution is sub-dominant at production statistics. The
F200 mask already does the heavy lifting. Notebook: `nb07e`.

### F. Centering

**F1 5-center sweep** —
Tested HST_mean (paper choice), F140W only, F200LP only, IFU peak,
arc-suppressed IFU peak. σ_e spread = 3.7 km/s. Carried as ±4 km/s in
the error budget. Notes: `NOTES_centering_investigation_2026-04-27.md`.

**F2 F140W vs F200LP centroid** —
Δ = 0.36" — driven by the F200LP arc dragging that filter's centroid
toward the arc. F140W is the cleaner bulge centroid; HST-mean is the
robust compromise.

### G. Bootstrap and polynomial degree

**G1 Wild bootstrap** —
Hybrid Rademacher sign-flip × local-residual scaling (75-pix rolling
window). For each iteration: residual = (galaxy − bestfit) × sign,
scaled by local σ. Re-fit, record V/σ. Implemented in
`scripts/bootstrap_ppxf.py:run_bootstrap_single_degree`.

**G2 Parallel bootstrap** —
joblib + threadpoolctl pinning BLAS=1 per worker (otherwise BLAS thread
oversubscription pathology — saw 4-hour runtimes). 8 workers gives
~6 min for 6 fits at N=500. Script:
`scripts/bootstrap_ppxf_parallel.py:run_bootstrap_single_degree_parallel`.

**G3 Polynomial degree sweep** —
Degrees 15-29 (15 values) per fit. Diagnostic: σ vs degree should be
flat within bootstrap envelope; monotonic trend indicates the polynomial
is absorbing LOSVD signal. Plotted at `nb09 §6.5` and saved to
`results/figures/nb09_sigma_vs_degree.png`. Saturation confirmed at
this degree range — picked from earlier `notebooks/03c` analysis.

**G4 N=50 vs N=500** —
N=50 smoke and N=500 production agree on σ_e to within 1 km/s. N=50 is
~6 min, N=500 is ~50 min. Production headline is N=500; the smoke is for
quick rebuilds.

### H. Cross-checks

**H1 §7 discrete Gültekin** — `notebooks/07c §7`. Annular sum
σ_e²(<R) = Σ F_j (V_j² + σ_j²) / Σ F_j with arc-filter on outer annulus.
Result: 254.99 −24.2/+28.4 km/s. <1σ from §6cum.

**H2 §7b flat-σ extrapolation** — `notebooks/07c §7b`. Push the §7 sum
into the outer annulus assuming σ_outer = σ at the last clean annulus.
Result: 271 −33/+35 km/s. Bracket on the §7 result.

**H3 F200 mask sensitivity (independent)** — `results/annular_bootstrap_07c_nomask/`.
Same §6cum pipeline with mask off. σ_e(<R_e) = 250.96 km/s. Δ = −16.4
km/s. Independently confirms the nb09 mask-weight sweep at w=1.0.

### I. Earlier sigma-discrepancy work (`notebooks/05x`)

The 21-test diagnostic notebook resolving the σ ≈ 301 km/s (190-spaxel
box, NB01) vs σ ≈ 210 km/s (inner circular aperture, NB05) discrepancy.
**Resolution:** σ depends on aperture; the 190-spaxel box is contaminated
by the lensed arc. This work motivated the cumulative-I path (Tests 7-13
explored masking variants, Tests 14-16 went radial, Tests 18-19 tried
binning approaches). PowerBin (Test 19) was rejected. The headline
analysis moved to circular apertures + cumulative I-weighting (nb06 →
nb07c → nb09).

---

## 3. File index

### Documents (top-level)

| File | Purpose |
|---|---|
| `TESTS_AND_DIAGNOSTICS.md` | **This file** — master test index |
| `NOTES_nb09_frame_fix_and_final_sigma_e_2026-04-28.md` | Frame fix derivation + addenda (a)–(d) |
| `HANDOVER.md` | Apr 27 σ_e analysis snapshot (pre-frame-fix; superseded by nb09) |
| `HANDOFF_Ir_sweep_2026-04-28.md` | I(r) sweep design doc (consumed; sweep complete) |
| `NOTES_methodology_2026-04-27.md` | Cube provenance + earlier methodology |
| `NOTES_centering_investigation_2026-04-27.md` | 5-center sweep result + ±4 km/s budget |
| `CLAUDE.md` | Repo conventions, key results summary |
| `README.md` | High-level project description |

### Notebooks (paper-ready / current)

| Notebook | Purpose |
|---|---|
| `notebooks/09_final_sigma_e_paper.ipynb` | **Paper-ready σ_e** (headline 268 ± 32 km/s) |
| `notebooks/04_redshift_verification.ipynb` | z = 0.67564 line-fit |
| `notebooks/02_streamlined_Bagpipes_SED.ipynb` | log M★ = 11.33 |
| `notebooks/08_sersic_total_photometry.ipynb` | Sersic-total photometry cross-check |
| `notebooks/07c_sigma_e_equalN.ipynb` | §6cum + §7 equal-N annular cross-check |
| `notebooks/07e_sigma_e_arc_subtract.ipynb` | Arc-subtraction sibling cross-check |
| `notebooks/05x_sigma_discrepancy_diagnostic.ipynb` | Tests 1-21 (foundational diagnostics) |

### Notebooks (history / reference, not paper-driving)

`01_*`, `02_*`, `03_bootstrap_ppxf_errors`, `03b_*`, `03c_*`, `05_*`,
`06_*`, `07_*`, `07a_*`, `07b_*`, `07d_*` — kept for the history of the
analysis.

### Scripts (production)

| Script | Purpose |
|---|---|
| `scripts/final_sigma_e.py` | nb09 production driver — runs all 3 mask tracks at N=500 |
| `scripts/build_nb09.py` | Generates the paper notebook from the saved artifacts |
| `scripts/bootstrap_ppxf.py` | Frame-aware ppxf input prep + wild bootstrap |
| `scripts/bootstrap_ppxf_parallel.py` | joblib parallel bootstrap (BLAS=1 per worker) |
| `scripts/measure_Re.py` | R_e via masked CoG (multi-band, multi-strategy) |
| `scripts/redshift_verify.py` | Line-fitting redshift module |
| `scripts/ifu_spectral_resolution.py` | LSF audit |

### Scripts (sweeps and audits)

| Script | Purpose |
|---|---|
| `scripts/audit_ppxf_methodology.py` | 4 orthogonal correctness audits |
| `scripts/sigma_inst_sensitivity.py` | DISPSCAL × {0.5..2.0} sweep |
| `scripts/run_isource_shape_sweep.py` | 10-shape I(r) sweep |
| `scripts/analyze_isource_shape_sweep.py` | I-shape sweep analysis + figure |
| `scripts/soft_mask_track.py` | Mask-weight runs (CLI `--weight`) |
| `scripts/analyze_mask_weight_sweep.py` | 5-point sweep analysis + linearity fit |
| `scripts/regen_s6cum_nomask_diagnostic.py` | Centering investigation |
| `scripts/relock_nomask_Re_N500.py` | One-off N=500 relock for the no-mask track |

### Result files (key)

| File | Contents |
|---|---|
| `results/final_sigma_e_paper.npz` | nb09 headline + per-SPS + 3-track summaries |
| `results/final_sigma_e_paper/` | Per-(aperture, SPS, mask_weight) caches |
| `results/sigma_e_ishape_sweep.npz` | I-shape sweep summary |
| `results/sigma_e_ishape_sweep_oldframe.npz` | Pre-frame-fix I-shape (audit) |
| `results/annular_bootstrap_07c_ishape/` | Per-(shape, SPS) I-shape caches (N=500) |
| `results/annular_bootstrap_07c_ishape_oldframe/` | Pre-frame-fix I-shape caches (N=250) |
| `results/annular_bootstrap_07c_nomask/` | §6cum nomask caches |
| `results/mask_weight_sweep.npz` | 5-point mask-weight sweep summary |
| `results/ppxf_methodology_audit.npz` | 4-audit results |
| `results/figures/nb09_*.png` | nb09 paper figures |
| `results/figures/nb07c_ishape_sweep.png` | I-shape sweep figure |

---

## 4. Reproduction recipes

```bash
conda activate ISMgas

# === Headline pipeline (~50 min at N=500) ===
python scripts/final_sigma_e.py --n_bootstrap 500
python scripts/build_nb09.py
jupyter nbconvert --to notebook --execute --inplace notebooks/09_final_sigma_e_paper.ipynb

# === Methodology audits (~3 min) ===
python scripts/audit_ppxf_methodology.py
python scripts/sigma_inst_sensitivity.py

# === I-shape sweep (~50 min at N=500) ===
python scripts/run_isource_shape_sweep.py --n-bootstrap 500 --n-jobs 8
python scripts/analyze_isource_shape_sweep.py

# === Mask-weight sweep (~25 min total) ===
# (w=0 and w=1 are produced by final_sigma_e.py; only need 0.25/0.5/0.75)
python scripts/soft_mask_track.py --weight 0.25
python scripts/soft_mask_track.py --weight 0.50
python scripts/soft_mask_track.py --weight 0.75
python scripts/analyze_mask_weight_sweep.py

# === Smoke run (~6 min, useful for rebuilds) ===
python scripts/final_sigma_e.py --n_bootstrap 50
```

Caches are keyed by N_bootstrap and mask_weight in the filename — re-runs
at the same N skip cached fits unless `--force` is passed.

---

## 5. What's intentionally NOT done

- **Re-run nb09 at higher N** (e.g. N=2000): the pooled posterior is already
  smooth; SPS systematics dominate the residual statistical variance.
- **Re-fit Sersic2D F200LP**: the n=0.30 fit is non-physical due to UV-arc
  domination; the F140W Sersic fit is the trustworthy one. F200LP_Sersic2D
  is the I-shape outlier *by design* and is excluded from the budget.
- **Voronoi binning** — Test 18 / nb05x — abandoned; not pursued at N=500.
- **PowerBin spatial binning** — Test 19 / nb05x — rejected (outer-bin
  σ > 800 km/s); not pursued.
- **Per-spaxel rotation map at production statistics** — no rotation
  detected at the seeing limit; not paper-driving.
- **Full nb07a Sersic-I revisit at the post-frame-fix headline** —
  superseded by the I-shape sweep (D3/D4) which covers the Sersic2D shape.

---

*See `~/.claude/projects/.../memory/MEMORY.md` for compressed cross-session
references, including:*
- *`reference_cumulative_vs_annular_sigma_e.md` — §6cum vs §7 design choice*
- *`reference_sps_systematic.md` — per-SPS V_sys + frame-fix details*
- *`reference_ppxf_vacair_handling.md` — vac/air conversion methodology*
- *`reference_kcwi_data_properties.md` — KCWI cube + LSF + paper boilerplate*
- *`project_nb09_final_sigma_e.md` — current headline status*
