# Tests & Diagnostics — AGEL0206 σ_e ApJL paper

**Last updated:** 2026-06-11 (M13 "best mask throughout" — best-mask R_e=2.097" cascade)

**Headline (NEW `_mtwdo_` + bad-pixel mask + Balmer-unmasked + He I 3819 mask + M10 sky audit + M11 systematic re-derivation + M13 best-mask R_e=2.097" aperture, wide arc-masked window):**
σ_e(<R_e) = **267 ± 13 km/s** = 267.31 ± 12.79 km/s
(stat ±4.6 ⊕ I-shape ±2.27 ⊕ mask ±6.65 ⊕ frame ±5 ⊕ centering ±4 ⊕ window ±3.82 ⊕ reduction ±3.45 **⊕ R_e-source ±5.13**)
Asymmetric form: 267.31 −12.98 / +12.61 km/s.
**All headline numbers are emitted deterministically by `scripts/paper_values.py` → `results/PAPER_VALUES.json` (single source of truth; do not hand-edit).**
Source (headline aperture): `results/run_wide_sigma_e/resys_best_mean/` (R_e=2.097"), via `scripts/run_sigma_e_Re_grid.py --n_bootstrap 500`. **M13:** adopting the best-mask R_e=2.097" (was 2.305") moved σ_e −2.3 km/s along σ(R) and re-derived R_e-source ±6.13→±5.13 (best-mask CoG family).

**Second-reduction cross-check (OLD cube + bad-pixel mask + Balmer-unmasked + He I + M10 sky masks):**
σ_e(<R_e) = 262.72 ± 17.99 km/s (asymmetric −18.10 / +17.88). The 6.90 km/s Δ
between the two CLEANED + He-I + M10-sky-masked reductions sets the refined ±3.45 km/s
reduction-pass systematic (was ±3.86 from clean+He-I-only cubes; M10 sky masks
tighten the inter-reduction gap further).

**Narrow-window cross-check (Ca H+K + G, OLD cube):**
σ_e(<R_e) = **268 ± 30 km/s** = 267.95 ± 30.10 km/s <!-- pv-skip: narrow-window cross-check, not the live headline -->
(stat ±24 ⊕ I-shape ±3.7 ⊕ mask ±7.5 ⊕ frame ±5 ⊕ centering ±4 ⊕ window ±15)

**Three-way consistency:** the headline (271.87, new clean + He I) and narrow-window cross-check (267.95, old narrow) are
within 4 km/s — both inside each other's 1σ stat error. The pre-He-I new-cube wide-window value
(268.98) was within 1 km/s of the narrow window; adding the He I mask moves the wide headline
~3 km/s redder of the narrow value, but still well within consistency.

**Method (wide, headline):** Single ppxf fit on the I-weighted aperture
spectrum at R<R_e=2.097" over rest 3800–5400 Å with explicit z=1.302
source-emission masks (MgII λ2800, [OII] λ3727, **He I 3819** (M8),
[NeIII]/Mg b cluster) + bad-pixel mask + Balmer kept in fit, F200LP
spatial arc mask, frame-aware (per-SPS vac/air), SPS-pooled (FSPS+EMILES+XSL)
wild bootstrap × polynomial-degree posterior.

**Source-emission verification figure (appendix candidate):** integrated arc
spectrum at `AGEL_0206_ApJL_Figures/AGEL0206_arc_source_spectrum.pdf` — direct
extraction from the 69 arc-flagged spaxels (deflector continuum subtracted)
showing [OII] 3727 + H9 + [NeIII] + He I 3819 in source-rest frame and the
same lines mapped to def-rest (the headline fit window).

**Method (narrow, cross-check):** Same machinery but on obs-frame
6500–7500 Å (rest 3879–4476 Å, anchored on Ca H+K and G-band only) —
the §6cum cumulative I-weighted ppxf path retained from nb07c.

This document is the canonical index of every test, sweep, audit, and
sensitivity check run for the σ_e measurement. Each row points to the
script/notebook that runs it and the result file or note it produced.

---

## 0. Quick-reference test matrix

| # | Test | Where | Result | Status |
|---|---|---|---|---|
| **A. Foundational** | | | | |
| A1 | KCWI cube provenance — Keck observing logs cross-checked | `NOTES_methodology_2026-04-27.md`, `reference_kcwi_data_properties.md` (updated 2026-05-26) | Verified contributing nights: K409 UT 2025 Aug 30 (12 RED + 4 BLUE DESJ0206 frames) + U002 PI Jones UT 2025 Nov 17 (20 RED + 5 BLUE DESJ0206 frames). Dec 29 2024 RED frames flagged pending FITS-header verification. Sept 29 K409 has zero DESJ0206 entries — earlier note was wrong. | ✓ |
| A2 | Systemic redshift via line fitting | `notebooks/04_redshift_verification.ipynb`, `scripts/redshift_verify.py` | z = 0.67564 ± 0.00033 (6 abs lines; vs DR2 0.67511). Source z=1.30263 ([O II] doublet). DRAFTING §2.2.2 full table | ✓ |
| A2b | Air↔vacuum audit of redshift features (2026-06-11) | DRAFTING §2.2.2 audit box; `scripts/redshift_verify.py`, nb04 | **PASS.** All rest λ NIST air <0.01 Å (G-band=CH bandhead, flagged). Red cube vacuum (`WAVE`)→air via Ciddor96 (same as σ_e); z=λ_obs,air/λ_rest,air−1 consistent → unbiased (else +90 km/s). **Blue cube is air (`AWAV`)** ≠ red; nb04 blue cells double-converted (latent bug, avoided in new scripts). | ✓ |
| A2c | Companion-galaxy redshifts (group membership) | `notebooks/04a_companion_redshifts.ipynb`, `scripts/redshift_verify.bootstrap_redshift` | NE (z=0.6758) + SW (z=0.6759) ellipticals, ~4.5″≈31 kpc, both at deflector z (\|Δcz\|≲60 km/s) → compact group/triple at z≈0.676. These are the masked field ellipticals (confirmed not interlopers). | ✓ |
| A2d | Source-arc rest-UV resonance lines | `scripts/source_uv_redshift.py`, `results/source_uv_redshift.npz` (2026-06-11) | **Mg II λλ2796,2803 (red arm): 16.7σ/18.4σ**, z=1.30113, −195 km/s vs systemic [O II] — reported as a velocity offset, **NOT claimed as outflow** (emission-infill/lensing/profile degenerate). **Fe II λλ2344,2382 (blue arm): 7.9σ/3.0σ** but +470 km/s cross-arm offset = blue-arm calib (TODO). Deflector rest-UV non-detected (too faint). | ✓ |
| A3 | R_e from best-mask curve-of-growth (F140W+F200LP) | `scripts/final_sigma_e.py:curve_of_growth`, `scripts/Re_bestmask_reconciliation.py` | **R_e = 2.097" = 14.76 kpc** (best-mask CoG, 2026-06-11 "best mask throughout"; was expert-mask 2.305"). Method sys ±0.100" (raw 2.097/sky-sub 1.922/Sérsic 1.897) | ✓ |
| A3d | R_e/CoG validation vs photutils | `scripts/validate_Re_photutils.py`, `results/validate_Re_photutils.npz` | **Our azimuthal-fill CoG reproduced by photutils `RadialProfile` to ±0.002"** (best-mask 2.097→2.077); unmasked matches `CurveOfGrowth` (2D aperture sum) ±0.004". Naive masked aperture-sum biases +0.25–0.41" high → fill is correct. centroid_2dg robust; centroid_quadratic fails on extended light | ✓ |
| A4 | HST-mean center via centroid_2dg | `scripts/final_sigma_e.py:find_center` | F140W & F200LP agree to 0.36"; centroid_2dg vs centroid_com 0.24–0.26" (A3d) | ✓ |
| A5 | Bagpipes SED fitting (M★) — aperture vs Sersic-total | `notebooks/02_streamlined_Bagpipes_SED.ipynb`, `08_sersic_total_photometry.ipynb` | log M★ = 11.33 +0.07/−0.09 (old expert-aperture); 11.40 +0.11/−0.15 (Sersic-total). **Superseded by N-series principled masking → headline 11.16 +0.31(sys)** | ✓ |
| A5b | Aperture-corrected photometry validation vs established codes | `scripts/validate_apcorr_established.py`, `results/validate_apcorr_established.npz` (2026-06-11) | **Full apcorr chain validated** (petrofit/statmorph absent → photutils + scipy incomplete-gamma refs). **b_n** Ciotti&Bertin vs exact gammaincinv ≤0.06%; **Sérsic total-flux formula** (Graham&Driver05 eq.4) vs render-to-20R_e ≤0.03%; **enclosed-light** rendered model vs analytic gammainc(2n,b_n·2^{1/n}) Δ≤0.0007; **F_raw** hard-edge `in_aperture` vs photutils `EllipticalAperture(exact)` ≤0.09%; **cutout truncation** of sum(model_full) ≤0.19% (≤0.002 mag). No production bug. | ✓ |
| A5c | Single-Sérsic ellipticity railing bug + multi-start fix | `scripts/sersic_total_photometry.fit_sersic2d`, `photometry_systematics.fit_sersic` (2026-06-11) | Single low-ellip start railed to a **spurious circular minimum** (b/a=1 all bands). **Multi-start** (moment seed + PA grid + circular fallback) recovers **F140W b/a=0.86, F150W2 0.80, F322W2 0.85, F200LP≈1 (faint)**. Cascade benign: masks <0.05% px, R_e 2.097→2.093 (kept 2.097), **M★ headline 11.47 unchanged**; Sérsic-only budget ±0.17→±0.19. | ✓ |
| A5d | Sérsic fitter validation vs published reference (synthetic) | `scripts/validate_sersic_fitter_synthetic.py`, `results/validate_sersic_fitter_synthetic.npz` | Injection-recovery vs astropy `Sersic2D` (peer-reviewed ref; GALFIT/imfit not installed). Validated for the deflector regime n≈1.2–1.6: strong ellip (b/a≤0.80) clean at all n; mild (b/a~0.85) clean at n≥2, borderline at n=1 (detectability limit). Realistic **b/a uncertainty ~±0.06** (bootstrap too tight). | ✓ |
| A5e | M★ low-mass tail = outshining (not the estimator) | `/tmp` analysis on `aperture_matched_photometry.npz` + Bagpipes posterior (2026-06-11) | Longer low tail present in **all 5 estimators incl. empirical raw** → not the total-Sérsic default. Age–dust–M/L degeneracy: low-M★ tail corr with age +0.90, dust −0.58, sSFR −0.83. Keep exponential-SFH prior, report as-is. | ✓ |
| A5f | Sérsic appendix parameter table (4 bands) | `scripts/sersic_parameter_table.py`, `results/sersic_parameter_table.{md,tex,npz}` | Per-band r_eff, n, b/a, PA, μ_e, m_tot; b/a/PA/n errors floored by A5d scatter; F200LP flagged circular. Deflector mildly elliptical b/a≈0.80–0.86, PA≈4–11° E of N, n≈1.2–1.6. DRAFTING §3.2. | ✓ |
| **B. Pipeline correctness audits** | | | | |
| B1 | Instrumental LSF audit (DISPSCAL=0.294) | `scripts/ifu_spectral_resolution.py` | FWHM=0.692 Å; σ_v_inst≈12-14 km/s | ✓ |
| B2 | σ_inst sensitivity (×0.5 to ×2.0 LSF) | `scripts/sigma_inst_sensitivity.py`, audit 3 | max \|Δσ\| = 0.83 km/s | ✓ |
| B3 | SPS template wavelength frames | `scripts/audit_ppxf_methodology.py`, NOTES Test 2 | FSPS=vacuum, EMILES=air, XSL=air (Ca K minimum + V_sys) | ✓ |
| B4 | V_sys air vs vacuum × 5 polynomial degrees | `scripts/audit_ppxf_methodology.py` audit 1 | ΔV consistent within ±2 km/s across degrees | ✓ |
| B5 | z × air-vac differential (obs vs rest application) | `scripts/audit_ppxf_methodology.py` audit 2 | 1.83 km/s differential — sub-budget | ✓ |
| B6 | fwhm_gal_dict frame consistency check | `scripts/audit_ppxf_methodology.py` audit 4 | dict in REST frame, matches Cappellari pattern | ✓ |
| **C. SPS templates & frame-fix** | | | | |
| C1 | Frame-aware ppxf: per-SPS native frame | `scripts/bootstrap_ppxf.py:SPS_NATIVE_FRAME` + `frame_galaxy='auto'` | V_sys split collapsed from ~110 to ~15 km/s | ✓ |
| C2 | End-to-end V_sys closure (frame swap) | NOTES Test 3, addendum 2026-04-29 | Frame fix matches diagnosis | ✓ |
| C3 | σ shift from frame fix | NOTES addendum, `error_budget()` | ≤5 km/s across SPS → carried as ±5 km/s budget | ✓ |
| C4 | 3-SPS pooled posterior | `scripts/final_sigma_e.py:pool_sps` | Per-SPS V_sys subtracted before pooling | ✓ |
| **D. Aperture & I(r)** | | | | |
| D1 | §6cum cumulative I-weighted ppxf | `scripts/final_sigma_e.py:run_aperture_sps` | σ_e(<R_e) = 267.82 km/s headline | ✓ |
| D2 | §6cum vs §7 (annular) cross-check | `notebooks/07c`, `~/.claude/.../reference_cumulative_vs_annular_sigma_e.md` | §7=255, §7b=271 — all <1σ from §6cum | ✓ |
| D3 | I(r) shape sweep (10 shapes, fixed mask) | `scripts/run_isource_shape_sweep.py`, `scripts/analyze_isource_shape_sweep.py` | Range = 12.9 km/s (excluding F200LP_Sersic2D outlier) → ±13 budget | ✓ |
| D4 | I(r) shape sweep refresh post-frame-fix at N=500 | NOTES addendum (d), 2026-04-29 | Frame fix has <0.5 km/s impact; ±13 budget validated | ✓ |
| D5 | Aperture choice: R<R_e/8, R<R_e/2, R<R_e | `scripts/final_sigma_e.py:APERTURE_FRACS` | R<R_e/8 dropped (inside seeing FWHM/2 = 0.64") | ✓ |
| D6 | Equal-N vs equal-width annular binning | `notebooks/07b` (5-bin), `notebooks/07c` (equal-N) | Equal-N inside R_safe=3R_e/4 chosen | ✓ |
| D7 | R_e source systematic — best-mask grid (was 4-estimator D7) | **`scripts/run_sigma_e_Re_grid.py`, `results/sigma_e_Re_grid_N500.npz` (2026-06-11)**; superseded `run_sigma_e_Re_systematic_wide.py` | **7-point σ_e-vs-R_e grid bracketing the 2.097" best-mask headline (N=500):** best_F140W 1.912"=259.73 / best_mean 2.097"=**267.31** / best_F200LP 2.281"=269.99 / exp 2.168/2.305/2.441 = 267.44/269.62/272.44 / CaHK+G 2.902"=279.69. **ADOPTED = best-mask CoG light family ±5.13** (1.912/2.097/2.281); CaHK+G & full grid (±9.98) are cross-checks, NOT folded. Was ±6.13 at the old expert headline. | ✓ folded (best-mask) |
| **E. Mask treatment** | | | | |
| E1 | F200LP arc mask reprojected to IFU grid | `scripts/final_sigma_e.py:load_setup` | 0.7% of all spaxels flagged (~38 inside R<R_e) | ✓ |
| E2 | Hard mask (w=0.0) headline | `scripts/final_sigma_e.py`, track A | σ_e = 267.82 ± 24 km/s | ✓ |
| E3 | No-mask sensitivity (w=1.0) | track B | σ_e = 252.82, Δ = −15.0 km/s | ✓ |
| E4 | Soft-mask single-point (w=0.5) | `scripts/soft_mask_track.py`, track C, NOTES addendum (b) | σ_e = 258.54, Δ = −9.29; super-linear | ✓ |
| E5 | 5-point mask-weight sweep w∈{0,0.25,0.5,0.75,1} | `scripts/soft_mask_track.py --weight ...`, `scripts/analyze_mask_weight_sweep.py`, NOTES addendum (c) | Per-step drops 5.0→4.3→3.5→2.2; quadratic c=+7.3 (concave-up); threshold-dominated | ✓ |
| E6 | Mask dilation (over-masking diagnostic) | `notebooks/07d_sigma_e_forceful_mask.ipynb`, `~/.claude/.../project_nb07d_overmasking_finding.md` | Dilation inflates σ via noise; do NOT dilate (negative finding) | ✓ |
| E7 | Arc-spectrum subtraction sibling | `notebooks/07e_sigma_e_arc_subtract.ipynb` | Matches §6cum within 0.1 km/s — residual arc dilution sub-dominant | ✓ |
| E8 | Arc-as-sky mechanism test (no-mask + arc-sky template) | `scripts/run_07f_arc_sky.py`, `notebooks/07f_arc_sky_subtract.ipynb` | σ_D=274.6 vs σ_A=267.8 / σ_B=252.8; recovery=145% — dilution explains the full no-mask Δ | ✓ |
| **F. Centering** | | | | |
| F1 | 5-center sweep (±0.4" perturbations) | `NOTES_centering_investigation_2026-04-27.md`, `scripts/regen_s6cum_nomask_diagnostic.py` | σ_e spread = 3.7 km/s → ±4 km/s budget | ✓ |
| F2 | HST F140W vs F200LP centroid offset | `scripts/final_sigma_e.py:find_center` | 0.36" — driven by F200LP arc; HST-mean robust | ✓ |
| **G. Bootstrap & polynomial degree** | | | | |
| G1 | Wild bootstrap (Rademacher × local residual) | `scripts/bootstrap_ppxf.py:run_bootstrap_single_degree` | 75-pix rolling window for σ_loc | ✓ |
| G2 | Parallel bootstrap (joblib, BLAS=1 per worker) | `scripts/bootstrap_ppxf_parallel.py` | ~6 min for 6 fits at N=500 | ✓ |
| G3 | Polynomial degree sweep (deg 15-29) | `notebooks/03c`, `scripts/build_nb09.py:§6.5` | σ flat within bootstrap envelope; polynomial saturated | ✓ |
| G4 | N=50 smoke vs N=500 production agreement | `scripts/final_sigma_e.py --n_bootstrap` | Within 1 km/s | ✓ |
| H | **Cross-checks (independent estimators)** | | | |
| H1 | §7 discrete Gültekin annular sum (frame-aware, refresh 2026-05-01) | `scripts/refresh_07c_gultekin.py`, `notebooks/07c §7` | **256.17 −13.0/+12.7 km/s** (post-fix; was 254.99 −24/+28 pre-fix) | ✓ |
| H2 | §7b flat-σ extrapolation (frame-aware, refresh 2026-05-01) | `scripts/refresh_07c_gultekin.py`, `notebooks/07c §7b` | **274.37 −16.2/+17.4 km/s** (post-fix; was 271 −33/+35 pre-fix) | ✓ |
| H3 | F200 mask sensitivity at N=500 | `results/annular_bootstrap_07c_nomask/` | Δ_mask = −16.4 km/s → ±16 budget | ✓ |
| **I. Earlier sigma-discrepancy work (Tests 1-21)** | | | | |
| I1 | NB01 vs NB05 σ reproduction | `notebooks/05x` Tests 1-4 | Confirmed σ depends on aperture | ✓ |
| I2 | Source mask + PSF-aware masking | `notebooks/05x` Tests 7-10 | F140W mask flags deflector core; F200LP mask is the right arc mask | ✓ |
| I3 | Threshold sweep + contamination weighting | `notebooks/05x` Tests 11-13 | Pre-Gültekin masking strategy work | ✓ |
| I4 | V_rms profile + bootstrap | `notebooks/05x` Tests 14-16 | Pre-Gültekin radial work | ✓ |
| I5 | Redshift sensitivity (z=0.67511 vs 0.67564) | `notebooks/05x` Test 17 | <1 km/s effect on σ | ✓ |
| I6 | Voronoi binning attempt | `notebooks/05x` Test 18 | Failed — abandoned | ✓ |
| I7 | PowerBin spatial binning | `notebooks/05x` Test 19, `scripts/run_powerbin_test19.py` | Rejected as binning scheme — outer-bin σ>800 km/s | ✓ |
| I8 | R_e four-way comparison | `notebooks/05x` Test 20, `scripts/measure_Re.py` | F140W+F200LP masked CoG mean = 2.305" (expert mask); **superseded by best-mask 2.097" (A3/A3d, 2026-06-11)** | ✓ |
| I9 | Per-spaxel vs PowerBin rotation map | `notebooks/05x` Test 21 | No rotation on Re-scale; σ-dominated | ✓ |
| **J. Wide arc-masked window (2026-05-08 → 2026-05-13)** | | | | |
| J1 | Wavelength-window sweep (15 windows, narrow → wide) | `notebooks/09a_wavelength_window_sweep.ipynb`, `scripts/run_window_sweep.py` | wR3800_5400 wins on stat precision; XSL template floor at def-rest 3675 Å sets blue edge | ✓ |
| J2 | Discovery of z=1.302 source-emission contamination | `notebooks/09a §3a`, `notebooks/09b §8` | Bimodal ppxf posterior in wide window (σ≈250 vs 100–150) caused by source emission lines mapped to deflector rest frame | ✓ |
| J3 | Spectral-mask catalog: 3 source-emission + 1 telluric, def-rest | `scripts/run_window_sweep.py:ARC_MASKS_REST` | Four bands (def-rest 3835-3855 = Mg II z=1.302; 4525-4545 = O₂ A-band telluric edge at obs 7593-7626 Å; 5115-5135 = [O II] z=1.302; 5260-5340 = [Ne III] z=1.302); 212 of 2374 pixels (8.9%) masked. Telluric ID via `NOTES_4534A_spike_investigation_2026-05-18.md` | ✓ |
| J4 | Arc-mask + ppxf clean=True diagnostic | `notebooks/09b §8` | Bimodality collapses; σ vs polynomial degree FLAT across deg 15-29; clean=True drops 0 pixels | ✓ |
| J5 | Per-SPS posterior collapse at wide arc-masked | `notebooks/09d §1.3` | SPS spread (FSPS / EMILES / XSL): 26.0 km/s (narrow) → 4.2 km/s (wide arc-masked); SPS systematic essentially disappears | ✓ |
| J6 | N=500 production at wide arc-masked window | `results/nb09a_wavelength_sweep/wR3800_5400_arcmask_*_N500.npz` | σ_e = 254.9 −7/+5 km/s (stat-only, all-3 SPS pool) | ✓ |
| J7 | Sersic2D bound-fix (n ∈ [1.0, 6.0]) | `scripts/run_isource_shape_sweep.py:181-207`, 2026-05-11 | F200LP previously escaped to n=0.30 (unphysical flat-disk), inflating I-shape budget; multi-init grid + tight bounds fix it. I-shape ±5.4 → ±1.5 | ✓ |
| J8 | I-shape sweep at wR3800_5400_arcmask (10 shapes × 3 SPS × N=250) | `results/ishape_sweep_wR3800_5400_arcmask/` | std = 1.5 km/s — 9× tighter than narrow's ±13 | ✓ |
| J9 | Mask-weight sweep at wR3800_5400_arcmask (w ∈ {0, 0.5, 1.0} × 3 SPS × N=250) | `results/maskweight_sweep_wR3800_5400_arcmask/` | peak-to-peak/2 = 3.8 km/s — 4× tighter than narrow's ±16; spectral arc mask absorbs the worst | ✓ |
| J10 | Three-window cross-check at N=500 | nb09 §9: w6500_7500 vs wR3800_5400_arcmask vs wR4000_5400_arcmask | Spread = 15.0 km/s → window systematic ±15 (dominant residual). wR4000_5400_arcmask = orthogonal Hβ + Mg b feature set (no Ca H+K) | ✓ |
| J11 | Both-windows side-by-side budget | `notebooks/09d_final_systematics_both_windows.ipynb`, `results/sigma_e_final_systematics_nb09d.npz` | WIDE 254.85 ± 17.87 vs NARROW 267.95 ± 30.10; consistent at 0.4σ | ✓ |
| **K. Resolved kinematics (now feasible with wide arc-masked window)** | | | | |
| K1 | Per-spaxel ppxf inside R < 1.5 R_e at S/N≥{2,3,5,10} | `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` | S/N≥5 gives 17 spaxels, σ ∈ [144, 251] km/s with 0 outliers — first clean per-spaxel σ map for this target | ✓ |
| K2 | PowerBin spatial binning at wide arc-masked window | `notebooks/11 §PowerBin` | 7 bins at target S/N=15; 2 outer bins still hit σ > 400 (irreducible KCWI S/N limit at aperture edge) — improvement vs Test 19 narrow-window failure but not fully cured | ⚠ partial |
| K3 | Wide vs narrow window for resolved maps | `notebooks/11` | Wide arc-masked unlocks per-spaxel + PowerBin maps; narrow window was per-spaxel-infeasible | ✓ |
| **M. Reduction-pass systematic (2026-05-20)** | | | | |
| M1 | New `_mtwdo_` red-side reduction PROMOTED to headline 2026-05-26 | `notebooks/09e_new_red_reduction.ipynb`, `scripts/run_wide_sigma_e.py --cube new`, `raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits` | **σ_e(<R_e) = 265.76 −18.80 / +18.33 km/s** (with new reduction-pass component). Per-SPS spread 2.67 km/s. Δ vs old reduction = +10.92 km/s → carried as **±5.5 km/s reduction-pass systematic** (M2). Cube co-aligned to <0.002″ with the old reduction (centering ruled out), shift consistent across all 15 polynomial degrees + all 3 SPS templates (mechanism is distributed across the wide window, not localized to one λ boundary — confirmed by 8000-Å boundary-mask cross-test M3). Caches: `results/run_wide_sigma_e/new/`. Figures: `results/figures/nb09e_reduction_comparison.png`, `nb09e_grating_consistency_check.png`. Audit: `NOTES_nb09e_audit_2026-05-20.md` (Addendum 2026-05-26). | ✓ HEADLINE |
| M2 | Reduction-pass systematic budget component | M1 vs the old-reduction headline | half-Δ = (265.76 − 254.85) / 2 = **±5.46 km/s**, carried as **±5.5 km/s** in the systematic budget (matches our peak-to-peak/2 convention from E5 / J10). **Flag**: only 2 reductions available; revisit / refine if a 3rd lands. Documented in `METHODS_AND_SYSTEMATICS.md` Part I.9 + project memory `project_nb09e_reduction_systematic.md`. | ✓ |
| M3 | 8000-Å boundary-mask cross-test to localize the +10.9 km/s shift | `/tmp/test_8000A_boundary_mask.py`, caches at `results/run_wide_sigma_e/new_8000A_masked/` | σ_e = 268.73 km/s after adding extra mask at def-rest 4715–4834 Å (= obs 7900–8100 Å) → **hypothesis REJECTED** (shift NOT localized to this one transition). Recovery toward headline = −27%. Confirms the +10.9 km/s shift is distributed across the wide window — likely from the 2-tier per-night flux scaling tilting the entire continuum, not from a single λ boundary. | ✓ negative |
| M4 | ppxf clean=True sigma-clip cross-test on new cube (N=100) | `/tmp/test_sigma_clip_new_cube.py`, caches at `results/run_wide_sigma_e/new_sigma_clip_N100/` | σ_e = 265.80 km/s (Δ vs headline = +0.04 km/s). **0 pixels rejected** at every (sps, degree) of 45 fits (of ~2161 good pixels each). Reason: ppxf clean=True scales against the noise array; the noise array on this dataset is overestimated, so a 3σ threshold in noise units sits well above the actual outliers (see M5). Test by itself is uninformative; superseded by M5. | ✓ negative |
| M5 | Local-MAD bad-pixel cleaning (rolling 75-pix window, 3σ threshold) on new cube (N=100) | `/tmp/test_local_mad_clip.py`, caches at `results/run_wide_sigma_e/new_local_mad_clip_N100/` | **52 outlier pixels flagged** (of 2161 good); biggest is 6-pix cluster at def-rest 5232–5236 Å (obs 8768–8774 Å) at **26σ in local-MAD units** — clear unresolved cosmic-ray residual that the noise-scaled clean=True missed. σ_e shifts to **267.18 km/s** (Δ vs headline = **+1.42 km/s**), FSPS-dominated (+3.67). Carried as a stated sensitivity in METHODS Part I.9 ("What is not in the budget"); not folded into formal budget because the shift is well below stat 1σ (±6.1). | ✓ sensitivity |
| M6 | OLD-cube replication of M5 clip (same 52 def-rest λ) | `/tmp/test_local_mad_clip_OLD_cube.py`, caches at `results/run_wide_sigma_e/old_local_mad_clip_N100/` | σ_e shifts 254.85 → **256.01 km/s** (Δ = +1.16 km/s on OLD cube). Both cubes shift by ~+1.3 km/s → **52 bad pixels are intrinsic to the data**, not reduction-specific. Reduction-pass gap preserved (10.91 un-clipped → 11.17 clipped). | ✓ replicates |
| M7 | Noise-scaling + ppxf clean=True scan on new cube | `/tmp/test_noise_scaled_clean.py`, caches at `results/run_wide_sigma_e/new_noise_scaled_clean_N100/` | Scanned noise×{1.0, 0.5, 0.3, 0.2, 0.1, 0.05} at FSPS deg=22 → ppxf clips {0, 1, 18, 95, 591, 1247} pixels respectively. Full pipeline at noise×0.3 (19 pix/fit, closest to M5's 52) gives σ_e = **271.04 km/s** (Δ = +5.28 km/s) — larger shift than local-MAD because ppxf clean=True is iterative and identifies a different (smaller) set of pixels. Sensitivity envelope: cleaning shifts range +1 to +5 km/s depending on method. | ✓ sensitivity |
| M8 | He I 3819 source-emission mask added to ARC_MASKS_REST (5237–5253 Å) (2026-05-27) | `scripts/run_wide_sigma_e.py --cube {new,headline}_clean_hei`, caches `results/run_wide_sigma_e/{new,headline}_clean_hei/` | Identified 3-pixel +ve residual cluster at def-rest 5244–5248 Å (= source-rest 3818–3820 Å, matches **He I 3819.60** emission from z=1.302 source). Consistent across NEW + OLD reductions → astrophysical, not CR. Adding mask shifts NEW: 268.98 → **271.87** (+2.89), OLD: 260.44 → **264.16** (+3.72). Reduction-pass gap shrinks 8.54 → 7.71 → reduction-pass systematic refined ±4.27 → **±3.86**. New TOTAL sys = ±17.25; **new paper headline σ_e(<R_e) = 271.87 ± 17.86 km/s** (asym −17.99/+17.74). | ✓ HEADLINE update |
| M9 | def-rest 5193–5210 Å bump investigation — does NOT mask (it's the Mg b LOSVD wing) | `/tmp/smoke_extra_mask_5190_5210.py` (N=100 smoke); cube-arc spectrum extraction in `AGEL0206_arc_source_spectrum.pdf`; sky-line noise diagnostic | Visible +5–7% bump in deflector data above ppxf fit at def-rest 5193–5204 Å. **NOT source emission** (arc-spaxel integrated spectrum in source-rest 3782–3788 Å is flat, no emission line — see fig). **NOT sky residual** (cube `noise_sky` shows 0.6–0.9× median std in this range — well below the strong OH airglow forest at obs 8760–8790 Å that's already covered by BAD_PIXELS + the He I mask). **Diagnosis: Mg b LOSVD wing** — the bump is mirrored on the blue side of Mg b (def-rest 5160–5170 also has +ve residuals). Mg b at 5183 Å convolved with σ ≈ 270 km/s LOSVD broadens 4–5 Å → wings extend out to ~5190–5210 Å. Smoke test: masking 5190–5210 drops σ_e to **264.83 km/s (−7.04 km/s)** — the LARGEST single mask effect seen. **Verdict: DO NOT MASK** — this region is the *signal* ppxf uses to constrain σ_e, not contamination. Documenting so we don't re-investigate. | ✓ DO NOT MASK |
| M10 | Full sky-line audit across fit window 3800–5400 Å (2026-05-27) | `BAD_PIXELS_REST` updated to 35 entries (26 original + 9 added). Caches `results/run_wide_sigma_e/{new,headline}_clean_hei/`; per-band rationale inline in `scripts/run_window_sweep.py`. Identification via cube `noise_sky` thresholded at >2.5× median; each added band cross-checked against arc spectrum (`AGEL0206_arc_source_spectrum.pdf`) to confirm no source-emission counterpart. | Found 9 unmasked OH airglow / sky-residual bands in the fit window after user-flagged def-rest ~4600 Å structure: (4380.0-4384.0), (4553.0-4554.6), **(4602.0-4610.0)** ← the 4600 Å feature, (4687.3-4688.3), (4691.5-4693.6) merged, (4770.8-4771.8), (4951.6-4955.0), (5011.9-5015.3), (5029.8-5033.8). All in def-rest Å, sky_noise 2.5–4× median. NEW cube σ_e: 271.87 → **269.62** (−2.25 km/s). OLD cube σ_e: 264.16 → **262.72** (−1.44 km/s). Reduction-pass gap shrinks 7.71 → 6.90 → refined ±3.86 → **±3.45**. New TOTAL sys = ±17.16; new paper headline **σ_e(<R_e) = 269.62 ± 17.78 km/s** (asym −17.92/+17.65). | ✓ HEADLINE update |
| M11 | Rigorous re-derivation of carried I-shape, F200-mask, fit-window systematics on NEW cube + M10 masks (2026-05-27) | Caches `results/ishape_sweep_wR3800_5400_arcmask_M10/` (10 I-shapes × 3 SPS × N=250), `results/maskweight_sweep_wR3800_5400_arcmask_M10/` (3 weights × 3 SPS × N=500), `results/nb09a_wavelength_sweep_M10/` (3 windows × 3 SPS × N=500). Total compute ~1.5 h via `/tmp/run_full_systematics_barrage.sh`. | **I-shape:** range 266.83–271.37 km/s → peak-to-peak/2 = **±2.27** (was ±1.5 carried). **F200 mask:** w00=269.69, w50=261.95, w100=256.39 → peak-to-peak/2 = **±6.65** (was ±3.8 carried; bigger because cleaned NEW cube + Balmer-unmasked is more arc-sensitive). **Fit-window:** wR3800_5400_arcmask=269.66, wR4000_5400_arcmask=268.58, w6500_7500=276.22 → peak-to-peak/2 = **±3.82** (was ±15 carried; cleaned NEW cube + M10 brings wide/narrow into agreement). **Net:** sys total ±17.16 → **±10.81**; TOTAL ±17.78 → **±11.77** (drop of −6.0 km/s). **NEW HEADLINE: σ_e(<R_e) = 269.62 ± 11.77 km/s (asym −11.98/+11.57).** | ✓ HEADLINE update |
| M12 | R_e-source fold-in + masking-approach + Hδ cross-checks (2026-06-08, N=500) | `scripts/run_sigma_e_Re_systematic_wide.py` (nb15), `run_sigma_e_mask_systematic.py` (nb13), `run_sigma_e_hdelta_test.py` (nb16); registry `scripts/paper_values.py`→`results/PAPER_VALUES.json` | **(1) D7 R_e-source folded in: ±6.13** (wide-window, 4 estimators; user chose full spread 2026-06-08) → sys ±10.81 → **±12.43**; **HEADLINE σ_e = 269.62 ± 13.27 km/s (asym −13.45/+13.10).** **(2) Arc-masking-approach systematic = ±5.85** (expert/sersic/perband/global reprojected to IFU grid; spectroscopic analogue of the M★ ±0.16 dex term) — overlaps F200-mask ±6.65, NOT added (larger-of-two). **(3) Hδ: KEEP UNMASKED** — local-MAD shows Hδ well-fit (max\|resid/noise\| 0.44 < global 0.81), masking shifts σ_e +6–8 km/s = LOSVD information not contamination (M9 pattern); TODO in `bootstrap_ppxf.py` closed. | ✓ HEADLINE update |
| **M13** | **"Best mask throughout" — best-mask R_e=2.097" cascade (2026-06-11, N=500)** | `scripts/run_sigma_e_Re_grid.py`, `validate_Re_photutils.py`, `Re_bestmask_reconciliation.py`, `aperture_{matched_photometry,2re_companions}.py`; registry `paper_values.py` | **Adopt the best-mask (single-Sérsic + color/morph + WCS-reg) CoG R_e=2.097" (was 2.305") as the headline everywhere.** **(1) σ_e re-measured at R_e=2.097" → 267.31 −12.98/+12.61 (sym ±12.79)** (was 269.62 ± 13.27; −2.3 along rising σ(R)). **(2) R_e-source re-derived** on a 7-pt grid → **adopted best-mask CoG light family ±5.13** (CaHK+G/full-grid ±9.98 = cross-checks, not folded; user 2026-06-11). **(3) R_e/CoG photutils-validated ±0.002"** (A3d). **(4) M★ re-measured at 2 R_e=4.19" → total 11.47 (unchanged, R_e-robust).** R_e method sys ±0.100". | ✓ HEADLINE update |
| **L. Figure preparation (2026-05-13)** | | | | |
| L1 | Figure 2 (narrow) — single posterior inset | `AGEL_0206_ApJL_Figures/figures.ipynb` cell `a546db7f`, `Mbh_sigma_SED_final.pdf` (was AGEL0206_sigma_e_SED_final.pdf) | Title shows σ_e = 267 ± 24(stat) ± 18(sys) km/s; no-arc-mask overlay removed | ✓ |
| L2 | Figure 2 (wide, headline) — single posterior inset | figures.ipynb cell `fig2_wide`, `AGEL0206_sigma_e_SED_final_wide.pdf` | σ_e = 255 ± 6(stat) ± 17(sys); residuals panel; 9 absorption features labeled; 4 arc-emission masks shaded | ✓ |
| L3 | Figure 3 (M_BH-M_star) cleanup | figures.ipynb cell `fig3_clean_no_err`, `Mbh_Mstar_relation_clean.pdf` | No per-object error bars; uniform color for local sample; filtered to Greene+2020 E+S0 list (60 unique, 4 KH13 suffix variants matched; 19 post-KH13 extras dropped) | ✓ |
| L4 | Figures 2 + 4 updated to principled-mask M★ (2026-05-29) | figures.ipynb cells 19/35/36; `AGEL0206_sigma_e_SED_final_wide.pdf`, `Mbh_Mstar_relation{,_clean}.pdf` | Fig2 inset posterior → perband_raw_10pct (median 11.16) + asymmetric title 11.16 +0.32/−0.08 + shaded fill-in reach to 11.46; Fig4 AGEL point → 10^11.16, asymmetric xerr −0.08/+0.32. Backup `figures.ipynb.bak.Mstar_principled_2026-05-29` | ✓ |

| **N — Photometry arc-mask verification + principled masking (2026-05-29)** | | | | |
| N1 | F200LP hand mask reproducible from objective selectors | `scripts/arc_mask_verification.py`, `notebooks/12` Part I | Color (m_F200LP−m_F140W) + Sersic-residual reproduce the expert F200LP mask; Sersic-residual (k=3) matches expert photometry to **0.016 mag**, R_e to 3% | ✓ |
| N2 | S/N-regime sweep (over-masking quantified) | `scripts/arc_mask_verification.py`, fig `arcmask_snr_sweep.png` | snr∈{2..20}, k∈{2..8}; F200LP model-flux-fraction stays 0.03–0.10 (clean), Δmag saturates at k≈3 = natural stopping point | ✓ |
| N3 | Independent per-band Sersic OVER-MASKS F140W/JWST | `scripts/principled_mask_photometry.py`, figs `principled_*.png` | Single Sersic under-fits bright IR galaxy (median resid +1.0) → over-masks 0.21–0.36 mag; low mdlfrac (0.02–0.03) is a FALSE guard. Use F200LP-reproj instead | ✓ negative |
| N4 | Reprojected F200LP footprint = transferable standard | `scripts/mask_method_comparison.py`, figs `maskcompare_*.png` | Reproducing expert JWST mags to **0.01–0.02 mag** by reprojecting F200LP (NOT union with F140W — F140W mask covers core, 17px at r<0.4″) | ✓ |
| N5 | 2-component (bulge+disk) Sersic model | `scripts/photometry_systematics.py:fit_2comp` | Single Sersic under-fits JWST; 2-comp cuts galaxy-body RMS **3.3–3.5×** (F150W2 3.25×, F322W2 3.46×). No PSF (env lacks webbpsf) — flagged | ✓ |
| N6 | IR source extension + per-band vs global mask | `scripts/photometry_systematics.py`, figs `photsys_*.png` | Deeper IR extends arc (F150W2 +9119px). **raw is mask-size-dependent** (up to 0.78 mag per-band↔global); **filled is mask-def-independent** (≤0.06 mag). Fill-in correction +0.18–0.96 mag | ✓ |
| N7 | M★ budget — raw-central, one-sided-up systematic | `results/photometry_systematics_Mstar.npz`, fig `Mstar_budget.png`, `notebooks/12` Part II | **log M★ = 11.16 ± 0.08(stat) +0.31(sys)** at 10% [**11.04 ± 0.14 +0.32** at 20%]; sys = Sersic fill-in (under-arc light), reach 11.46/11.36. per-band vs global negligible for filled. **NEW HEADLINE M★** (supersedes 11.33) | ✓ HEADLINE update |
| N8 | Explicit masking-approach systematic on M★ | `results/Mstar_masking_systematic.npz`, fig `Mstar_masking_systematic.png` | **±0.16 dex** (10%: ±0.160, 20%: ±0.162) = peak-to-peak/2 across {expert-aper, per-band/global × raw/filled} spanning 11.15–11.47. Decomp: under-arc(raw↔filled) ±0.15, mask-def(per-band↔global) ±0.004, mask-extent(expert↔IR-ext) 0.18. Symmetric statement of the headline's one-sided +0.31. **TODO: σ_e analogue** | ✓ |
| N9 | PSF-convolved fill model (firm up the fill leg) | `scripts/psf_fill_model.py`, `results/psf_fill_model.npz` | Gaussian PSF at instrument FWHM (F200LP 0.08″/F140W 0.14″/F150W2 0.06″/F322W2 0.13″; env lacks webbpsf). **ΔPSF(filled) ≤ 0.004 mag** all bands → ΔM★ ≪0.01 dex, negligible vs stat ±0.08 / masking ±0.16. Fill is PSF-robust (arc outside PSF core) | ✓ |
| A3b | R_e pixscale fix in measure_Re.hst_Re (diagnostic) | `scripts/measure_Re.py` (now reads WCS; F200LP literal 0.08→0.05) | F200LP cutout is 0.05″/pix (was hard-coded 0.08) → diagnostic F200LP R_e 3.05→**1.91″**; F140W unchanged. **Headline R_e=2.305″ UNAFFECTED** (final_sigma_e.py already WCS-correct: F200LP=2.52″, F140W=2.16″, mean≈2.34″). Flag: two CoG algos differ (measure_Re 1.91 vs final_sigma_e 2.52 for F200LP) — pre-existing | ✓ fix |
| A3c | CoG-algorithm reconciliation (measure_Re vs final_sigma_e) | `notebooks/14_CoG_reconciliation.ipynb`, `scripts/measure_Re.py` (`hst_Re` now `r_max_arcsec=6.0` default), `results/Re_cog_reconciliation.npz` | The two CoG algos are the **same up to `r_max`**: measure_Re capped at 4.0″ (`arange(1,80,1)` px), final_sigma_e integrates to 6.0″. **Matched r_max → agree <0.04″** (residual = centroid_2dg vs fixed-RA/Dec center). **Headline R_e=2.305″ = canonical** (final_sigma_e, r_max=6″). **New flag:** raw CoG R_e is r_max-non-convergent (non-zero sky pedestal); headline sits at the top of a 2.1–2.5″ method family (sky-sub 2.14″, Sérsic 2.18″) → **R_e method systematic ±0.08″**. Not folded (would shift σ_e aperture + M★; flagged for user). | ✓ reconciled |

---

## 1. Headline numbers

### σ_e at three apertures (paper headline = wide arc-masked window)

| Aperture | σ_e [km/s] (WIDE, headline) | σ_e [km/s] (NARROW, cross-check) | Notes |
|---|---|---|---|
| R<R_e (=2.305") | **255 ± 18** (254.85 ± 17.87) | 268 ± 30 (267.95 ± 30.10) | **paper headline** — KH13/Greene+20 σ_e at WIDE; NARROW cross-check at 0.4σ |
| R<R_e/2 (=1.152") | ~225 (nb07c §6cum NARROW) | 224 ± 18 | gradient diagnostic (still cited from nb07c) |
| R<R_e/8 (=0.288") | dropped | dropped | inside seeing FWHM/2 = 0.64" |

### Error-budget composition (R<R_e) — both windows side by side

| Component | WIDE wR3800_5400_arcmask | NARROW w6500_7500 | Source |
|---|---|---|---|
| Statistical (bootstrap pooled 1σ) | ±6.1 | ±23.9 | wild-bootstrap N=500, FSPS+EMILES+XSL pool |
| I-shape spread (10 shapes, post-Sersic2D-bound-fix) | ±1.5 | ±3.7 | `results/ishape_sweep_wR3800_5400_arcmask/` + `results/annular_bootstrap_07c_ishape/` |
| F200 mask on/off Δ (peak-to-peak / 2) | ±3.8 | ±7.5 | `results/maskweight_sweep_wR3800_5400_arcmask/` + `results/mask_weight_sweep.npz` |
| Frame fix (max σ shift across SPS, carried) | ±5 | ±5 | NOTES addendum + audit 1 |
| Centering (5-center sweep ±0.4", carried) | ±4 | ±4 | `NOTES_centering_investigation_2026-04-27.md` |
| Fit-window (3-window spread from §9) | ±15 | ±15 | nb09 §9: w6500_7500 / wR3800_5400_arcmask / wR4000_5400_arcmask (orthogonal Hβ+Mg b) |
| **Quadrature total** | **±17.87** | **±30.10** | `results/sigma_e_final_systematics_nb09d.npz` |

### Where each component changed between windows

- **Statistical ±24 → ±6.1** (4× tighter): wide window has ~4× more spectral pixels (2161 vs 555 good pixels). This is the headline reason to elevate the wide window as primary.
- **I-shape ±13 → ±1.5** (9× tighter): driven mostly by the Sersic2D bound-fix (J7), which removed the n=0.30 pathology that contributed ±5.4 of the old narrow ±13. The wider window also stabilizes σ across I-shape choices.
- **Mask ±16 → ±3.8** (4× tighter): the spectral arc-emission masks (J3) absorb most of what the F200 spatial mask used to absorb; the F200 mask becomes a sub-budget perturbation rather than a primary uncertainty.
- **SPS template spread**: unstated as a separate line in the budget because it's contained in the statistical pool. But it collapsed 26 → 4 km/s between windows (J5) — was a major budget item at the narrow window, sub-statistical at the wide.
- **Frame fix, centering**: unchanged (carried from the narrow-window era's careful audits — physically independent of the fit window).
- **Fit-window ±15**: NEW budget component (didn't exist before nb09d because there was only one window). Spread between three physically meaningful windows: Ca H+K + G (narrow), full canonical range (wide arc-masked), and Hβ + Mg b only (wide minus Ca H+K). This is now the dominant residual systematic.

### Sensitivity values (for paper text)

WIDE arc-masked window (headline):
- No-spatial-mask (w=1.0) at R<R_e: σ_e = 247.28 km/s — Δ from headline (w=0) = −7.6 km/s
- Half-spatial-mask (w=0.5) at R<R_e: 250.15 km/s — Δ = −4.7 km/s
- Per-SPS at wR3800_5400_arcmask: FSPS=253.9, EMILES=253.3, XSL=257.5 → spread 4.2 km/s

NARROW Ca H+K + G window (cross-check, retained from nb07c era):
- Soft-mask (w=0.5) at R<R_e: 258.5 km/s — Δ from narrow-headline = −9.3 km/s
- No-mask (w=1.0) at R<R_e: 252.8 km/s — Δ = −15.0 km/s
- σ_e(w) is essentially linear with weak super-linearity (quadratic c=+7.3)
- Per-SPS at w6500_7500: FSPS=253.7, EMILES=267.9, XSL=279.6 → spread 26.0 km/s (the SPS systematic that motivated the move to wide)

---

## 2. Detailed test catalog

This expands each row of the matrix above. Subsections reference scripts,
notebooks, and result files — open the linked file for the detailed code
and the inline plots/output.

### A. Foundational data and geometry

**A1 KCWI cube provenance** —
The cube `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` is the final
multi-night reduction stacked by K. R. Gupta on UT 2026-01-07. The
header's `DATE-OBS = 2025-11-17`, `PROGNAME = U002`, `PROGPI = Jones`
describes only the **first input frame** (`kr251117_00129`), not the
full input set.

**Verified contributing nights** (cross-checked 2026-05-26 against the
Keck observing logs):

| Night (UT) | Program | PI | DESJ0206 frames | On-source RED + BLUE |
|---|---|---|---|---|
| 2025 Aug 30 | K409 | TBD | 12 RED `kr250830_00090–00101`, 4 BLUE `kb250830_00052–00055` | 60 min + 66 min |
| 2025 Nov 17 | U002 | Jones | 20 RED `kr251117_00129–00148`, 5 BLUE `kb251117_00087–00091` | 90 min + 110 min |
| 2024 Dec 29 | TBD | (Yuguang Chen) | 4 RED `kr241229_00092–00095` per local stack list | **NOT confirmed as DESJ0206 pointings** — only the 44 MB image-only PDF `kcwi-241228-written-log.pdf` (Drive id `16zWO8yaCIvRtoCGenNTndvB7dYdrCkMH`) exists; no machine-readable log. Pending raw-frame header dump. |

Confirmed total on-source (Aug 30 + Nov 17): ≈ 170 min RED + 176 min BLUE.

Two earlier-recorded contributors are **incorrect**: ~~Sept 29 2025 K409~~
has zero DESJ0206 entries in its log; ~~Dec 29 2024~~ remains
unverified. Cite K409 + U002 in the paper, flag Dec 29 as pending. See
`reference_kcwi_data_properties.md` (memory) for full Drive IDs and
log-file pointers.

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
Bagpipes SED fitting on HST + JWST aperture photometry (notebook 02,
N_posterior = 500). Bagpipes config: exponential-τ SFH, BPASS+CLOUDY
templates, Calzetti dust, free metallicity, redshift fixed at 0.67564.

Per-filter observed flux + best-fit modelled flux:

| Filter | λ_pivot [Å] | AB obs | AB model (p50) | F_obs [10⁻¹⁸ cgs] | F_model (p50) | residual |
|---|---|---|---|---|---|---|
| F200LP | 4 972 | 22.613 | 22.563 | 3.97 ± 0.06 | 4.16 | +4.7% |
| F140W | 13 923 | 19.134 | 19.264 | 12.47 ± 0.05 | 11.06 | −11.3% |
| F150W2 | 16 720 | 18.942 | 19.033 | 10.31 ± 0.005 | 9.49 | −8.0% |
| F322W2 | 32 470 | 18.604 | 18.505 | 3.73 ± 0.001 | 4.09 | +9.6% |

Posteriors:

| Parameter | p16 | p50 | p84 |
|---|---|---|---|
| **log M★/M☉** | **11.24** | **11.33** | **11.41** |
| mass-weighted age [Gyr] | 3.69 | 5.10 | 6.11 |
| exponential SFH age [Gyr] | 4.46 | 5.98 | 6.97 |
| exponential τ [Gyr] | 0.43 | 0.68 | 1.21 |
| metallicity Z/Z☉ | 0.17 | 0.89 | 1.87 |
| A_V [mag] | 0.34 | 0.82 | 1.55 |
| z_phot | 0.674 | 0.675 | 0.676 |

→ M★ = 2.15 × 10¹¹ M☉; z_phot agrees with the line-fit z = 0.67564
within ±0.001.

The model under-predicts the two NIR bands (F140W, F150W2) by ~8-11%
and over-predicts the bluest+reddest by ~5-10%. Classic "redder NIR
slope than a single-τ SFH prefers" tension; not catastrophic, but
flags that a more flexible SFH (delayed-τ or non-parametric) might
tighten the fit if revisited.

**Sersic-total cross-check** (`notebooks/08_sersic_total_photometry.ipynb`,
`scripts/measure_Re.py`-style 2D Sersic fits per filter, extrapolated
to total flux):

| Filter | AB aperture | AB Sersic-total | Δmag | Sersic n | r_eff [″] |
|---|---|---|---|---|---|
| F200LP | 22.613 | 20.672 | −1.94 | 1.42 | 2.49 |
| F140W | 19.134 | 18.282 | −0.85 | 1.54 | 1.88 |
| F150W2 | 18.942 | 18.154 | −0.79 | 1.40 | 1.72 |
| F322W2 | 18.604 | 17.633 | −0.97 | 1.97 | 2.03 |

Re-fed through Bagpipes with the same priors:

| Quantity | aperture (paper) | Sersic-total | Δ |
|---|---|---|---|
| log M★/M☉ p50 | 11.33 +0.07/−0.09 | 11.40 +0.11/−0.15 | +0.065 dex |
| M★ [M☉] | 2.15 × 10¹¹ | 2.49 × 10¹¹ | +16% |

The aperture under-counts outer-profile flux by ~0.8-1.9 mag depending
on filter, but at the integrated-mass level the Sersic correction is
+0.07 dex — within the aperture-photometry quoted uncertainties. The
paper cites the aperture value as the headline; Sersic-total is the
sensitivity floor.

Cache: `results/bagpipes_sed_results.npz` (aperture, 500-sample
posterior), `results/bagpipes_sersic_refit.npz` (Sersic-total),
`results/sersic_total_photometry.npz` (Sersic fit parameters).
Figure: `results/AGEL0206_spectra_SED_fit.pdf`.

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
single 1-D spectrum, fit ppxf. This is the headline path — the co-add-then-fit
method of Cappellari et al. 2006 (SAURON IV §2.3), as KH13 / Greene+20 / SAURON /
ATLAS3D / MaNGA compute σ_e (the σ_e definition integral is Gültekin 2009 eq. 1;
**NOT** Cappellari 2006 eq. 1, which is the aperture correction). Cross-references in
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

**D7 R_e source systematic (mean vs F140W vs F200LP vs Ca H+K + G)** —
The headline integrates the I-weighted aperture spectrum inside the
**mean** F140W+F200LP masked CoG R_e (=2.305"). Three alternative R_e
definitions test sensitivity to that choice (Track A masked, N=500):

| R_e source | R_e ["] | σ_e [km/s] | Δσ_e |
|---|---|---|---|
| mean (paper)         | 2.305 | 267.82 | — |
| F140W only           | 2.168 | 264.87 | −2.96 |
| F200LP only          | 2.441 | 275.86 | +8.04 |
| Ca H+K + G-band depth| 2.866 | 281.74 | +13.92 |

Spread = 16.9 km/s — at the mask-budget level (±16) but sub-budget
vs the total ±32. σ_e rises monotonically with R_e (more outer-bulge
spaxels in the I-weighted aperture). The Ca H+K + G-band absorption-
depth I-map (rest 3925-3942, 3960-3976, 4297-4313 Å, summed
`(continuum − flux)` per spaxel — see `cahk_g_line_depth_map` in
`scripts/final_sigma_e.py`) is the most stellar-only of the four;
its larger R_e reflects that the deflector light extends further than
the F200LP UV-leaning band suggests.
Code: `scripts/final_sigma_e.py` §7 + `cahk_g_line_depth_map`.
Display: `notebooks/09 §7b` + `results/figures/nb09_re_source_systematic.png`.
Caches: `results/final_sigma_e_paper/Re_{F140W,F200LP,CaHK}_{sps}_N500.npz`.

### E. Mask treatment

**Masking vs no-masking — pro/con summary**

| Track | Pro | Con |
|---|---|---|
| Hard mask (w=0.0, headline) | Removes arc contamination; matches literature convention; F200LP mask is purpose-built for this lens | Drops 38 of 184 spaxels inside R<R_e (~21% by count, ~27% by I-weight) → lower S/N; hard boundary may over-clip a few edge spaxels |
| No mask (w=1.0, sensitivity) | Maximum S/N; no boundary artifacts | Arc adds low-velocity continuum dilution → σ biased low by 15-17 km/s (Δ_mask = −16 km/s at R<R_e) |
| Soft mask (w=0.5) | Smooth middle ground | Super-linear in w (quadratic c=+7.3); threshold-dominated bias → first 25% of arc weight already captures 34% of total mask sensitivity. Arbitrary weight choice; not a literature convention |
| Mask dilation (E6) | — | Inflates σ via noise (dropped spaxels are S/N contributors, not bias contributors). **Negative finding** — do not dilate |

The 16 km/s no-mask Δ is propagated into the error budget as
`σ_mask = ±16`, so the choice is bounded, not assumed.

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
*residual* arc dilution (the small fraction that leaks through the
F200 mask) is sub-dominant at production statistics. The F200 mask
already does the heavy lifting. Note: this test compares
"masked + arc-subtract" against "masked alone" — it does NOT directly
validate the no-mask vs masked Δ mechanism (that's E8). Notebook:
`nb07e`.

**E8 Arc-as-sky mechanism test (rigorous test of the no-mask Δ)** —
Run §6cum at the *no-mask* R<R_e aperture with the bright outer-arc
spectrum passed to ppxf as a `sky` template (free-amplitude additive
component, NOT convolved with the deflector LOSVD — physically
appropriate since the arc is at z=1.302, decoupled from the deflector
kinematics). At N=500, FSPS+EMILES+XSL pooled:

| Track | σ_e [km/s] | Δ from σ_A |
|---|---|---|
| A: masked headline                 | 267.82 | — |
| B: no-mask                         | 252.82 | −15.00 |
| **D: no-mask + arc-sky (E8)**      | **274.64** | **+6.82** |

Recovery fraction (σ_D − σ_B)/(σ_A − σ_B) = **145%** — i.e. arc-sky
**overshoots** the masked headline by ~7 km/s. Interpretation:

- Continuum dilution explains the FULL no-mask Δ. There is no
  detectable kinematic-blend or template-mix mechanism contributing in
  the *opposite* direction.
- The 7 km/s overshoot is consistent with the known caveat that the
  outer-arc-mask spaxels (R > 3R_e/4) still receive seeing-blurred
  deflector light at the few-percent level, so subtracting α × arc
  oversubtracts a small amount of *real* deflector flux too. This
  inflates σ slightly.
- The masked headline (Track A, paper number) remains the cleaner
  measurement because it doesn't have this oversubtraction artefact.
- α_arc ≈ 0.31 across all three SPS — internally consistent.

Code: `scripts/run_07f_arc_sky.py` (uses `setup_ppxf_inputs_from_spectrum`
+ ppxf `sky=` argument). Display: `notebooks/07f_arc_sky_subtract.ipynb`.
Caches: `results/arc_sky_07f/{sps}_N500.npz` + `results/arc_sky_07f.npz`.
Figures: `results/figures/nb07f_{arc_template,recovery,posteriors}.png`.

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

**H1 §7 discrete Gültekin (refreshed 2026-05-01)** — `notebooks/07c §7`,
recomputed with frame-aware ppxf via `scripts/refresh_07c_gultekin.py`.
Annular sum σ_e²(<R) = Σ F_j (V_j² + σ_j²) / Σ F_j with arc-filter on
outer annulus.

| | Pre-frame-fix (Apr 24) | Post-frame-fix (May 1) | Δ |
|---|---|---|---|
| σ_e (arc-filtered) | 254.99 −24.2/+28.4 | **256.17 −13.0/+12.7** | +1.18 |

Δ < ±5 km/s frame budget. Errors **halved** because the V_sys split
collapsed from ~99 km/s (pre-fix) to ~8 km/s (post-fix per-SPS V_sys:
fsps −11.2, emiles −7.4, xsl −3.1), shrinking the V² contribution to
the Gültekin sum's MC scatter.

**H2 §7b flat-σ extrapolation (refreshed 2026-05-01)** —
`notebooks/07c §7b`, same refresh. Push the §7 sum into the outer
flagged bin assuming σ_outer = σ at the last clean annulus.

| | Pre-frame-fix | Post-frame-fix | Δ |
|---|---|---|---|
| σ_e | 271 −33/+35 | **274.37 −16.2/+17.4** | +3.37 |

Both H1 and H2 remain consistent with the §6cum headline (267.82) at
<1σ. The point of these cross-checks is preserved by the refresh, and
the post-fix bands are tighter, strengthening rather than weakening
the consistency.

Note: nb07e (E7, arc-spectrum-subtraction sibling) was *not* refreshed
because its result is a **relative comparison** (matches §6cum within
0.1 km/s) — a uniform per-SPS σ shift moves both arms of the
comparison together and the relative match is unchanged.

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

### J. Wide arc-masked window (current paper headline, May 2026)

This section documents the May 8–13 work that elevated the rest 3800–5400 Å
arc-masked window above the narrow Ca H+K + G-band window as the paper
headline. The wide window had been off-limits in the nb07c era because of
source-emission contamination from the z=1.302 spiral and template
issues — both are now resolved.

**J1 Wavelength-window sweep (nb09a)** —
`notebooks/09a_wavelength_window_sweep.ipynb` + `scripts/run_window_sweep.py`.
Bootstrap σ_e at 15 candidate windows from def-rest 3500–5500 down to
narrow Ca H+K only. Result: wR3800_5400 wins on statistical precision
(~4× more spectral pixels than narrow). Blue edge set by **XSL template
floor at def-rest 3675 Å with 5% velocity padding** (J2 below); red edge
set by KCWI red-channel coverage at z=0.67564.

**Provenance of the 3800–5400 range:** NOT from a specific literature
prescription — it's set by the data (XSL template floor + retaining Ca H+K
to cross-anchor with the narrow w6500_7500 result). Closest literature
precedent for going as blue as 3800 Å is **LEGA-C (Bezanson+2018,
van der Wel+2021)** at z~0.7 using ~3500–5500 rest. SAURON/ATLAS3D's
canonical kinematic window is much narrower (4800–5380; Mg b + Fe5270 +
Fe5335 only). The Lick-tradition "4000–5400 (no Ca H+K)" range is
captured by the `wR4000_5400_arcmask` sensitivity check (orthogonal
Hβ + Mg b feature set) — NOT the headline.

**J2 Source-emission contamination discovery (nb09a §3a, nb09b §8)** —
Wide-window ppxf wild-bootstrap posterior was bimodal (σ ≈ 250 alongside
σ ≈ 100–150; frac<200 km/s = 60–80%) at deg ≥ 23. Caused by source-frame
emission lines from the z=1.302 spiral arc contaminating the deflector
spectrum. Empirically identified as >4σ excess above local continuum:

- **def-rest 3835–3855 Å** = Mg II λλ2796/2803 doublet from z=1.302
- **def-rest 4525–4545 Å** = **O₂ A-band telluric leading-edge residual at
  obs 7593–7626 Å** (NOT source emission — originally flagged as a
  ~7.9σ wild-bootstrap excess and provisionally tagged "source rest
  3300 Å feature"; relabeled 2026-05-18 after the spike was confirmed
  present at the same wavelength in deflector R<1.5″, the CLAUDE.md sky
  box, and a 4–8″ off-deflector annulus, with deflector/off-source
  amplitude ratio 1.11× despite a 17.5× continuum-brightness ratio. See
  `NOTES_4534A_spike_investigation_2026-05-18.md` and
  `results/figures/spike_4534A_diagnostic.png`.)
- **def-rest 5115–5135 Å** = source [O II] λλ3727/3729 — sits **on the
  Mg b absorption**; was distorting the FSPS Mg b depth ratio to 1.18
- **def-rest 5260–5340 Å** = Mg b / [Ne III] λ3869 cluster from source

Total: 212 of 2374 pixels (8.9% of fit window) masked.

**J3 Effect of the arc mask (nb09b §8)** —
With arc-mask + ppxf clean=True:
- σ vs polynomial degree FLAT across deg 15–29 (was bimodal at 250↔175
  with a break at deg 23)
- N=200 bootstrap collapses to a single sharp peak (frac<200 km/s = 0%
  across all 3 SPS)
- clean=True drops 0 pixels at deg=21 — already clean after arc-mask
- σ_e (FSPS+EMILES+XSL pool) = **254.9 −7/+5 km/s** at wR3800_5400_arcmask

**J4 Per-SPS template spread collapse (nb09d §1.3)** —

| window | FSPS | EMILES | XSL | all-3 pool | SPS spread |
|---|---|---|---|---|---|
| NARROW w6500_7500 | 253.7 | 267.9 | 279.6 | 268.0 | **26.0 km/s** |
| WIDE wR3800_5400_arcmask | 253.9 | 253.3 | 257.5 | 254.8 | **4.2 km/s** |

At the narrow window the 3 SPS libraries fundamentally disagreed by
26 km/s with strong ordering (FSPS < EMILES < XSL — see
`reference_sps_systematic.md`). At the wide arc-masked window they agree
to 4 km/s, well below per-SPS bootstrap stat error. The SPS-template
uncertainty essentially disappears once ppxf has access to Mg b + Fe5270
and the broader feature set. Strong paper argument for elevating the
wide window. Figure: `results/figures/nb09d_per_sps_both_windows.png`.

**J5 Sersic2D bound-fix (`scripts/run_isource_shape_sweep.py`, 2026-05-11)** —
The I-shape sweep was inflated by a non-physical Sersic2D fit. Original
script allowed `n ∈ [0.3, 8.0]` and `ellip ∈ [0.0, 0.95]`. The F200LP fit
escaped to **n=0.30, ellip=0.00** — a degenerate flat-disk minimum
(astropy flagged it as unsuccessful). The flat profile under-weighted the
center and yielded σ = 271.8 km/s, a +20 km/s outlier inflating the
I-shape std to ±5.4. Fix:

- Tighten bounds: `n ∈ [1.0, 6.0]`, `ellip ∈ [0.0, 0.6]`,
  `r_eff ∈ [r_eff_pix × 0.5, r_eff_pix × 2.0]` (was ×0.3–×3.0)
- Try 3 initial-condition grids: (n_init, ellip_init) ∈
  {(1.5, 0.05), (2.0, 0.2), (3.5, 0.2)}; keep lowest-χ² fit
- F200LP now converges to a physical n that gives σ = 251.4 km/s
- F140W was already converging at n=2.14, fix did not change it materially

After the fix, I-shape std collapsed: **±5.4 → ±1.5 km/s** at wide
arc-masked. Universal fix applied to the I-shape sweep code path — all
downstream sweeps now use these bounds.

**J6 I-shape sweep at wR3800_5400_arcmask (10 shapes × 3 SPS × N=250)** —
`results/ishape_sweep_wR3800_5400_arcmask/{shape}_{sps}_N250.npz`. Ten
I(r) weight maps: IFU_band, IFU_wl, F140W/F200LP × {raw, arcmasked,
1Dcog, Sersic2D}. Std across the 10 shapes, pooled across 3 SPS = 1.5 km/s.
N=100 → N=250 changed result by 0.05 km/s — fully converged.

**J7 F200LP spatial-mask sensitivity at wR3800_5400_arcmask
(3 weights × 3 SPS × N=250)** —
`results/maskweight_sweep_wR3800_5400_arcmask/{w00,w50,w100}_{sps}_N250.npz`.
w=0 (hard mask) = 254.84, w=0.5 = 250.15, w=1.0 (no spatial mask) = 247.28.
Peak-to-peak / 2 = 3.8 km/s. **4× tighter than narrow's ±16 km/s budget**
because the spectral arc-emission masks now absorb the worst of the arc
contamination — the spatial mask is no longer the primary defense.

**J8 Three-window cross-check at N=500 (nb09 §9)** —
Final fit-window systematic. Three windows at full N=500:

| window | rest pixels | σ_e (all-3 pool) | role |
|---|---|---|---|
| w6500_7500 | obs 6500–7500 Å (rest 3879–4476) | 269.7 km/s | narrow Ca H+K + G only |
| wR3800_5400_arcmask | rest 3800–5336 | **254.8 km/s** | HEADLINE |
| wR4000_5400_arcmask | rest 4000–5400 | ~265 km/s | orthogonal Hβ + Mg b (no Ca H+K) |

Spread = 15 km/s → fit-window systematic = ±15 km/s. The
wR4000_5400_arcmask window is the Lick-tradition feature-set
cross-check: removes Ca H+K entirely, so any agreement between the
3800-anchored headline and 4000-anchored alternative is a check that
Ca H+K is NOT artificially driving the σ_e measurement. The 10 km/s
gap between WIDE and wR4000 is the dominant residual systematic.

**J9 Both-windows side-by-side budget
(`notebooks/09d_final_systematics_both_windows.ipynb`)** —
Cache-reuse: narrow I-shape at N=500 reused from `annular_bootstrap_07c_ishape/`
(only F200LP_Sersic2D and F140W_Sersic2D refit 2026-05-11 with new
Sersic bounds); narrow mask-weight at N=500 reused from `mask_weight_sweep.npz`;
both N=500 stat pools from `nb09a_wavelength_sweep/`. Total final budget
table (also `results/sigma_e_final_systematics_nb09d.npz`):

| component | NARROW w6500_7500 | WIDE wR3800_5400_arcmask |
|---|---|---|
| stat (N=500) | ±23.9 | ±6.1 |
| I-shape (10 shapes, N=250, Sersic2D bound-fix applied) | ±3.7 | ±1.5 |
| F200 mask (peak-to-peak / 2) | ±7.5 | ±3.8 |
| frame (vac/air, carried) | ±5.0 | ±5.0 |
| centering (HST WCS, carried) | ±4.0 | ±4.0 |
| fit-window (3-window from §9) | ±15.0 | ±15.0 |
| **TOTAL** | **±30.10** | **±17.87** |

**Headline: σ_e(<R_e) = 254.85 ± 17.87 km/s (WIDE) vs 267.95 ± 30.10
km/s (NARROW).** Difference 13.1 km/s, well within both budgets at 0.4σ.

### K. Resolved kinematics (per-spaxel + PowerBin at wide arc-masked window)

The wide arc-masked window unlocks resolved-kinematics measurements that
were previously infeasible at the narrow window. Documented in
`notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` and memory
`project_nb11_resolved_kinematics.md`.

**K1 Per-spaxel ppxf inside R < 1.5 R_e** —
Per-spaxel ppxf fits at S/N floors of {2, 3, 5, 10}. Results inside
R < 1.5 R_e = 3.46" aperture (358 spaxels total):

| S/N floor | N spaxels | % σ in [150,350] | median σ | usable? |
|---|---|---|---|---|
| ≥ 10 | 0 | -- | -- | max per-spaxel S/N is 8.6 |
| **≥ 5** | **17** | **94%** | **201 km/s** | **YES — production map** |
| ≥ 3 | ~76 | 68% | ~240 | usable as soft map with outlier flagging |
| ≥ 2 | ~134 | 53% | ~243 | mostly noise-dominated beyond 1 R_e |

S/N ≥ 5 is the cleanest result — 17 central spaxels with σ ∈ [144, 251]
km/s, 0 outliers, spanning the central ~1 R_e. σ(R) drops with radius —
consistent with elliptical bulge morphology. Aperture-edge spatial
resolution is fundamentally constrained by KCWI per-spaxel S/N, not
ppxf or fit window.

**K2 PowerBin spatial binning at wide arc-masked window** —
PowerBin (Cappellari 2025) at target capacity S/N=15. 7 bins inside
R < 1.5 R_e; median σ = 294.5 km/s. **2 of 7 outer bins still hit
σ > 400 km/s** — irreducible KCWI per-spaxel S/N limit at the aperture
edge. Major improvement over Test 19 narrow-window failure (which had
σ > 800 outer bins) but not fully cured. Flag outer bins as
noise-limited rather than attempting physical interpretation.

**K3 Wide vs narrow for resolved maps** —
At the narrow w6500_7500 window, both per-spaxel and PowerBin attempts
failed (Test 19 in nb05x — outer bins σ > 800 km/s, no per-spaxel S/N
above floor). At wR3800_5400_arcmask the same machinery works because
the per-pixel S/N quadruples (more line-strength signal per spaxel from
the broader feature set + arc-emission masking). This is a direct
qualitative gain enabled by the methodology change.

**Caches & figures (nb11):**
- `results/nb11_perbin_perspaxel_wide.npz` — per-bin σ/V/χ², per-spaxel
  σ/V/χ², bin_map, V_sys
- `results/figures/nb11_perspaxel_sn.png` — S/N map + S/N histogram
- `results/figures/nb11_powerbin_map.png` — PowerBin bin map + per-bin S/N
- `results/figures/nb11_powerbin_kinematics.png` — σ map, V map, σ(R), V(R)
- `results/figures/nb11_perbin_perspaxel_sigma_vs_r.png` — combined σ(R)

### L. Paper-figure preparation (2026-05-13)

Final figure cleanup for the ApJL submission, in
`../AGEL_0206_ApJL_Figures/figures.ipynb`.

**L1 Figure 2 narrow (cell `a546db7f`)** —
Updated title to show σ_e = 267 ± 24 (stat) ± 18 (sys) km/s
(stat ± sys separately rather than a single combined ±30). Removed the
"no arc mask" blue overlay from the inset; left inset now shows a
single red posterior. Saves to `AGEL0206_sigma_e_SED_final.pdf`.

**L2 Figure 2 wide (cell `fig2_wide`, NEW, paper headline)** —
σ_e = 255 ± 6 (stat) ± 17 (sys) km/s. I-weighted aperture spectrum over
rest 3800–5400 Å with 9 absorption features labeled (Ca K, Ca H, Hδ,
G-band, Hγ, Hβ, Mg b, Fe5270, Fe5335) and 4 source-emission arc masks
hatched. Residual sub-panel underneath. Inset shows single red posterior
(no-arc-mask comparison removed). Saves to
`AGEL0206_sigma_e_SED_final_wide.pdf`.

**L3 Figure 3 M_BH–M_star cleanup (cell `fig3_clean_no_err`, NEW)** —
Cleaned version of the comparison plot. NO per-object error bars, uniform
medium-gray color for all local early-types (no E vs S0 / KH13 vs
post-KH13 subdivision). Filtered to **Greene+2020 E+S0 list (60 unique
galaxies)**: dropped 19 post-KH13 additions that were not in Greene+2020
Fig 10 (mostly Saglia+2016 / Thater+2019 fills), matched 4 KH13 suffix
variants (NGC 1316b / 2778c / 3998c / 4486A) to their Greene+2020 names
(NGC 1316 / 2778 / 3998 / 4486A). All literature relations (Greene+20,
Reines & Volonteri+15, Farrah+25, Sahu+24) preserved. AGEL0206 point at
M★=10^11.3, M_BH=6.5e8 preserved. Saves to `Mbh_Mstar_relation_clean.pdf`.
The legacy cell 33 with per-object error bars and E vs S0 color subdivision
is preserved as `Mbh_Mstar_relation.pdf`.

**Sanity / withdrawn:** during the wide-window writeup an erroneous
"Cappellari 2017 recommends extending blueward" citation was drafted in
chat — **retracted**, no such recommendation exists in the ppxf paper.
The wavelength-range choice is data-driven (XSL template floor + retain
Ca H+K), with LEGA-C (van der Wel+21) as the closest literature
precedent (see J1 above).

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
| `notebooks/09_final_sigma_e_paper.ipynb` | **Paper-ready σ_e** at both windows; §9 three-window cross-check, §10 I-shape N=250, §10.1 mask-weight N=250 |
| `notebooks/09a_wavelength_window_sweep.ipynb` | 15-window scan; identifies wR3800_5400 sweet spot |
| `notebooks/09b_lsf_and_sigma_clip.ipynb` | Arc-mask + ppxf clean=True diagnostic with 3×3 zoom plots |
| `notebooks/09d_final_systematics_both_windows.ipynb` | **Side-by-side WIDE vs NARROW final budget** (paper headline writeup) |
| `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` | Per-spaxel + PowerBin σ map at wR3800_5400_arcmask (NEW — previously infeasible) |
| `notebooks/04_redshift_verification.ipynb` | z = 0.67564 line-fit |
| `notebooks/02_streamlined_Bagpipes_SED.ipynb` | log M★ = 11.33 |
| `notebooks/08_sersic_total_photometry.ipynb` | Sersic-total photometry cross-check |
| `notebooks/07c_sigma_e_equalN.ipynb` | §6cum + §7 equal-N annular cross-check (narrow window, retained for cross-check) |
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
| `scripts/run_isource_shape_sweep.py` | 10-shape I(r) sweep — **Sersic2D bound-fix applied 2026-05-11** (n ∈ [1.0, 6.0] + multi-init grid) |
| `scripts/analyze_isource_shape_sweep.py` | I-shape sweep analysis + figure |
| `scripts/run_window_sweep.py` | Wavelength-window driver with arc-mask flag; underpins nb09a + nb09 §9 |
| `scripts/soft_mask_track.py` | Mask-weight runs (CLI `--weight`) |
| `scripts/analyze_mask_weight_sweep.py` | 5-point sweep analysis + linearity fit |
| `scripts/regen_s6cum_nomask_diagnostic.py` | Centering investigation |
| `scripts/relock_nomask_Re_N500.py` | One-off N=500 relock for the no-mask track |

### Result files (key)

| File | Contents |
|---|---|
| **WIDE arc-masked (current headline)** | |
| `results/sigma_e_final_systematics_nb09d.npz` | **Both-windows side-by-side budget** (paper headline writeup) |
| `results/nb09a_wavelength_sweep/wR3800_5400_arcmask_*_T*_N500.npz` | N=500 stat pool at headline window |
| `results/nb09a_wavelength_sweep/wR4000_5400_arcmask_*_T*_N500.npz` | N=500 orthogonal-feature-set cross-check |
| `results/ishape_sweep_wR3800_5400_arcmask/` | 10 shapes × 3 SPS × N=250 (I-shape ±1.5) |
| `results/maskweight_sweep_wR3800_5400_arcmask/` | 3 weights × 3 SPS × N=250 (mask ±3.8) |
| `results/sigma_e_window_sweep_09a.npz` | 15-window posterior summary |
| `results/nb11_perbin_perspaxel_wide.npz` | Per-bin + per-spaxel kinematics at wide |
| **NARROW Ca H+K + G (cross-check)** | |
| `results/final_sigma_e_paper.npz` | nb09 narrow headline + per-SPS + 3-track summaries |
| `results/final_sigma_e_paper/` | Per-(aperture, SPS, mask_weight) caches |
| `results/sigma_e_radial_07c.npz` | nb07c §6cum + §7 + §7b posteriors |
| `results/annular_bootstrap_07c/` | Per-aperture cumulative I-weighted caches |
| `results/annular_bootstrap_07c_ishape/` | Per-(shape, SPS) I-shape caches (N=500, F140W/F200LP Sersic2D refit 2026-05-11) |
| `results/annular_bootstrap_07c_nomask/` | §6cum nomask caches |
| `results/mask_weight_sweep.npz` | 5-point mask-weight sweep summary (narrow) |
| **Audits & figures** | |
| `results/ppxf_methodology_audit.npz` | 4-audit results |
| `results/figures/nb09_*.png` | nb09 paper figures (incl. fit windows + I-shape comparison) |
| `results/figures/nb09d_per_sps_both_windows.png` | Per-SPS posterior collapse (26→4 km/s) |
| `results/figures/nb11_*.png` | Resolved-kinematics maps (per-spaxel S/N, σ map, V map, σ(R), V(R)) |
| `results/figures/nb07c_ishape_sweep.png` | I-shape sweep figure |
| `../AGEL_0206_ApJL_Figures/AGEL0206_sigma_e_SED_final_wide.pdf` | Paper Figure 2 (wide, headline) |
| `../AGEL_0206_ApJL_Figures/AGEL0206_sigma_e_SED_final.pdf` | Paper Figure 2 (narrow, cross-check) |
| `../AGEL_0206_ApJL_Figures/Mbh_Mstar_relation_clean.pdf` | Paper Figure 3 (clean, no error bars, Greene+2020-matched) |

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

**Updated 2026-05-13** with status changes after the wide arc-masked window
work and the Sersic2D bound-fix.

- **Re-run nb09 at higher N** (e.g. N=2000): the pooled posterior is already
  smooth; at the wide arc-masked window the stat error is ±6 km/s while
  the fit-window systematic dominates at ±15 km/s. More N is wasted compute.
- **~~Re-fit Sersic2D F200LP~~ DONE 2026-05-11** — the n=0.30 fit was a
  parameter-bound pathology, not a UV-arc artifact. Bound-fix (J7) tightened
  bounds to n ∈ [1.0, 6.0] with multi-init grid. F200LP_Sersic2D now
  converges to a physical fit (σ=251.4 km/s) and is **kept in** the I-shape
  budget, not flagged as an outlier.
- **Voronoi binning** — Test 18 / nb05x — abandoned; not pursued at N=500.
- **~~PowerBin spatial binning~~ PARTIAL** — Test 19 narrow-window rejected
  was correct at the time. **At the wide arc-masked window (K2 / nb11)
  PowerBin now works** for the inner ~5 of 7 bins; 2 outer bins still hit
  σ > 400 (irreducible KCWI per-spaxel S/N at aperture edge). Not used as
  headline binning scheme but published as supplementary σ map.
- **~~Per-spaxel rotation map at production statistics~~ DONE (K1 / nb11)**
  — at the wide arc-masked window, S/N ≥ 5 gives 17 central spaxels with
  clean σ values and a visible σ(R) gradient. No coherent rotation detected
  above the noise (consistent with elliptical-bulge expectation). Published
  as supplementary kinematic map.
- **Full nb07a Sersic-I revisit at the post-frame-fix headline** —
  superseded by the I-shape sweep (D3/D4/J6) which covers the Sersic2D
  shape with the bound-fix applied.
- **Re-quote σ_e on KH13/Bell+2003 stellar masses for the M_BH–M★ plot** —
  Figure 3 keeps the KH13 dual-M/L dynamical bulge masses for the local
  comparison sample. The Greene+2020 list match (L3) restricts to the same
  60 galaxies as Greene+2020 Fig 10 but with our (better) bulge mass values.
- **Cross-match nb09a "narrow" window to other AGEL targets** — paper is
  single-target. The arc-emission-mask methodology (J2/J3) is reusable for
  future AGEL targets where source z is known, but applying it elsewhere
  is not in scope for the ApJL.

---

*See `~/.claude/projects/.../memory/MEMORY.md` for compressed cross-session
references, including:*
- *`reference_cumulative_vs_annular_sigma_e.md` — §6cum vs §7 design choice*
- *`reference_sps_systematic.md` — per-SPS V_sys + frame-fix details*
- *`reference_ppxf_vacair_handling.md` — vac/air conversion methodology*
- *`reference_kcwi_data_properties.md` — KCWI cube + LSF + paper boilerplate*
- *`project_nb09_final_sigma_e.md` — current headline status*
