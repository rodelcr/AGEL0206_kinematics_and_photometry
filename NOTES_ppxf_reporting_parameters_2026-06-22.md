# ppxf reporting parameters — stock-take vs Looser et al. 2024 (2026-06-22)

Goal: match the ppxf-reporting specificity of **Looser et al. 2024, Nature 629, 53**
(DOI 10.1038/s41586-024-07227-0; Zotero USS66QPG), so our methods section has every number on hand.

**Key difference in PURPOSE:** Looser use pPXF as one of four SED/SFH codes on a low-res (R100) prism
and **FIX** the velocity dispersion to a virial estimate (σ ≈ √(GM⋆/5R_e) = 50 km/s) because the prism
can't resolve it. **We use ppxf to MEASURE** the stellar velocity dispersion σ_e (our headline), on
R≈3600 KCWI spectra. So our SSP-grid/SFH details are lighter (kinematics, not SFH — SFH/M⋆ is our
Bagpipes job), but our kinematic-fit details (window, degree, masking, bootstrap, multi-SPS pooling)
are heavier. Below: what Looser report ▸ our value (✓ on hand / ⚠ surface in METHODS).

## What Looser et al. 2024 report for pPXF (their bar)
- Code + refs: pPXF (Cappellari 2017; Cappellari 2023), χ²-minimization.
- Templates: C3K atmospheres + MIST isochrones, solar abundances; **62 ages × 10 metallicities**;
  age 10^6.0–10^9.2 yr; log(Z/Z⊙) = −2.5…+0.5.
- Velocity dispersion: **fixed** at virial 50 km/s (prism too low-res to fit it).
- Resolution: SSP convolved to the wavelength-dependent prism LSF.
- Normalization: spectrum + templates renormalized by median flux per pixel.
- Outlier rejection: initial fit + σ-clip.
- Errors: residual-based bootstrap of the best fit, **without regularization, 1,000 iterations**;
  bootstrapped SSP grids averaged → non-parametric SFH.
- Dust: SSP multiplicatively coupled to an attenuation curve; A_V = 0.4 ± 0.1 (note: driven by UV slope).
- Outputs: age–metallicity weight grid, mass-weighted fractions, M⋆, SFR, Z, t_quench, t_form, A_V.
- Data/code availability statement (GitHub + MAST).

## Our ppxf kinematics setup — every number on hand

**A. Code & references** ✓
- **ppxf v9.4.5** (`ppxf.__version__`). Integrated-spectrum variant uses **veldis** (degree=[4,30]).
- Refs: Cappellari & Emsellem 2004 (PASP 116, 138); Cappellari 2017 (MNRAS 466, 798); Cappellari 2023
  (MNRAS 526, 3273). σ_e definition: Gültekin et al. 2009 eq. 1; aperture co-add: Cappellari et al. 2006.

**B. Template libraries — THREE, pooled** ✓ (ppxf SPS data files v9.0, `spectra_{sps}_9.0.npz`)
| library | ref | N_age × N_Z (=N_templ) | age range (Gyr) | log(Z/Z⊙) range | native frame | native FWHM (median, full range) |
|---|---|---|---|---|---|---|
| **FSPS** | Conroy, Gunn & White 2009 | 43 × 9 = **387** | 0.001–15.8 | −1.75…+0.25 | **vacuum** | 2.51 Å |
| **EMILES** | Vazdekis et al. 2016 | 25 × 6 = **150** | 0.063–15.8 | −1.71…+0.22 | **air** | 4.94 Å |
| **XSL** | Verro et al. 2022 | 26 × 8 = **208** | 0.050–15.8 | −2.20…+0.20 | **air** | 0.80 Å |
- (Median FWHM is over the FULL stored range; in our fit window rest 3800–5400 Å the optical resolution
  is finer and is matched per-wavelength by `sps_lib` — see D. Template grids restricted to the fit range
  via `lam_range` in `sps_lib`.)

**C. Velocity dispersion — MEASURED (our headline)** ✓
- **σ_e(<R_e) = 267.31 −11.98/+11.58 km/s** (stat −5.0/+4.1 ⊕ sys ±10.87; sym ±11.77), at R_e=2.097″.
- `moments=2` (fit V and σ only; no h3/h4 → ppxf BIAS penalization irrelevant, default). `trig=False`.
- Per-SPS V_sys subtracted before pooling (frame-aware; collapses V_sys split-track ~110→~15 km/s).
- 7-component systematic budget (stat ⊕ I-shape ⊕ arc-mask ⊕ centering ⊕ fit-window ⊕ reduction-pass ⊕
  R_e-source) — all in `results/PAPER_VALUES.json` / CLAUDE.md.

**D. Spectral resolution / LSF** ✓
- KCWI `DISPSCAL = 0.294` → **FWHM_inst = 0.692 Å** (constant in observed Å); σ_v,inst ≈ 12–14 km/s.
- Templates convolved to the data resolution via `fwhm_dict = {lam, fwhm}` passed to `sps_lib`, in each
  library's native frame. LSF audited: TESTS B1 (FWHM=0.692 Å) + B2 (×0.5–2.0 LSF → max |Δσ| 0.83 km/s).

**E. Wavelength frame (air/vacuum)** ✓
- KCWI native = **vacuum**. Galaxy expressed in each library's native frame:
  `SPS_NATIVE_FRAME = {fsps: vacuum, emiles: air, xsl: air}`. Frame test D4 (<0.5 km/s σ impact).

**F. Fit window & sampling** ✓
- Headline: **rest 3800–5400 Å** (`wR3800_5400_arcmask`; obs = 3800/5400 × (1+0.67564)). **2161 good pixels.**
- velscale set from the data log-step (C_KMS · Δln λ). z_systemic = 0.67564.
- Cross-check windows: `wR4000_5400_arcmask` (Hβ+Mg b, no Ca H+K); narrow Ca H+K+G `w6500_7500`.

**G. Polynomials** ✓
- Additive degree **swept 15–29 (15 degrees), pooled**; **mdegree = 0** (no multiplicative polynomial).

**H. Masking / goodpixels** ✓
- No-Balmer mask (Hδ, Hγ, Hβ kept in fit); **35-entry BAD_PIXELS_REST** (26 CR + 9 OH/sky, M10 audit);
  **ARC_MASKS_REST** source-emission catalog (Mg II 2796/2803, [O II] 3727/29, **He I 3819**, [Ne III] 3869,
  mapped to deflector rest frame); spatial arc mask (F200LP-located, reprojected to IFU grid).

**I. Normalization** ✓ — ppxf-internal (galaxy normalized to median); templates normalized by `sps_lib`.

**J. Error method** ✓ (heavier than Looser)
- **Wild bootstrap**: hybrid Rademacher sign-flip × local-residual scaling (rolling 75-pix window).
- **N = 500** production per (SPS, degree) [N=50 smoke]. **Pooled** across 3 SPS × 15 degrees →
  ~22,500 samples → asymmetric 16/50/84 percentiles. This pooled width marginalizes the SPS-library +
  polynomial-degree systematics (between-SPS ±2.04 ⊕ within-SPS ±4.22 ≈ pooled ±4.64).
- `regul` not used (**regul = 0**; kinematics, no SFH regularization) — same "without regularization"
  spirit as Looser, but for σ not SFH.

**K. Aperture / extraction** ✓
- I-band-weighted **R < R_e aperture spectrum** (Cappellari+2006 co-add; Gültekin 2009 eq. 1 σ_e
  definition). Cross-checks: §6cum cumulative I-weighted, §7 discrete annular.

**L. Production drivers / data-code availability** ✓
- `scripts/bootstrap_ppxf.py` (frame-aware prep + wild bootstrap), `scripts/final_sigma_e.py` +
  `scripts/run_wide_sigma_e.py` (3-SPS × 15-degree × N=500 pooled driver). Cube:
  `raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`.

## To surface explicitly in the METHODS paragraph (we have all of these; just write them down)
1. **ppxf v9.4.5** + the three SPS refs (currently only the SPS names are in CLAUDE.md).
2. The **per-library grid dimensions / age & Z ranges / native frame** table above (NEW — was not written
   down anywhere before today; now captured here).
3. **moments=2, mdegree=0, regul=0, degree range 15–29, trig=False** in one sentence.
4. **2161 good pixels**, **N=500 × 3 SPS × 15 degrees = ~22,500 pooled bootstrap samples**.
5. FWHM_inst = 0.692 Å and the template→data convolution (already in METHODS/TESTS B1–B2; cross-cite).

Nothing is missing to reach Looser's specificity — items 1–5 are wording tasks, and the one genuine
data gap (the SPS grid table) is now recorded here.
