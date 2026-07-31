# HANDOFF — structural nature, relation-offset provenance, ppxf stock-take, tutorials, scatter-band figures (2026-07-06)

Supersedes `HANDOFF_mstar_cored_total_light_2026-06-17.md`. Everything below is on disk + self-consistent;
`scripts/paper_values.py --check` is GREEN on all 5 docs. **Single source of truth = `results/PAPER_VALUES.json`**
(written by `scripts/paper_values.py`). NEVER hand-type headline numbers.

## TL;DR — headline (unchanged this arc; Planck 2015, 7.2764 kpc/arcsec at z_l)
- **σ_e(<R_e) = 267.31 −11.98/+11.58 km/s** (267 +4/−5 stat ± 11 sys). ppxf **v9.4.5**, wild bootstrap N=500 × 3 SPS × 15 deg pooled.
- **log(M⋆/M☉) = 11.46 +0.07/−0.06 (stat) ± 0.17 (sys)** → reported **11.5 ± 0.1 (stat) ± 0.2 (sys)** at 1 dp.
- **M_BH = 5.2 +1.7/−1.4 ×10⁸ M☉** (log 8.72) — **PROVISIONAL** (Ferrami free-BH 1σ envelope).
- **R_e = 2.097″ = 15.26 kpc**. z_l = 0.67564, z_s = 1.30263. Passive SED: age 3.8 Gyr, A_V 0.21, sSFR −10.5.
- **GUIDING PRINCIPLE (keep):** independent probe — never choose estimators by relation-consistency.

## What changed this arc

### 1. M⋆ robustness + quiescent-prior write-up
- **Drop-F200LP test (A5g, `scripts/mstar_drop_f200lp.py` → `results/mstar_drop_f200lp.npz`):** refit on F140W+F150W2+F322W2 only. Quiescent prior → log M⋆ = **11.56** (Planck), +0.10 dex vs 11.46 (within ±0.17 sys) → **robust to dropping the bluest band**. Fiducial prior → 11.29 (SFR≈56) = the young+dusty branch collapses harder without F200LP, reconfirming the spectrum (not any single band) selects the solution.
- **Quiescent-prior definition** documented (DRAFTING §3.2 table + TESTS A5h): same exponential-τ model, 3 bounds tightened vs fiducial — **age (0.1,15)→(4,15) Gyr, tau (0.3,10)→(0.1,1.5) Gyr, Av (0,2)→(0,0.6) mag**; massformed/metallicity/Calzetti/z unchanged. Motivated by the **independent** KCWI absorption spectrum (not circular). Fiducial↔quiescent spread = SFH-prior systematic (0.08 dex).

### 2. PSF models + core-Sérsic (structural)
- **PSF models built, all 4 bands, in-env (A6, `scripts/psf_star_census.py`, `scripts/build_psf_models.py` → `results/psf_models/<band>_psf.npz`):** env has NO webbpsf/grizli (TinyTim uncompiled). **F140W** = STScI `PSFSTD_WFC3IR_F140W.fits`; **F200LP/F150W2/F322W2** = empirical EPSF from field stars (HST full-frame drc/drz 146/102 stars; JWST i2d 3 SW / 13 LW; the 25″ lens cutout has 0). FWHM 0.103/0.085/0.049/0.112″ — match instrument. F200LP needs ov=2 + quadratic kernel (UVIS undersampled).
- **Core-Sérsic vs single-Sérsic (A7, `scripts/core_sersic_test.py` → `results/core_sersic_test.npz`):** PSF-convolved on F150W2+F140W. Free core-Sérsic "prefers" ΔBIC 386/43 **but it's a global-fit artifact** — r_b=2.5–3.6 kpc (10–50× too big for a depleted core), no inner residual deepening. **Core UNRESOLVED; r_b < 0.024″ = 0.17 kpc; <1% light → negligible for M⋆.** Report: `REPORT_cored_vs_coreless_2026-06-17.md`.
- Notebook `notebooks/18_cored_vs_coreless_sersic.ipynb` synthesis updated.

### 3. DRAFTING §3.4 — "Structural nature of the deflector"
New section: **(a)** massive quiescent early-type (E/S0 not cleanly separated); **(b) classical bulge/elliptical, NOT pseudobulge** (mass 10¹¹·⁵, σ 267, **V/σ ≈ 0.15** dispersion-support from `radial_sigma_combined_posterior.npz` V≈40, passivity; the low single-Sérsic n 1.2–1.6 is a box-limited artifact) → valid clean M•–σ probe; **(c) cored/coreless indeterminate** (A7). Numbers fact-checked (see §5).

### 4. Relation-offset significance + provenance (O1)
- **`scripts/relation_offset_significance.py` → `results/relation_offset_significance.npz`:** offset + pythagorean N_σ vs Greene+2020. **Distance to mean: −0.26 (M•–σ) / −0.49 (M•–M⋆) dex** (=1.8σ measurement-only); **N_σ incl. intrinsic scatter = 0.57σ / 0.69σ → fully consistent.** DRAFTING §3.3.5a.
- **Provenance VERIFIED** (`NOTES_relation_offset_provenance_2026-06-17.md` + `.pdf`): Greene+2020 Suppl. Table 5 **"Early"** — M•–σ **8.03±0.06 / 4.24±0.25 / ε0.43±0.04** (pivot 160); M•–M⋆ **7.89±0.09 / 1.33±0.12 / ε0.65±0.05** (pivot 3e10). Confirmed vs in-repo `../AGEL_Mbh_sigma/Greene20_Supple_5.csv` + `aa58_greene_supmat.pdf`. Methodology cross-checked vs Pacucci+2023/Pacucci&Loeb2024 + Melo-Carneiro+2025 (field-standard: offset ÷ combined intrinsic-scatter+measurement).
- **Greene sample max redshift (§5, `scripts/greene_sample_redshift.py`):** A1836-BCG 152.4 Mpc → **z≈0.034** for both samples; AGEL0206 z=0.676 is ~20× higher.
- ⚠️ **Registry `mbh_sigma_offset` = −0.69 is WRONG** (added ε to the mean); correct mean-offset = −0.26. Per user the figure band stays; registry reconciliation (−0.69→−0.26, add `mbh_mstar_offset`=−0.49, add n_sigma leaves) **deferred pending authorization**.

### 5. fact-critic + Zotero curation on §3.4 (2026-07-02)
- **`FACT_CRITIC_DRAFTING_FACTS_paper_2026-05-29_2026-07-02.md`:** 15/17 numbers PASS vs registry/npz; **2 auto-fixed** in §3.4 (F140W box n "0.9→2.0"→**1.2→2.0**; "0.02–0.1 kpc = a few–40 mas"→**≈3–14 mas**).
- **Curation:** added 5 CrossRef-verified refs to Zotero **"SMBH from Lensing" (RQSB3CH9)**: Faber+1997 (10.1086/118606), **Lauer+2007 cusp/core ApJ 664,226 (10.1086/519229) — distinct from the selection-bias Lauer+2007 already in library (3TUYW9FH)**, Kormendy+2009 (10.1088/0067-0049/182/1/216), Fisher&Drory 2008 (10.1088/0004-6256/136/2/773), Fisher&Drory 2016 (10.1007/978-3-319-19378-6_3). §3.4 prose now cites these; Lauer disambiguation noted inline.

### 6. ppxf reporting stock-take vs Looser et al. 2024 (`NOTES_ppxf_reporting_parameters_2026-06-22.md`)
Every number on hand at Looser's specificity. Newly captured (was undocumented): **ppxf v9.4.5**; SPS grids **FSPS 43×9=387 (vac), EMILES 25×6=150 (air), XSL 26×8=208 (air)**, age/Z ranges, native FWHM; call params **moments=2, mdegree=0, regul=0, trig=False, degree 15–29, 2161 good pixels, FWHM_inst 0.692 Å**. (Their ppxf FIXES σ for SFH; ours MEASURES it.)

### 7. ppxf/redshift TUTORIALS (01–04) — headline alignment + audit + Drive
- z → **0.67564** across 01–04; **SPS dispersion now POOLED** in 02 §10 (concatenate FSPS+EMILES+XSL bootstrap → one posterior = recommended error; removed the old WRONG "don't add template spread" line); chi2 prose fixed (box actually runs **chi2/DOF≈0.04**, noise over-estimated); RuntimeWarnings suppressed. Independent audit found + fixed 2 BLOCKERs. All 4 re-executed clean.
- **Drive updated:** `My Drive/AGEL/ppxf_redshift_tutorials/` (ready-to-run `ppxf_tutorials/` folder incl. cube + `setup.sh` one-command + `START_HERE.txt/.md`, + `ppxf_tutorials_code.zip`). Md5-identical local↔Drive. See `[[project_ppxf_tutorial_suite]]`.

### 8. NEW ApJL figures — intrinsic-scatter Figs 3 & 4
New section at the bottom of `AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb` (cells ~24–26). Copies of the final M•–σ (cell 16) and M•–M⋆-clean (cell 21) but the Greene band is the **true ±1σ INTRINSIC SCATTER** (0.43 / 0.65 dex) not the mean-locus fit-uncertainty; **solid `#eb4934` central line inside the band**; both merged into **one legend entry** (custom `HandlerBandLine` = band rectangle + line through centre). Saved as **`Mbh_sigma_relation_scatterband.pdf`** + **`Mbh_Mstar_relation_clean_scatterband.pdf`**. Originals (`*_relation.pdf`, `*_clean.pdf`) untouched. **Render recipe:** the two new cells only need setup cells [1,10,11,13,14,19]; run in env **`agel0206_figs`** (a throwaway subset notebook avoids re-running Fig1/2 and the random-scatter cells).

## Outstanding TODOs
- **Registry reconcile (needs user OK):** `mbh_sigma_offset` −0.69 → mean −0.26; add `mbh_mstar_offset` −0.49 + `n_sigma` leaves.
- **M_BH PROVISIONAL** until Ferrami resolves free-vs-fixed. `DLOGM_COSMO` cleanup (regen `aperture_matched_photometry.npz` natively under Planck).
- §3.4/REPORT expected-core-radius: reconcile 0.02–0.1 kpc vs 0.005–0.05″ to one cited M•–r_b scaling (Lauer+2007/Rusli+2013) — both below PSF so conclusion unaffected.
- DRAFTING §2.4.2 pipeline walkthrough still pre-audit.

## Resume / verify
```bash
conda activate ISMgas; cd ~/Documents/AGEL/AGEL0206_kinematics_and_photometry
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md TESTS_AND_DIAGNOSTICS.md METHODS_AND_SYSTEMATICS.md PROJECT_BRIEF.md   # GREEN
python -m scripts.relation_offset_significance    # 0.57σ / 0.69σ
python -m scripts.greene_sample_redshift          # A1836-BCG z≈0.034
python -m scripts.core_sersic_test                # core unresolved
# figures: AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb (env agel0206_figs); new cells at bottom
```
Provenance PDFs: `NOTES_relation_offset_provenance_2026-06-17.pdf`. All 5 state docs linter-GREEN & headline-consistent.
