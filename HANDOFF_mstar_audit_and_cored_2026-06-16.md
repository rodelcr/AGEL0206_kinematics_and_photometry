# HANDOFF — M⋆ audit (apcorr + SED prior), M_BH uncertainty, figs/tables/docs, cored investigation — 2026-06-16

Resume point after a long arc. Everything below is on disk and self-consistent; `paper_values.py --check` is GREEN.
Single source of truth = `results/PAPER_VALUES.json` (written by `scripts/paper_values.py`). NEVER hand-type headline numbers.

## TL;DR — current headline (Planck 2015, H₀=67.7, Ω_m=0.302; 7.2764 kpc/arcsec at z_l)
- **σ_e(<R_e) = 267.31 −11.98/+11.58 km/s** (sym ±11.77; stat −5.0/+4.0 ⊕ sys ±10.87). UNCHANGED this arc.
- **log(M⋆/M☉) = 11.46 +0.07/−0.06 (stat) ± 0.17 (sys)** → reported one-decimal **11.5 ± 0.1 (stat) ± 0.2 (sys)**.
  sys = masking 0.086 ⊕ apcorr-model 0.118 ⊕ SFH-prior 0.080. **(Was 11.50 pre-audit.)**
- **R_e = 2.097″ = 15.26 kpc** (method sys ±0.10″). UNCHANGED.
- **M_BH = 5.2 +1.7/−1.4 ×10⁸ M☉** (log 8.72 +0.12/−0.14) — **PROVISIONAL** (Ferrami). M•–σ offset −0.69 +0.12/−0.14 dex (vs Greene+2020).
- **Passive SED (quiescent fit):** mass-weighted age 3.8 Gyr, SFR 8 M☉/yr, sSFR −10.5, Z 1.25 Z☉, A_V 0.21.
- z_l=0.67564, z_s=1.30263.

## What changed this arc (2026-06-16)
1. **figures_paper4 re-run + cleaned** (`AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb`): inline outputs synced to registry; cell[8] stale σ_e comment block removed; legacy `figures.ipynb` marked **SUPERSEDED** (banner). 4 ApJL PDFs regenerated.
2. **M⋆ AUDIT (two coupled fixes; supersedes 11.50 → 11.46):**
   - **(a) Validated-fit aperture correction** (`scripts/aperture_correction_validated.py`): the auto `fit_sersic`
     (photometry_systematics.py, R_E_INIT=2.305) was biased (r_eff 2.2–2.7″, ellip→0 vs photutils-validated
     1.6–2.0″ & CoG R_e=2.097″) → over-corrected beyond-aperture light ~0.15 mag/band → old total inflated to 11.50.
     FIX: Sérsic SHAPE fixed to the photutils-validated per-band table fit, amplitude+sky fit to data; auto mask kept for F_raw.
   - **(b) Quiescence-constrained SED prior** (`scripts/mstar_headline_quiescent.py`): the 4-band SED is age–dust–M/L
     DEGENERATE (ΔlnZ=−0.18, identical χ²; `scripts/sed_quiescence_check.py`) and the flat-age prior slid onto a
     spurious young+dusty branch (SFR~57). The KCWI absorption-line spectrum (no emission) forces the passive solution.
     Headline adopts an old/short-τ/low-dust prior; fiducial↔quiescent spread carried as the SFH-prior systematic.
   - Net: photometry −0.11 dex, old-population M/L +0.10 dex → **11.46**. Ladder now MONOTONIC: raw 11.28 < rawapc 11.34
     < filled 11.42 < total 11.46 < sersic 11.49 (old auto had total 11.50 *above* sersic 11.43 — the inversion that flagged it).
3. **M_BH uncertainty fixed** (`scripts/bh_mass_combine.py`): iterated ±1.4e8 (within-free mixture-1σ, too tight) →
   free+fixed envelope 3.8e8–2.3e9 (too wide; fixed-BH not in final selection) → **ADOPTED = free-BH selected-EPL
   1σ model envelope (min med−σ to max med+σ over the 3 free-BH models) = 3.8–6.9e8 → 5.2 +1.7/−1.4 ×10⁸.**
4. **Fig 2** (cell 8): M⋆ inset now shows **stat ± sys** (mirrors σ_e); SED panel (spectrum+photometry+M⋆) now from the
   validated+quiescent headline run (`results/mstar_headline_quiescent.npz`), superseding the old `bagpipes_sed_results.npz`.
   **Fig 3/4** (cells 16/19/20): M_BH & M⋆ point read from registry (joint-1σ M_BH; M⋆ stat⊕sys); hardcoded values/stale comments removed.
5. **key_results_table.tex** → plain **single-column `table`/`tabular`** (NOT deluxetable; user: it spanned 2 columns),
   **no tablenotes** (user deleted them — do not re-add), σ_e integer, M⋆ one-decimal stat±sys. Compiles in aastex631.
   Generator: `scripts/make_latex_tables.py`. Offset row dropped per user slim edit.
6. **Docs synced + PDFs rebuilt:** DRAFTING §3.2 (final M⋆ procedure + budget, auto-rendered `PV:auto:mstar_budget`),
   §3.3.1/§3.3.5 (M_BH envelope), low-tail bullet; CLAUDE/METHODS/PROJECT_BRIEF headlines; TESTS L3 value + L0 supersede
   note. **PDFs rebuilt** via `./build_docs_pdf.sh` (pandoc in /opt/anaconda3/bin; needs it on PATH). Glyph fixes:
   `𝓔`→`Z` (evidence), `⊙`(U+2299)→`☉`(U+2609). `paper_values.py --check` GREEN on all 5 docs.
7. **Cored-vs-coreless investigation STARTED** — see next section.

## Cored-vs-coreless investigation (NEW, in progress)
**Notebook:** `notebooks/18_cored_vs_coreless_sersic.ipynb` (ISMgas; runs clean end-to-end). Figs:
`results/figures/cored_test_{A_boxsweep,B_profile}.png`.
- **Physics:** σ=267 (>>240 divide) → EXPECTED to be a **core** elliptical (Faber97/Lauer07/Kormendy09), BCG-like
  (~100 kpc from ACT-CL J0206). BUT expected core r_b~0.1–0.3 kpc = 0.014–0.041″ < every PSF (JWST SW 0.05″, F140W 0.13″)
  → **core UNRESOLVED** at z=0.68; <1% of light; M_BH–core-scouring cross-check NOT feasible.
- **The real bite = Test A (box-size/n):** single-Sérsic n & total light KEEP RISING with fit box and don't plateau
  (F140W: n 0.9→2.0, total **+0.73 mag** box 4→10″ and still rising; F322W2 +0.43 mag). Classic **box-limited** cD
  signature → the single-Sérsic total is **non-convergent**, so our M⋆ (4″-box-based) is a **LOWER BOUND** (could be
  ~+0.3 dex). F200LP (faint rest-UV) & F150W2 (huge mosaic, field-contaminated) rail at large box → flagged/excluded.
- **Implication:** the headline M⋆ total-light is **method-sensitive for this cD** (galaxy↔ICL). The follow-up is the
  TOTAL-LIGHT METHOD, not the core (which is unresolved).

## Files modified / created this arc
- NEW scripts: `aperture_correction_validated.py`, `sed_quiescence_check.py`, `mstar_headline_quiescent.py`,
  `sed_properties_final.py`. Edited: `paper_values.py` (M⋆ from quiescent run + 3-comp sys + sed block + M_BH envelope +
  render_mstar_budget), `make_latex_tables.py` (plain table), `bh_mass_combine.py` (envelope), `bagpipes_sersic_refit.py`
  (fit_and_extract returns mass_weighted_age).
- NEW notebook: `notebooks/18_cored_vs_coreless_sersic.ipynb`.
- NEW results: `mstar_headline_quiescent.npz`, `aperture_correction_validated.npz`, `sed_properties_final.npz`,
  `figures/cored_test_*.png`. Updated: `PAPER_VALUES.json`, `key_results_table.tex`.
- Docs: DRAFTING_FACTS/CLAUDE/METHODS/PROJECT_BRIEF/TESTS (+ rebuilt PDFs). Figures notebook + 4 ApJL PDFs + figures.ipynb banner.

## Outstanding TODOs (tasks #9–#11 + carry-overs)
- **#11 (BIG):** decide the M⋆ **total-light method** given Test A (single-Sérsic underestimates). Options: curve-of-growth
  total, larger fixed metric aperture, or bulge+envelope/ICL-aware decomposition. If M⋆ revises, cascade via registry.
- **#10:** PSF model setup (webbpsf JWST / TinyTim or grizli HST / empirical from field stars) — prerequisite for any
  rigorous inner-profile/core-Sérsic work. Env (ISMgas) lacks webbpsf.
- **#9:** core-Sérsic vs single-Sérsic (likely inconclusive — core unresolved; needs PSF).
- **External:** M_BH PROVISIONAL until Ferrami's final free-vs-fixed resolution (fixed-BH ~2e9 favored by his prose — UNRESOLVED).
  `DLOGM_COSMO`: regen `aperture_matched_photometry.npz` natively under Planck (drop +0.0282). DRAFTING §2.4.2 pipeline
  walkthrough still pre-audit (only §3.2 rewritten). METHODS `″` monospace-glyph warning (minor, pre-existing).

## How to resume / verify
```bash
conda activate ISMgas
cd ~/Documents/AGEL/AGEL0206_kinematics_and_photometry
python scripts/paper_values.py                         # regen registry
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md TESTS_AND_DIAGNOSTICS.md \
    METHODS_AND_SYSTEMATICS.md PROJECT_BRIEF.md         # drift linter → GREEN
python scripts/make_latex_tables.py                    # regen key_results_table.tex
PATH="/opt/anaconda3/bin:/Library/TeX/texbin:$PATH" ./build_docs_pdf.sh DRAFTING_FACTS_paper_2026-05-29.md TESTS_AND_DIAGNOSTICS.md
# figures: AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb (agel0206_figs kernel; Fig1 aplpy needs lens env)
# cored notebook: notebooks/18_cored_vs_coreless_sersic.ipynb (run from notebooks/, it chdir's to repo root)
```
First action next session: **the M⋆ total-light method decision (#11)** — Test A says the single-Sérsic total is a
lower bound; build a curve-of-growth (or ICL-aware) total and see if the headline M⋆ moves.
