# HANDOFF (comprehensive) — M⋆ audit, M_BH envelope, cored/total-light investigation — 2026-06-17

Supersedes `HANDOFF_mstar_audit_and_cored_2026-06-16.md`. Everything below is on disk + self-consistent;
`paper_values.py --check` is GREEN on all 5 docs. **Single source of truth = `results/PAPER_VALUES.json`**
(written by `scripts/paper_values.py`) — NEVER hand-type headline numbers; figures/tables/docs read the registry.

## TL;DR — current headline (Planck 2015, H₀=67.7, Ω_m=0.302; 7.2764 kpc/arcsec at z_l)
- **σ_e(<R_e) = 267.31 −11.98/+11.58 km/s** (sym ±11.77; stat −5.0/+4.0 ⊕ sys ±10.87). Reported integer in figs/table: **267 +4/−5 (stat) ± 11 (sys)**.
- **log(M⋆/M☉) = 11.46 +0.07/−0.06 (stat) ± 0.17 (sys)** → reported one-decimal **11.5 ± 0.1 (stat) ± 0.2 (sys)**.
  sys = masking 0.086 ⊕ apcorr-model 0.118 ⊕ SFH-prior 0.083. **CoG total-light cross-check = 11.69 (NOT folded — see below).**
- **R_e = 2.097″ = 15.26 kpc** (method sys ±0.10″).
- **M_BH = 5.2 +1.7/−1.4 ×10⁸ M☉** (log 8.72 +0.12/−0.14) — **PROVISIONAL** (Ferrami). M•–σ offset −0.69 +0.12/−0.14 dex (vs Greene+2020).
- **Passive SED (quiescent fit):** mass-weighted age 3.8 Gyr, SFR 8 M☉/yr, sSFR −10.5, Z 1.25 Z☉, A_V 0.21.
- z_l=0.67564, z_s=1.30263.

## GUIDING PRINCIPLE (set this arc — keep it)
**We are building an INDEPENDENT probe of the M•–σ / M•–M⋆ relations.** Do NOT choose any estimator
(M⋆, M_BH, …) because it is closer to / "method-matched" to the local relations — that is circular and
engineers the result. The relations are guide-the-eye only. Estimator choices are made on
**measurement-quality grounds**, and the reported errors must honestly span the systematics so the
placement on the relations is a genuine test.

## What changed this arc (2026-06-16 → 06-17)
1. **M⋆ AUDIT (supersedes 11.50 → 11.46), two coupled fixes:**
   - **Validated-fit aperture correction** (`scripts/aperture_correction_validated.py`): the auto `fit_sersic`
     was biased (r_eff 2.2–2.7″, ellip→0 vs photutils-validated 1.6–2.0″) → over-corrected beyond-aperture
     light ~0.15 mag/band. FIX: Sérsic SHAPE fixed to the validated per-band table fit, amplitude+sky to data.
   - **Quiescence-constrained SED prior** (`scripts/mstar_headline_quiescent.py`): 4-band SED is age–dust–M/L
     DEGENERATE (ΔlnZ=−0.18; `scripts/sed_quiescence_check.py`); flat-age prior slid to a spurious young+dusty
     branch (SFR~57) the KCWI absorption-line spectrum rules out. Headline uses old/short-τ/low-dust prior;
     fiducial↔quiescent spread = SFH-prior systematic. Net: −0.11 (photom) +0.10 (M/L) → 11.46. Ladder monotonic.
2. **M_BH uncertainty = free-BH (selected EPL) 1σ model envelope** (`scripts/bh_mass_combine.py`): min(med−σ)
   to max(med+σ) over the 3 free-BH EPL models = 3.8–6.9e8 → **5.2 +1.7/−1.4 ×10⁸**. (Iterated: ±1.4 mixture-1σ
   too tight → free+fixed 3.8e8–2.3e9 too wide → free-BH 2σ 2.6–8.3e8 too wide → adopted free-BH **1σ** envelope.
   Fixed-BH runs ~2e9 are the lower-evidence, UNRESOLVED alternative, NOT in the error.)
3. **Cored vs coreless + TOTAL-LIGHT investigation** (`notebooks/18_cored_vs_coreless_sersic.ipynb`):
   - σ=267 ⇒ physically a **core** elliptical, but core r_b~0.01–0.04″ < every PSF → **unresolved**; <1% of light;
     M_BH–core-scouring cross-check not feasible. Single-Sérsic is the pragmatic *inner* model.
   - **Test A (box-size/n):** single-Sérsic n & total light keep rising with fit box, don't plateau →
     **box-limited / non-convergent** (cD envelope). F140W n 0.9→2.0, total +0.73 mag (box 4→10″, still rising).
   - **Test A2 (curve-of-growth total, `scripts/cog_total_light.py`):** empirical CoG total = **log M⋆ 11.69**
     (CoG@8″ mags → quiescent Bagpipes), +0.22 dex above the single-Sérsic. BUT it lands only **+0.04 dex
     beyond the existing +1σ upper bound (11.65)** → **the total-light/cD-envelope ambiguity is ALREADY
     captured** by the systematic budget (the apcorr-model term spans the more-outer-light direction).
   - **DECISION (2026-06-17):** central M⋆ stays single-Sérsic **11.46** (well-defined/reproducible); CoG (11.69)
     recorded as a **cross-check ONLY** — NOT a separate systematic (a +0.22 one-sided term would double-count
     the apcorr-model contribution). CoG also non-convergent (sky/ICL) + strongly bandpass-dependent.
4. **Cascade + hygiene:** Fig 2 inset shows stat±sys (M⋆ & σ_e); Figs 3/4 read registry (M_BH joint-1σ, M⋆ stat⊕sys).
   `key_results_table.tex` = plain single-column `table`/`tabular` (NOT deluxetable; no tablenotes — user removed them,
   do NOT re-add), σ_e integer, M⋆ one-decimal. Docs synced; **glyph fixes** `𝓔`→`Z`, `⊙`(U+2299)→`☉`(U+2609).
   PDFs rebuilt (DRAFTING/TESTS/METHODS) via `./build_docs_pdf.sh` (needs `/opt/anaconda3/bin` pandoc on PATH).

## Files this arc
- NEW scripts: `aperture_correction_validated.py`, `sed_quiescence_check.py`, `mstar_headline_quiescent.py`,
  `sed_properties_final.py`, `cog_total_light.py`. Edited: `paper_values.py`, `make_latex_tables.py`,
  `bh_mass_combine.py`, `bagpipes_sersic_refit.py`.
- NEW notebook: `notebooks/18_cored_vs_coreless_sersic.ipynb` (+ figs `results/figures/cored_test_{A_boxsweep,B_profile}.png`, `cog_total_light.png`).
- NEW results: `mstar_headline_quiescent.npz`, `aperture_correction_validated.npz`, `sed_properties_final.npz`,
  `cog_total_light.npz`. Updated: `PAPER_VALUES.json`, `key_results_table.tex`, ApJL figure PDFs + doc PDFs.

## Outstanding TODOs
- **#10 PSF model** (webbpsf JWST / TinyTim or grizli HST / empirical) — prerequisite for any inner-profile/core-Sérsic
  work; ISMgas lacks webbpsf. **#9 core-Sérsic** — likely inconclusive (core unresolved); needs #10.
- **External:** M_BH PROVISIONAL until Ferrami resolves free-vs-fixed (fixed ~2e9 favored by his prose, UNRESOLVED).
  `DLOGM_COSMO`: regen `aperture_matched_photometry.npz` natively under Planck (drop +0.0282). DRAFTING §2.4.2 pipeline
  walkthrough still pre-audit (only §3.2 rewritten). METHODS `″` monospace-glyph warning (minor, pre-existing).
- **Possible next:** fixed physical-radius M⋆ (e.g. 30–50 kpc) as a clean alternative definition; better field-masking
  to extend the box/CoG test on the JWST mosaics (F150W2/F322W2 rail at large box now).

## How to resume / verify
```bash
conda activate ISMgas; cd ~/Documents/AGEL/AGEL0206_kinematics_and_photometry
python scripts/paper_values.py                          # regen registry
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md TESTS_AND_DIAGNOSTICS.md \
    METHODS_AND_SYSTEMATICS.md PROJECT_BRIEF.md          # GREEN
python scripts/make_latex_tables.py                     # key_results_table.tex
python -m scripts.cog_total_light                        # CoG cross-check (re-fits Bagpipes; ~min)
PATH="/opt/anaconda3/bin:/Library/TeX/texbin:$PATH" ./build_docs_pdf.sh DRAFTING_FACTS_paper_2026-05-29.md TESTS_AND_DIAGNOSTICS.md
# figures: AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb (agel0206_figs; Fig1 aplpy needs lens env)
# cored nb: notebooks/18_cored_vs_coreless_sersic.ipynb (run from notebooks/, it chdir's to repo root)
```
First action next session: **#10 PSF setup** (unblocks core-Sérsic / inner-profile), or the fixed-physical-radius M⋆
alternative. M⋆=11.46 and M_BH=5.2+1.7/−1.4e8 are settled for now; both stay registry-driven.
