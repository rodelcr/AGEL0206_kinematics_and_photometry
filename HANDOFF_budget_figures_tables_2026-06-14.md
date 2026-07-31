# HANDOFF — error-budget closeout + figure/table regen — 2026-06-14

Resume point after a planned computer restart. Everything below is on disk and
self-consistent (registry ↔ docs ↔ figures ↔ tables); `paper_values.py --check` is GREEN.

## TL;DR — current headline (single source of truth: `results/PAPER_VALUES.json`)

- **σ_e(<R_e) = 267.31 −11.98/+11.58 km/s** (sym **±11.77**; stat −5.04/+3.98 ⊕ sys ±10.87), at the
  best-mask R_e = 2.097″ aperture.
- **log(M⋆/M☉) = 11.50 +0.09/−0.14** (aperture-corrected total, matched 2 R_e, 10% floor) ± 0.086 (masking sys).
- **R_e = 2.097″ = 15.26 kpc** (Planck 2015; method sys ±0.10″).
- **M_BH = 5.2 −1.3/+1.5 × 10⁸ M☉** (log 8.72; θ_E≈1.36″, γ≈1.31±0.08) — **PROVISIONAL** (Ferrami draft).
- z_l = 0.67564, z_s = 1.30263. Cosmology: Planck 2015 (H₀=67.7, Ω_m=0.302), 7.2764 kpc/arcsec at z_l.

## What changed this arc (2026-06-14)

1. **M⋆ error budget consolidated + de-staled** in DRAFTING §3.2: added a two-framing PV:auto block
   (`mstar_budget`) — headline-total (posterior +0.09/−0.14 ⊕ masking ±0.086) and the Sérsic-only
   named-component ±0.19 decomposition (mask ±0.125 dom). Fixed §3.2's pre-Planck literals (11.47→11.50,
   estimators, sub-components). Wired the M⋆ budget into `scripts/paper_values.py` (new `logMstar_budget`,
   `render_mstar_budget`).
2. **Dropped the σ_e "frame (vac/air) ±5" systematic** (user decision). Rationale (documented in
   CLAUDE.md + TESTS): it's a DETERMINISTIC per-SPS native-frame correction with <0.5 km/s σ impact
   (TESTS D4), and the ±5 was the inter-SPS spread (TESTS C3) ALREADY in the pooled 3-SPS bootstrap stat
   → double-count. Relabeled "F200 mask"→"arc masking" (the general masking term; subsumes the ±4.58
   mask-approach via larger-of-two). Budget: sys ±11.97→±10.87, total ±12.79→**±11.77**. Central 267.31
   unchanged. Synced registry + DRAFTING/CLAUDE/TESTS/METHODS headline + rebuilt DRAFTING/TESTS PDFs.
3. **Figure 2 fixed** (`AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb` cell[8]): σ_e data path
   `new_clean_hei`(269.6→"270") → `resys_best_mean`(**267**); M⋆ inset switched from superseded
   `perband_raw`(11.39) to the headline `total` estimator (**11.50**, dropped the fill-band overlay);
   broadened the M⋆ x-range to med±~0.7. Regenerated all 4 ApJL PDFs (sigma_e_SED_final_wide, Mbh_sigma,
   Mbh_Mstar, Mbh_Mstar_clean).
4. **Two LaTeX tables** (deterministic, read from PAPER_VALUES.json):
   - `results/key_results_table.tex` ← NEW `scripts/make_latex_tables.py` (σ_e, R_e, fluxes, M⋆, M_BH,
     M_BH/M⋆, θ_E, γ, system params).
   - `results/sersic_parameter_table.tex` ← `scripts/sersic_parameter_table.py` (appendix). **Fixed its
     cosmology bug**: KPC_PER_ARCSEC was hardcoded 7.04 (old H₀=70) → now reads 7.2764 from the registry;
     kpc column rescaled (e.g. F200LP 14.09→14.56).
   - Extended the registry with `photometry` (matched-2R_e AB mags), `lens_model` (M_BH/θ_E/γ, PROVISIONAL),
     coords; set z_source 1.302→1.30263.
5. **Cosmology reference notebook** built+executed: `cosmology_reference_AGEL0206.ipynb`
   (+ `results/cosmology_reference_AGEL0206.csv`). Lookback: deflector 6.368 Gyr, source 9.030 Gyr.
6. **Consistency cleanup (was TODO A — DONE):** de-staled `METHODS_AND_SYSTEMATICS.md` headline banner +
   headline table (269.62→267.31, R_e 2.305→2.097/15.26, M⋆ 11.16/11.33→11.50, asym −11.98/+11.58),
   `PROJECT_BRIEF.md` headline + table + poster mockup, and removed the live "Frame (±5.0)" bullet from
   DRAFTING §2.4.3. Linter GREEN on all 5 docs (now incl. PROJECT_BRIEF). **NOTE:** METHODS/PROJECT_BRIEF
   *narrative body prose* may still cite older historical σ values — those are flagged historical, not headline.

## Confirmed methodology fact (asked this session)
There is **no bootstrap in the M⋆/SED estimation** — the M⋆ posterior is the **Bagpipes** Bayesian
nested-sampling result (`run_bagpipes_for_mags`, `sampler=multinest/nautilus`, `n_live=400`, 10% flux
floor). **N=500 is exclusively the ppxf σ_e wild bootstrap** (`bootstrap_ppxf.py`). The Sérsic-param
errors use a separate seeded parametric bootstrap (n_boot=80).

## Files modified / created this arc
- `scripts/paper_values.py` — +logMstar_budget, render_mstar_budget, dropped frame component,
  +photometry/lens_model/coords blocks, z_s=1.30263.
- `scripts/sersic_parameter_table.py` — KPC now read from registry (was 7.04).
- `scripts/make_latex_tables.py` — NEW (key-results table generator).
- `DRAFTING_FACTS_paper_2026-05-29.{md,tex,pdf}`, `CLAUDE.md`, `TESTS_AND_DIAGNOSTICS.{md,pdf}`,
  `METHODS_AND_SYSTEMATICS.md` — headline/budget synced (see TODO: METHODS/PROJECT_BRIEF still partly stale).
- `results/PAPER_VALUES.json`, `key_results_table.tex`, `sersic_parameter_table.{tex,md,npz}`,
  `cosmology_reference_AGEL0206.csv`.
- `cosmology_reference_AGEL0206.ipynb` — NEW.
- `AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb` cell[8] (source); 4 regenerated PDFs.

## Outstanding TODOs (carry into next session)

**A. Consistency cleanup — ✅ DONE this arc** (see "What changed" #6). Linter GREEN on all 5 docs.

**B. Optional polish:**
- Re-run+SAVE `figures_paper4_2026-06-08.ipynb` so its INLINE outputs match the regenerated PDFs (needs
  the `lens` env for the Fig 1 aplpy cells; non-aplpy cells run in `agel0206_figs`). Clean the stale
  commented σ_e block at the top of cell[8].
- Legacy `figures.ipynb` (RE_TAG='2p305', M⋆=11.47) — fix or mark superseded (paper PDFs come from figures_paper4).
- Add M•–σ offset + mass-weighted-age rows to the key-results table; test-compile both .tex in AASTeX.

**C. External / decisions (not ours to close):**
- Lens-model numbers PROVISIONAL until Ferrami's final MultiNest runs (M_BH/θ_E/γ/15-model table).
- GAP G2: no X-ray data for quiescence. Source-z error bar [TBD]. Wild-bootstrap citation TODO.
  R_e method-sys (±0.10″) fold decision. `DLOGM_COSMO`: regenerate `aperture_matched_photometry.npz`
  natively under Planck 2015 (then drop the +0.0282 analytic shift).

## How to resume / verify after restart

```bash
conda activate ISMgas              # σ_e/M⋆/registry/tables/cosmology nb
# (figures use the agel0206_figs kernel; Fig 1 aplpy needs the `lens` env)

cd ~/Documents/AGEL/AGEL0206_kinematics_and_photometry
python scripts/paper_values.py                       # regenerate registry, print headline table
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md \
    TESTS_AND_DIAGNOSTICS.md METHODS_AND_SYSTEMATICS.md   # drift linter → expect GREEN
python scripts/make_latex_tables.py                  # regenerate results/key_results_table.tex
python scripts/sersic_parameter_table.py             # regenerate appendix table (seeded; ~30 s)
./build_docs_pdf.sh DRAFTING_FACTS_paper_2026-05-29.md TESTS_AND_DIAGNOSTICS.md
```

First action next session: **TODO B** — re-run+save `figures_paper4` so its inline outputs match the
regenerated PDFs (and clean cell[8]'s stale commented σ_e block). Then C/D as you prioritize.
