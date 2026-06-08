# HANDOFF — σ_e systematics closeout + deterministic values registry — 2026-06-08

Covers the Tasks 1–8 arc requested 2026-06-08 ("go for tasks 1–7 … make new notebooks … then
build scripts that populate values deterministically, NOT LLM-typed"). New notebooks 13–16,
new drivers, and the `paper_values.py` registry.

## Headline change

**σ_e(<R_e) = 269.62 ± 13.27 km/s** (asym −13.45 / +13.10) — was ±11.77. The **D7 R_e-source
systematic (±6.13) was folded into the budget** (user decision: full 4-estimator spread). Central
value 269.62 unchanged. M★, R_e unchanged. **All headline numbers now emitted by
`scripts/paper_values.py` → `results/PAPER_VALUES.json` (single source of truth).**

## What each task produced

| Task | Result | Artifacts |
|------|--------|-----------|
| 1 σ_e masking systematic | **±5.85 km/s** (expert 269.62 / sersic 262.9 / perband 268.6 / global 274.6). Overlaps F200-mask ±6.65 → larger-of-two, **NOT added** (no budget change). | `scripts/run_sigma_e_mask_systematic.py`, nb13, `results/sigma_e_mask_systematic_N500.npz` |
| 2 figures copy | `figures_Mstar_principled_2026-06-08.ipynb` (figures repo), **registry-driven** (Fig 2/3/4 load PAPER_VALUES.json), Fig-1 errors cleared, stale σ_e renders cleared for `lens`-env re-run | (ApJL figures repo) |
| 4 CoG reconciliation | algos agree <0.04″ at matched r_max; gap was r_max (4″ vs 6″). `measure_Re.hst_Re` now `r_max_arcsec=6.0`. **Flag: raw CoG R_e non-convergent → R_e method sys ±0.08″** (headline 2.305″ = top of 2.1–2.5″ family) | nb14, `results/Re_cog_reconciliation.npz`, A3c |
| 5 D7 R_e-source @ wide | **±6.13 km/s** (mean 2.305"=269.62 / F140W=267.44 / F200LP=272.44 / CaHK+G 2.90"=279.69). Folded into budget (M12). Light-family-only ±2.50 | `scripts/run_sigma_e_Re_systematic_wide.py`, nb15, `results/sigma_e_Re_systematic_wide_N500.npz` |
| 6 Hδ decision | **KEEP UNMASKED.** Local-MAD: Hδ well-fit (0.44 < 0.81, not an outlier); masking shifts σ_e +6–8 km/s = LOSVD *information*, not contamination (M9 pattern). TODO in `bootstrap_ppxf.py` closed | `scripts/run_sigma_e_hdelta_test.py`, nb16, `results/sigma_e_hdelta_test_N500.npz` |
| 3 merge | branch `photometry-masking-mstar-2026-05-29` → main (see git) | — |
| 7 reduction-pass | **passive — no action.** ±3.45 rests on 2 reductions (NEW vs OLD `_mtwdo_`). Refine only if a 3rd reduction lands. No 3rd reduction exists; nothing to compute | (documented here + CLAUDE.md budget row) |
| 8 values registry | `scripts/paper_values.py` + `results/PAPER_VALUES.json` + drift linter | see below |

## Task 8 — deterministic values registry (design + status)

**Problem the user named:** headline numbers were hard-coded in the figures notebooks and repeated
across CLAUDE.md / TESTS / DRAFTING_FACTS / memory, re-typed by the LLM from local context — drift-prone.

**Architecture (single source of truth → consumers):**

1. **Registry — `scripts/paper_values.py` → `results/PAPER_VALUES.json`.** Reads ALL result caches
   and computes every headline number *from the data* with provenance + formula per entry:
   - σ_e central + asym stat (pools the 3-SPS `new_clean_hei` bootstrap caches)
   - σ_e budget: 7 named components → quadrature sys → asym totals (live: stat, R_e-source;
     carried-constant w/ provenance: I-shape, F200-mask, frame, centering, fit-window, reduction)
   - log M★ (10%/20%, fill reach, masking sys), R_e (arcsec/kpc, method sys), Hδ decision, constants
   - Regenerate: `python scripts/paper_values.py` (prints a table, writes the JSON).
2. **Consumer A — figures (DONE).** `figures_Mstar_principled_2026-06-08.ipynb` cells load
   `PAPER_VALUES.json` (`SIGMA_E_WIDE_SYS_ERR`, `sigma_e_p50/errup/errlo`, `Mstar_p50/errup/errlo`)
   instead of literals. **No LLM-typed values in the figure code.**
3. **Consumer B — docs (linter, DONE; auto-fill, PROPOSED).** `paper_values.py --check <files>` is an
   **anchored, precision-aware** drift linter: it parses only the canonical headline statements
   (`σ_e(<R_e) = N ± M`, `−lo/+hi`, `log(M⋆) = N`) and flags any that disagree with the registry.
   Skips cross-check/historical lines (marker heuristic + `<!-- pv-skip -->` escape hatch). Currently
   **GREEN** on CLAUDE.md + TESTS + DRAFTING_FACTS. Run it after any number change.
   - **Next step (proposed, not built):** sentinel-region auto-fill — wrap each doc's headline block
     in `<!-- PV:auto --> … <!-- /PV:auto -->` and add `paper_values.py --render` to overwrite the
     region's content from the registry. That makes the *block* generated, not just validated. (The
     fuzzy-proximity linter first attempted was abandoned — 155 false positives on legitimate nearby
     sweep numbers; the anchored+skip design replaced it.)

**Why not fully auto-generate the docs now:** the docs interleave many legitimate cross-check and
historical numbers with the live headline; a safe auto-fill needs the sentinel regions added first
(mechanical, low-risk follow-up). The linter already prevents silent drift in the meantime.

## New notebooks (all executed headless, ISMgas)

`13_sigma_e_masking_systematic` · `14_CoG_reconciliation` · `15_Re_source_systematic_wide` ·
`16_Hdelta_masking_decision`. Each loads its cache + figure under `results/figures/nb1{3,4,5,6}_*.png`.

## Open / flagged (not blocking)

- **R_e method systematic ±0.08″** (A3c) not folded into σ_e/M★ — would shift the aperture; flagged
  for a decision (a 2-component-Sérsic / PSF-deconvolved R_e would be the robust replacement).
- Doc sentinel-region auto-fill (Task 8 next step, above).
- Manuscript prose (G1 lens model / G2 X-ray / G3 peculiar velocities) — unchanged from prior handoff.
- Figures repo: re-run `figures_Mstar_principled_2026-06-08.ipynb` cells 8→end in the `lens` env to
  regenerate the σ_e-dependent Fig 2/Fig 3 renders (cleared here; Fig 1 needs aplpy interactively).

## Reproduce everything

```
conda activate ISMgas
python scripts/run_sigma_e_mask_systematic.py   --n_bootstrap 500   # Task 1
python scripts/run_sigma_e_Re_systematic_wide.py --n_bootstrap 500   # Task 5
python scripts/run_sigma_e_hdelta_test.py        --n_bootstrap 500   # Task 6
python scripts/paper_values.py                                       # registry → JSON
python scripts/paper_values.py --check CLAUDE.md TESTS_AND_DIAGNOSTICS.md DRAFTING_FACTS_paper_2026-05-29.md
```
