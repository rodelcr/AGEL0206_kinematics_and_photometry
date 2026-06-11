# HANDOFF — "Best mask throughout" + photutils CoG validation (2026-06-11, M13)

**TL;DR.** Per the user directive *"I think we need to use the best possible mask throughout
our analysis. period."*, the best-mask (single-Sérsic + color/morph gate + WCS-registered)
curve-of-growth **R_e = 2.097″** (was the expert-mask 2.305″) is now the headline R_e everywhere
— cascaded to R_e, σ_e, and the M★ aperture. Also validated our R_e/CoG machinery against
photutils. Headline moves are small and inside the systematics:

| quantity | before (M12) | after (M13) |
|---|---|---|
| **R_e** | 2.305″ = 16.23 kpc | **2.097″ = 14.76 kpc** (best mask; method sys ±0.100″) |
| **σ_e(<R_e)** | 269.62 ± 13.27 | **267.31 ± 12.79** (−12.98/+12.61), at R_e=2.097″ |
| **R_e-source sys** | ±6.13 (4-est, incl CaHK) | **±5.13** (best-mask CoG light family) |
| **log M★** | 11.47 +0.09/−0.15 | **11.47 +0.10/−0.13** (unchanged; R_e-robust) |

## What changed & why

1. **R_e → 2.097″ (best mask).** `re_mask_sensitivity.py`: the best masks (color/morph gate +
   WCS reg + companion) remove outer diffuse-arc/interloper light that inflated the expert-mask CoG.
   The HST 2-R_e companion mask is **empty** (companions live only in deep JWST) → the global
   color/morph mask alone gives **2.097″** (self-consistent fixed point). NB: an earlier value 2.078″
   was an artifact of a stale RE=2.305 companion mask; the fixed point is 2.097″.
   - Method systematic re-derived on best masks (`Re_bestmask_reconciliation.py`,
     `results/Re_cog_reconciliation_bestmask.npz`): raw CoG 2.097 / sky-sub 1.922 / Sérsic r_eff 1.897
     → **±0.100″**.

2. **σ_e re-measured at R_e=2.097″.** `run_sigma_e_Re_grid.py` (supersedes
   `run_sigma_e_Re_systematic_wide.py`): a 7-point σ_e-vs-R_e grid bracketing the new headline,
   reusing the 4 existing D7 caches + 3 new (resys_best_{F140W,mean,F200LP}). Headline σ_e =
   σ_e(best_mean) = **267.31 −5.04/+3.98 (stat)**; −2.3 km/s vs 2.305″ along the rising σ(R) profile.
   - **R_e-source systematic = ±5.13** — user decision: best-mask CoG **light family** only
     {1.912/2.097/2.281″}. CaHK+G (2.90″→279.7, different I-map definition) and the full grid (±9.98)
     are cross-checks, NOT folded. (Resolves the prior referee soft-spot of folding CaHK.)

3. **M★ re-measured at matched 2 R_e = 4.19″.** `aperture_2re_companions.py` RE 2.305→2.097;
   `aperture_matched_photometry.py` re-run (stale Bagpipes posteriors cleared first). **Total
   (headline) = 11.47, unchanged** — converges 11.49/11.47/11.45 across 1/2/2.5 R_e (R_e-robust by
   construction). raw 11.17 / raw+apcorr 11.35 / filled 11.35 / Sérsic-total 11.41. 5% floor
   re-checked at the new aperture: total 5%=11.49 vs 10%=11.47 (Δ+0.02) → 10% kept.

4. **photutils CoG/R_e validation** (`validate_Re_photutils.py`, `results/validate_Re_photutils.npz`):
   - UNMASKED: our R_e matches photutils `CurveOfGrowth` (direct 2D aperture sum) to **±0.004″**.
   - MASKED: our azimuthal-fill CoG reproduced by photutils `RadialProfile` (integrated) to **±0.002″**
     (best-mask mean 2.097→2.077).
   - A naive masked aperture-SUM (photutils dropping masked pixels) biases R_e **+0.25–0.41″ high** →
     our azimuthal-FILL is the correct masked-CoG treatment.
   - centroid_2dg robust (vs centroid_com 0.24–0.26″); centroid_quadratic fails on extended light
     (1.3–3.5″ — expected, it's a point-source tool). Headline uses centroid_2dg.

## Aperture-corrected photometry validation (`validate_apcorr_established.py`)
petrofit/statmorph are absent from ISMgas, so the apcorr chain was validated against **photutils**
(exact sub-pixel aperture photometry) + **scipy incomplete-gamma** (the analytic Sérsic enclosed-light
law, Graham & Driver 2005). ALL PASS — no production bug, M★ unchanged:
- b_n (Ciotti&Bertin99) vs exact `gammaincinv(2n,½)`: **≤0.05%**.
- Sérsic total-flux formula (`sersic_total_flux_analytic`, G&D05 eq.4) vs numeric render-to-20 R_e: **≤0.03%**.
- Aperture correction: rendered-model enclosed light vs analytic `γ(2n,b_n·2^{1/n})/Γ(2n)`: **Δ≤0.0007**.
- Empirical `F_raw` (hard-edge `in_aperture`) vs photutils `EllipticalAperture(method='exact')`: **≤0.09%**.
- `sum(model_full)` finite-cutout truncation vs to-∞ total: **≤0.19% (≤0.002 mag)**.
→ `results/validate_apcorr_established.npz`. (Validation-script gotchas, for re-runs: fit_sersic2d fits a
SUBIMAGE → render with the (x1,y1) offset; compare enclosed light at the model's OWN r_eff; photutils
needs `mask=`, not NaN injection.)

## Files

New scripts:
- `scripts/run_sigma_e_Re_grid.py` → `results/sigma_e_Re_grid_N500.npz`
- `scripts/validate_Re_photutils.py` → `results/validate_Re_photutils.npz`
- `scripts/validate_apcorr_established.py` → `results/validate_apcorr_established.npz`
- `scripts/Re_bestmask_reconciliation.py` → `results/Re_cog_reconciliation_bestmask.npz`
- `scripts/prep_fig2_data_bestmask.py` → `results/run_wide_sigma_e/resys_best_mean/figure2_data.npz`

Modified:
- `scripts/aperture_2re_companions.py` (RE=2.097), re-ran → `aperture_2re_masks.npz`,
  `aperture_matched_photometry.npz`, `aperture_floor5_check.npz`.
- `scripts/paper_values.py` — σ_e central → `resys_best_mean`; R_e → bestmask recon; R_e-source →
  grid `sys_re_bestlight`; aperture → 2×2.097; `render_headline` rewritten for the 5-estimator M★
  schema (old 10/20% fill_reach keys were stale). PAPER_VALUES.json regenerated; **drift-check OK**.
- Docs: `DRAFTING_FACTS_paper_2026-05-29.md` (+PDF rebuilt; auto-blocks rendered), `CLAUDE.md`,
  `TESTS_AND_DIAGNOSTICS.md` (A3/A3d/A4/D7 updated + M13 row).
- Figures (`../AGEL_0206_ApJL_Figures/figures.ipynb`, backup `.bak.bestmask_2026-06-11`): Fig 2
  (`AGEL0206_sigma_e_SED_final_wide.pdf`) repointed to resys_best_mean + sys 11.96; M•–σ
  (`Mbh_sigma_relation.pdf`) point → 267.31 ±12.61/−12.98. Re-executed; both regenerated. Fig 4 (M★)
  unchanged. (Fig 1 aplpy cells error — pre-existing env gap, unrelated.)
- Memory: `project_2026-06-11_best_mask_throughout.md` (+ MEMORY.md index; M12 demoted).

## Open / flagged
- σ_e cache caveat: headline now at `resys_best_mean` (R_e=2.097″), NOT `new_clean_hei` (2.305″).
  `paper_values.sigma_e_pool()` reads resys_best_mean; new_clean_hei is the 2.305″ cross-check.
- Raw CoG still sits at the top of the R_e method family (no sky pedestal subtraction); now bounded
  (±0.100″) and photutils-validated, not a flagged unknown.
- Fig 1 (lens cutouts) needs `aplpy` (absent in ISMgas) — separate env; unchanged by this work.
