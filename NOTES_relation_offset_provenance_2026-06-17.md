# Provenance & field-standard check — M•–σ / M•–M⋆ offset significance (2026-06-17)

Backs `scripts/relation_offset_significance.py` → `results/relation_offset_significance.npz`
(TESTS O1). Records (1) where the Greene+2020 relation parameters come from and how they were
verified, and (2) that the "how-many-σ-from-the-relation" methodology matches field standard.

## 1. Relation parameters — verified against Greene, Strader & Ho (2020)

**Reference:** Greene, Strader & Ho 2020, ARA&A 58, 257, "Intermediate-Mass Black Holes"
(arXiv:1911.09678; DOI 10.1146/annurev-astro-032620-021835). Zotero key **S7XU8UGY**.

**Functional form** (verified verbatim in the paper main text, §8 Scaling Relations):
> "Assuming that log(M_BH/M⊙) = α + β log(σ*/160 km s⁻¹) + ε, where ε is the intrinsic scatter,
> we present our fits in Supplemental Table 5 and Figure 3."

and the M•–M⋆ caption (Supplemental Table 5): "α + β log(M*/M0) + ε for **M0 = 3 × 10¹⁰ M⊙**".

**Adopted coefficients = Greene+2020 Supplemental Table 5, "Early" (early-type) subsample** — the
correct sample for AGEL J0206's deflector (massive, passive, E/S0). Confirmed against BOTH the in-repo
machine table `../AGEL_Mbh_sigma/Greene20_Supple_5.csv` AND the supplemental PDF
`../AGEL_Mbh_sigma/aa58_greene_supmat.pdf` (Supplemental Table 5, full uncertainty columns):

| Fit | Sample | α | β | ε (intrinsic scatter) | pivot |
|-----|--------|---|---|-----------------------|-------|
| M•–σ*  | **Early** | **8.03 ± 0.06** | **4.24 ± 0.25** | **0.43 ± 0.04** | σ/160 km s⁻¹ |
| M•–M*  | **Early** | **7.89 ± 0.09** | **1.33 ± 0.12** | **0.65 ± 0.05** | M*/3×10¹⁰ M⊙ |

(For context, the other Supplemental-Table-5 rows: M•–σ All-no-limits 7.88/4.34/0.53, All-limits
7.87/4.55/0.55, Late-no-limits 7.40/2.54/0.50, Late-limits 7.44/3.61/0.58; M•–M* All-no-limits
7.56/1.39/0.79, etc. We use **Early** throughout, matching the E+S0 comparison sample plotted in the
paper figures and the deflector's morphology.)

These are **exactly** the values hard-coded in the paper figure notebook
(`AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb`, M•–σ / M•–M⋆ cells) and in
`scripts/relation_offset_significance.py` — so the significance calc, the figures, and Greene+2020
are mutually consistent. ✅

**Verification trail (2026-06-17):** Zotero fulltext of S7XU8UGY confirmed the M•–σ functional form +
σ/160 pivot; web cross-check confirmed the M•–M⋆ Early row (7.89/1.33, pivot 3×10¹⁰) is the version
used across the recent literature; the in-repo `Greene20_Supple_5.csv` and `aa58_greene_supmat.pdf`
confirmed all central values **and** uncertainties to the quoted precision. No fabricated numbers.

## 2. Methodology — "how many σ from the relation" is field standard (VERIFIED against other analyses)

For a single object the offset from a scaling relation is judged against the **total** dispersion =
measurement uncertainty combined in quadrature with the relation's **intrinsic scatter** (the
population dispersion about the mean). An object is "consistent with the relation" if
|Δ| ≲ √(σ_meas² + ε_int²). Verified (2026-06-17, full text read in Zotero) against the closest analyses:

- **Pacucci, Nguyen, Carniani, Maiolino & Fan 2023** (ApJL 957, L3) + **Pacucci & Loeb 2024**
  (arXiv:2401.04159, read in full) — the reference high-z offset analysis. They quantify the offset of
  a point from the **local relation** and assess significance **against the intrinsic scatter combined
  with the measurement errors**: their ">3σ violation" of the M•–M⋆ relation holds "*unless… errors of a
  factor ∼60 in their black hole mass or their stellar mass all in the same direction*" — i.e. exactly
  offset ÷ (measurement ⊕ intrinsic scatter), the same construction we use. They fit the relation with a
  free **intrinsic-scatter** term (Hogg, Bovy & Lang 2010 likelihood; orthogonal scatter ν) and correct
  for Malmquist-type bias (Lauer+2007) — the population-fit generalization of our single-point test.
- **Melo-Carneiro et al. 2025** (Cosmic Horseshoe, MNRAS; the DIRECT analogue — a lensing+dynamics SMBH
  placed on M•–σ at z_l=0.44, read in full) — 5σ BH detection; measurement uncertainty from **200 Monte
  Carlo realizations (standard deviation as the error)**; the BH is reported as an over-massive M•–σ
  outlier with consistency quoted "within Nσ" against the relation, same σ-framework.
- **Reines & Volonteri 2015**, **Kormendy & Ho 2013** — report residuals about the relation relative to
  the (vertical) intrinsic scatter ε.

**Convention check / nuances (honest):**
- **Vertical vs orthogonal scatter.** Pacucci fits the *relation itself* and uses *orthogonal* intrinsic
  scatter (because slope+scatter+points are inferred jointly). We instead place ONE point on a **fixed,
  published** relation, so the correct convention is the **vertical** offset ΔlogM• at fixed x compared
  to the relation's **vertical** intrinsic scatter — which is exactly how Greene+2020 define ε (their
  form is logM• = α + β logx + ε, ε in the log M• direction). So vertical-offset ÷ vertical-ε is
  internally consistent; mixing in Pacucci's orthogonal ε would be the inconsistency, and we avoid it.
- **x-error propagation.** Measurement error in the independent variable (σ_e, M⋆) is propagated
  vertically through the slope, β·σ_x, then added in quadrature — standard.
- **Relation-parameter (α,β) errors.** Optional; for a single object the measurement ⊕ intrinsic-scatter
  terms dominate. We report N_σ both without and with them (0.57/0.69 → 0.56/0.67) — negligible.
- **Intrinsic scatter must be in the denominator for a single object.** This is the crux and it is
  standard: ε is the population dispersion, so a lone galaxy is "on the relation" if it lies within
  √(meas² + ε²). Reporting Δ/ε alone (M•–σ 0.60, M•–M⋆ 0.75) gives essentially the same answer here
  because ε dominates.

**Conclusion: our treatment is consistent with the field standard.** We use the vertical single-point
convention matched to Greene+2020's published vertical ε (the appropriate choice when placing a point on
a fixed relation), with measurement and x-errors propagated in quadrature — the same offset-vs-combined-
scatter logic Pacucci+2023/2024 and Melo-Carneiro+2025 apply.

**Our estimator** (`relation_offset_significance.py`):
```
N_σ = |Δ| / sqrt( σ_M•²  +  (β·σ_x)²  +  ε_int²  [+ σ_α² + (σ_β·log(x/pivot))²] )
Δ   = log M•(obs) − [α + β·log(x_obs/pivot)]          (offset from the MEAN relation)
```
- `σ_M•`     : our log M_BH measurement error (PROVISIONAL; Ferrami free-BH 1σ envelope)
- `β·σ_x`    : our σ_e / log M⋆ error propagated through the slope
- `ε_int`    : Greene+2020 intrinsic scatter (the dominant term)
- optional   : the relation-parameter (α, β) uncertainties — the "full" variant

This is the **pythagorean** combination the analysis asked for. We report three numbers per relation:
(i) the dex offset from the scatter-free mean line; (ii) the distance to that mean line in measurement-σ
only; (iii) the consistency N_σ including intrinsic scatter (the field-standard single-object test).

## 3. Results (with M_BH = 5.2 ₊₁.₇/₋₁.₄ ×10⁸, PROVISIONAL; σ_e=267, logM⋆=11.46)

| relation | Δ (mean) | meas-only | **N_σ (incl. ε_int)** | N_σ (+α,β) | Δ/ε_int |
|----------|----------|-----------|------------------------|-----------|---------|
| **M•–σ**  | −0.26 dex | 1.75σ | **0.57σ** | 0.56σ | 0.60 |
| **M•–M⋆** | −0.49 dex | 1.79σ | **0.69σ** | 0.67σ | 0.75 |

(`Δ/ε_int` = offset in pure intrinsic-scatter units, the simplest field metric — agrees with the full
quadrature because ε_int dominates.) **Verdict:** AGEL J0206 is a mildly under-massive BH relative to
the *mean* early-type relations but lands **well within their intrinsic scatter (<0.7σ)** — i.e.
consistent with both, as a genuine independent high-z (z=0.676) probe. The most distant calibrator in
the Greene+2020 Early sample is A1836-BCG at 152 Mpc (z≈0.034), so AGEL J0206 extends the test ~20× in
redshift.

## 4. Known inconsistency to reconcile (flagged, not yet changed)

The registry leaf `lens_model.mbh_sigma_offset` (−0.69) was computed against [mean + ε_int] — i.e. the
intrinsic scatter was **added to the mean relation** before differencing. That is not the offset from
the relation (the mean-line offset is **−0.26**), and it is not where the relation is drawn in the
figure. Per user (2026-06-17) the figure scatter band is intentionally the mean-locus 1σ (how we choose
to treat per-object scatter when plotting) and is kept; the registry offset reconciliation
(−0.69 → −0.26, add `mbh_mstar_offset` = −0.49, add `n_sigma` leaves) is deferred pending authorization.

## 5. Redshift reach of the Greene+2020 calibrating sample (the independent-probe point)

The Greene+2020 relations are calibrated on **local dynamical** BH masses. The supplemental tables list
**distances** (`Distance_Mpc`), not redshifts; converting with z ≈ H₀·D/c (H₀ = 67.7, Planck 2015) gives
the redshift reach. The most distant calibrator in **both** relation samples is the same object:

| relation | sample (as built in the figures) | N | most distant galaxy | distance | z ≈ H₀D/c |
|----------|----------------------------------|---|---------------------|----------|-----------|
| **M•–σ**  | Suppl-2 (E+S0) + Suppl-3 (E) + Suppl-4 (S0) | 83 | **A1836-BCG** | 152.4 Mpc | **≈ 0.034** |
| **M•–M⋆** | Greene Fig-10 E+S0 with logM⋆ | 61 | **A1836-BCG** | 152.4 Mpc | **≈ 0.034** |

Next most distant (both samples), all giant ellipticals/BCGs: NGC 6086 (138 Mpc, z≈0.031), NGC 7768
(116 Mpc, z≈0.026), NGC 4889 (102 Mpc, z≈0.023), NGC 3842 (92 Mpc, z≈0.021); the overwhelming majority
sit far closer (Virgo/Coma and nearer). The full 120-galaxy Greene+2020 dynamical compilation also caps
at A1836-BCG.

**Source:** `../AGEL_Mbh_sigma/Greene20_Supple_{2,3,4}.csv` (`Distance_Mpc` column), samples assembled
exactly as in the paper figure code (M•–σ = Suppl-2 HT∈{E,S0} ∪ Suppl-3 ∪ Suppl-4 HT=S0; M•–M⋆ =
Suppl-3 HT=E ∪ Suppl-4 HT=S0 with finite logM⋆). Reproduce: `python -m scripts.greene_sample_redshift`. **Caveat:** these are redshift-independent distances
(SBF/Cepheid for the nearest, Hubble-flow for the farthest); the z conversion is exact only in the
Hubble-flow limit — reliable for A1836-BCG (~3440 km/s) but a recession-velocity proxy for the nearest
members where peculiar velocities dominate.

**Why it matters:** the local calibration extends only to **z ≈ 0.034**. AGEL J0206 at **z = 0.67564**
is ~20× higher in redshift (and ~4.4× in comoving distance) than the most distant calibrator — the core
justification for treating this as a genuinely independent, high-z probe of the relations.
