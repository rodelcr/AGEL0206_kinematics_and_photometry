# Literature review — M_BH–σ slope & its evolution (for the Discussion)

**Date:** 2026-06-16
**For:** AGEL J0206 ApJL discussion section.
**Target object:** passive **elliptical** deflector, **σ_e ≈ 267 km/s** (log σ_e ≈ 2.43),
**z = 0.676** (lookback ≈ 6 Gyr), placed on the local M_BH–σ relation.

> **Provenance discipline:** every number below was read from the **fulltext PDF in
> Zotero** (section / table / equation cited), not from web summaries. The earlier
> deep-research web run had its verification phase rate-limited (all votes abstained),
> so its numbers were *unverified* — this doc supersedes it. Two papers are **abstract-only
> in the library** and are flagged ⚠ABSTRACT-ONLY; do **not** quote their slope/zero-point
> until the PDF is attached and the value re-checked.

All papers are tagged `Msigma-litreview-2026-06` in Zotero. Zotero item keys in `[...]`.

---

## Theme A — Local calibration & intrinsic scatter (the anchor AGEL0206 sits on)

The local M_BH–σ slope β is **not a single agreed number** — it ranges β ≈ 4.0–5.6
depending on sample composition and which galaxy types are included.

| Reference | Relation (form) | β (slope) | ε intrinsic (dex) | Notes |
|---|---|---|---|---|
| **Kormendy & Ho 2013** `[M27KJIWH]` | log(M•/10⁹)=−0.509±0.049 + (4.384±0.287)log(σ/200) [Eq. 3] | **4.38±0.29** | **0.29±0.02** | classical bulges + ellipticals ONLY; pseudobulges/bulgeless excluded. **High-end saturation:** "M• becomes essentially independent of σe … in ellipticals that have cores" (§6.7); NGC 3842/4889 excluded. |
| **McConnell & Ma 2013** `[QUJKBDGH]` | log M• = 8.32±0.05 + (5.64±0.31)log(σ/200) [Eq. 2] | **5.64±0.31** (full) | **0.38** | early-type 8.39+5.20logσ vs late-type 8.07+5.06logσ — similar slope, offset intercept (early ~2× M• at fixed σ). Core vs power-law ETGs offset ~2×. Steep 5.64 comes from MIXING types. MPFITEXY returns ~0 fitted scatter for σ>275 km/s (artifact of steep slope). 72 galaxies. |
| **Gültekin et al. 2009** `[UET4CWDU]` | log M• = α + β log(σ/200) [Eq. 2] | all **4.24±0.41**; E-only **3.96±0.42** | all **0.44±0.06**; E-only **0.31±0.06** | scatter larger than earlier work; **selection by R_infl/resolution biases to larger M• AND larger slope** (warns against sphere-of-influence cuts). 49 M• + 19 limits. |
| **van den Bosch 2016** `[LW5VV2S2]` | log M• = −4.00±0.51 + (5.35±0.23)log σ [Table 1, Fig. 1] | **5.35±0.23** | **0.49±0.03** | 230 galaxies, dwarfs→BCGs. BH tracks **global** σ, not bulge specifically; "universal," extends to bulgeless. Largest scatter (large sample). |
| **Saglia et al. 2016** `[6IXNFZJJ]` ⚠ABSTRACT-ONLY | M_BH–σ + BH fundamental plane | **NOT VERIFIED** | 0.26 (KH13 45-gal subsample, *bivariate FP*; 1-D M–σ not verified) | 97 galaxies (31 core E, 17 power-law E, 30 classical bulges, 19 pseudobulges); pseudobulges lower M• & decoupled. **PDF missing — slope/zero-point unverifiable.** |

**Take for the discussion:** AGEL0206 at σ_e ≈ 267 km/s sits **near the top of the local
σ baseline**, in the regime where (i) KH13 and McConnell & Ma report the relation
*saturating / returning anomalous fitted scatter* (≳275 km/s), and (ii) the slope itself
is most sample-dependent. So the local prediction for M• is sensitive to which calibration
is adopted (β ≈ 4.4 vs 5.6 → different M• at fixed σ).

---

## Theme B — Slope VARIES with galaxy type & mass (the high-σ steepening) — MOST RELEVANT

This is the strongest, most directly applicable theme: at the **high-σ / massive-elliptical
end where AGEL0206 lives, the M_BH–σ slope steepens well above the canonical ~4–5.**

- **Graham 2026, "Galaxy morphology dependent (M_bh)–σ₀ relations"** `[VX56BGZR]`
  (arXiv 2606.05808; MNRAS submitted; single author A.W. Graham) — **the key paper.**
  - Massive **E + ES,e** galaxies: **β = 7.83 ± 1.32** (Table 2 #1), at pivot **log σ₀ = 2.43**
    — i.e. the relation is pinned at exactly AGEL0206's σ_e.
  - Adding 16 ultramassive-BH BCGs steepens to **β = 8.62 ± 1.55** (#2).
  - S0s: β ≈ 4.4–4.6 (wet) down to 2.5–3.1 (primeval, dust-poor); spirals β ≈ 2.6.
  - Mechanism ("Virial Mirror," §4.1): "dry mergers drive galaxies almost vertically upward
    in the M_bh–σ₀ diagram, steepening the … slope of ∼4–5 into the observed ∼7–9 regime."
    SMBH is a "passenger," AGN feedback in "maintenance mode" at high mass.
  - **Prescription (§3.5): "if σ₀ > 230 km s⁻¹, use Equation 1"** (the steep β=7.83 relation)
    → applies to AGEL0206.
  - Explicit **caution against a monolithic relation** (abstract; §3.3.1 flags Greene+2020-style
    IMBH extrapolation of single power laws as "under-predicting the ultra-massive black holes").

- **Dullo, Gil de Paz & Knapen 2020, "Ultramassive black holes … M_BH–σ vs M_BH–R_b"** `[3VFSNNB9]`
  (arXiv 2012.04471).
  - Sérsic + normal-core: **β = 4.88 ± 0.29**; all incl. core-Sérsic: **β = 5.76 ± 0.37**;
    core-Sérsic only: **β = 10.67 ± 4.90** (steep, baseline-limited). [Table 3, BCES bisector]
  - **Large-core (ultramassive) galaxies offset 2.5–4σ ABOVE M_BH–σ**; Holm 15A ~1.1 dex above;
    M_BH–σ underpredicts the most massive BHs "up to a factor of 40."
  - **M_BH–R_b (break radius) is the tighter relation:** β = 1.20 ± 0.14, ~62% less scatter.
  - Cause: repeated **dry major mergers add M_BH while keeping σ ~unchanged** → vertical scatter;
    σ–L saturates at M_V ≲ −23.5. Conclusion 4: "the assumption that all … elliptical galaxies
    obey a single M_BH–σ relation may be invalid. Extremely massive galaxies … treated separately."

- **Sahu, Graham & Davis 2019, early-type M_BH–M\*** `[EKMRYCCG]` (ApJ 876, 155).
  - M_BH–M*,sph β = 1.27 ± 0.07; M_BH–M*,gal β = 1.65 ± 0.11.
  - **E (no-disk) vs ES/S0 (disk) are two offset relations: intercepts differ by 1.12 dex.**
  - No Sérsic/core-Sérsic bend WITHIN ETGs (historical "bend" reattributed to ETG-vs-LTG morphology).
  - (M_BH–σ deferred to a follow-up; no σ slope here.)

- **Hartmann et al. 2014, "The effect of bars on the M•–σ_e relation"** `[8BQZQ5IA]` (N-body sims).
  - Bars inflate central σ_e by ~12% (up to 20–40%) → barred galaxies shift **below** the relation
    (δα ≈ −0.19 to −0.20 dex) and inflate scatter.
  - **Contrast argument for AGEL0206:** "SMBHs in elliptical galaxies and in classical bulges follow
    a similar M•–σ_e relation … lack of an offset between ellipticals and unbarred classical bulges"
    — a passive elliptical is on the CLEAN, unoffset relation, free of bar-driven σ bias.

- **King 2026, "Ultramassive black holes and the three M–σ relations"** `[B2QKSH47]` ⚠ABSTRACT-ONLY
  (MNRAS 547, stag295). From abstract: **same slope M ∝ σ⁴ in spirals, field ellipticals, cluster
  ellipticals, but three different normalizations**; a mixed-type sample yields a fitted power slightly
  >4. Interpretation: longer accretion episodes needed to expel gas in ellipticals/cluster ellipticals.
  **Body unavailable — per-type normalization constants NOT verified.**

**Take:** the literature robustly documents slope *steepening at the massive-elliptical end*
(Graham 2026 β≈7.8 at the AGEL0206 pivot; Dullo+2020 ultramassive offsets ≤40×; cited Graham &
Scott 2013 core-Sérsic ~7). Whether this is genuine slope change or **parallel relations with
type-dependent normalization** (King 2026) is the open question AGEL0206 can speak to.

---

## Theme C — Redshift evolution of the scaling relations (AGEL0206 is at lookback ~6 Gyr)

- **Shen et al. 2015, "SDSS-RM: No Evidence for Evolution in M•–σ\* to z~1"** `[BZ4GU5GU]`.
  - **No intrinsic evolution of M•–σ* to z~1**: intercept at σ=200 stable (α ≈ 8.38, consistent
    with KH13) across 0.1<z<1 (46 objects at z>0.6). "a constant M•–σ* relation is favored to z∼1."
  - The *apparent* flattening (fitted slope 1.54→1.08) is **luminosity-threshold selection bias
    (Lauer+2007)**, not intrinsic. 88 broad-line quasars, single-epoch virial masses.
  - **Directly brackets AGEL0206's redshift** — and finds the M•–σ plane stable there.

- **Farrah et al. 2025, "Assembly of SMBHs at z<1 in ETGs"** `[ITQYBMLC]`.
  - **M•–M\***: z~0.8 ETG relation lies ~1 dex below local (z<0.2) ETGs in M• (net 0.4–1 dex; ~0.6 dex
    param uncertainty), ~0.5 dex above local AGN; **slope consistent** (1.25±~0.5 vs local 1.10).
  - Mg II virial masses, bias-corrected. Flags **M•–σ as "more fundamental and less prone to bias
    (Shankar et al. 2017)"** than M•–M*.

- **Pacucci & Loeb 2024, "Redshift Evolution of M•–M\* for JWST SMBHs at z>4"** `[MEFG9T7M]`
  (arXiv 2401.04159).
  - Model: **M•/M\* ∝ (1+z)^(5/2)** → z=5 ≈ 1.74 dex overmassive. Only the **normalization** evolves
    (refit slope 1.12±0.08). Notes **M•–σ and M•–M_dyn HOLD at 4<z<7** (Maiolino+2023) while M•–M* evolves.
  - Milder literature alternatives cited: Wyithe & Loeb 2003 (1+z)^1.5; Bennert+2011 (1+z)^1.15.

- **Bhattacharyya & Mangalam 2018, "Evolution of the M•–σ relation"** `[7CZKFVMU]` (theory, arXiv 1808.04536).
  - Predicts the **slope p of M•∝σ^p evolves as p ∝ (1+z)^(−α), α ≈ 0.24–0.34 up to z≈1** (p between 4–5);
    slope *decreases* toward higher z (modest at z=0.68).
  - Proposes **TMT** to probe M•–σ evolution at high z.

- **Reines & Volonteri 2015, local M•–M\*** `[EIVT6ZVS]` (the z=0 benchmark; no evolution measured).
  - Local AGN: log M• = 7.45±0.08 + (1.05±0.11)log(M*/10¹¹). Dynamical-BH ellipticals: 8.95±0.09 +
    (1.40±0.21), scatter 0.47 dex. **AGN normalization >1 dex below dynamical-BH ellipticals** — the
    offset that frames every high-z M•–M* comparison.

**Take:** at AGEL0206's redshift the *M•–σ normalization is observed to be stable* (Shen+15), while
*M•–M\* shows a ~dex-level offset* (Farrah+25) — consistent with M•–σ being the more fundamental,
less-evolving plane. A dynamical-side (lensing+kinematics) M•–σ point at z=0.68 is a clean test of this.

---

## Theme D — Future experiments to measure M• at distance / test slope & z-evolution

- **Greene, Strader & Ho 2020 (ARA&A 58, 257)** `[S7XU8UGY]` — **IMBH review** (M• ≈ 10²–10⁵ M⊙).
  Use only for the **low-mass / IMBH end** framing — NOT a high-mass M•–σ calibration.
  - §8.1: "we see no evidence for a change in M_BH–σ* relations at low σ*" and "we do not see evidence
    for flattening at 10⁵ M⊙" (once **upper limits** are included) — disfavors pure heavy-seed models.
  - **Selection-bias direction (corrected):** detection-only fits give a **shallow/flat** late-type
    M•–σ slope (bias toward most massive BHs at fixed galaxy property); including limits restores
    consistency. (The "slightly steeper" phrase in the review is for **M•–M\*** vs Reines & Volonteri.)
  - §8.3: "extreme caution in using AGN signatures to infer scaling relations … in this low-mass regime."

- **The GRAVITY+ Project** `[NXHXS7QP]` (arXiv 2301.08071).
  - **NIR interferometric spectroastrometry of the BLR → dynamical-style M• for "hundreds of AGN
    across cosmic time"** (a few 100 with GRAVITY-Wide + LGS AO). Reaches z>2 (Hβ in K to z~3).
    Independent of RM-calibrated scaling relations. Completion ~2026.

- **Golubchik et al. 2024, "Reverberation mapping of high-mass and high-z quasars using gravitational
  time delays"** `[F7RHI3SH]` (arXiv 2408.00073) — **most relevant to a strong-lensing paper.**
  - Use **strong-lensing time delays between multiple images of cluster-lensed quasars** to make RM
    feasible at high z / high L (where rest-frame lags are ~300–1000 d): the leading image's line
    response is seen ~one lag earlier; magnified multi-images give denser sampling.
  - Forecast yields: z>4 lensed broad-line AGN detectable — **dozens of thousands with JWST, hundreds
    with Euclid, thousands with Roman**, all-sky. Calibrates the R–L relation / BH-mass scaling in the
    early Universe.
  - Needs **galaxy clusters** (galaxy-scale lens TDs are too short). Proof of concept: Williams et al.
    2021, lensed quasar SDSS J2222+2745, C IV RM lag from image time delays.

**Take:** the natural forward-looking arc for an AGEL/strong-lens paper is (1) dynamical M• at
cosmological distance via interferometry (GRAVITY+) and (2) **lensing-time-delay RM of cluster-lensed
quasars (Golubchik+2024) — a strong-lensing route to high-z M•** that complements the dynamical
deflector-σ point this paper provides.

---

## Suggested narrative arc for the Discussion

1. **Place the point:** σ_e=267 km/s puts AGEL0206 at the high-σ end of the local baseline, where
   the M•–σ slope is most uncertain (β from 4.24 Gültekin → 4.38 KH13 → 5.35–5.64 McConnell&Ma /
   vdBosch) and where the relation may saturate (KH13 §6.7; McConnell&Ma σ>275 km/s).
2. **Slope-by-type framing:** cite Graham 2026 (β≈7.8 at log σ=2.43, the AGEL0206 pivot) and
   Dullo+2020 (ultramassive offsets ≤40×, M_BH–R_b tighter) for the massive-elliptical steepening;
   note the King 2026 alternative (parallel σ⁴ relations, different normalizations). Use Hartmann+2014
   to argue the *elliptical* deflector is on the clean, bar-unbiased relation.
3. **Redshift framing:** Shen+15 (M•–σ stable to z~1) and Farrah+25 (M•–M* offset, M•–σ "more
   fundamental") → a dynamical-side M•–σ point at z=0.68 tests normalization evolution where it's
   predicted small (Bhattacharyya & Mangalam 2018; Pacucci & Loeb 2024 for the contrasting high-z M•–M*).
4. **Future:** GRAVITY+ and lensing-time-delay RM (Golubchik+2024) as the path to populating M•–σ at
   cosmological distance — directly in the strong-lensing wheelhouse.

---

## Data-quality flags / TODO before submission
- ✅ **PDFs now attached for the whole batch** (2026-06-16) — Saglia 2016 `[6IXNFZJJ]` and King 2026
  `[B2QKSH47]` had only paywalled publisher links, so the free arXiv PDFs (1601.00974, 2602.13382)
  were merged onto the journal records. **Still TODO:** their slope/normalization numbers were read
  from the *abstract* only (the body wasn't available when the readers ran) — re-read the now-attached
  fulltext before quoting any Saglia/King value (marked ⚠ABSTRACT-ONLY in Themes A & B).
- **Re-run a retraction/correction check** (`scite_check_retractions tag=Msigma-litreview-2026-06`) —
  the API was unreachable on 2026-06-16.
- **Greene+2020 is the IMBH review**, scope-limited to the low-mass end — do not cite it for a
  high-σ M•–σ calibration.
- Numbers here are verified against fulltext; the earlier web deep-research numbers were NOT (rate-
  limited verification) — prefer this doc.

## Reference list (verified; Zotero keys + identifiers)
| Paper | Zotero key | DOI / arXiv |
|---|---|---|
| Kormendy & Ho 2013, ARA&A 51, 511 | M27KJIWH | 10.1146/annurev-astro-082708-101811 |
| McConnell & Ma 2013, ApJ 764, 184 | QUJKBDGH | 10.1088/0004-637X/764/2/184 |
| Gültekin et al. 2009, ApJ 698, 198 | UET4CWDU | 10.1088/0004-637X/698/1/198 |
| van den Bosch 2016, ApJ 831, 134 | LW5VV2S2 | 10.3847/0004-637X/831/2/134 |
| Saglia et al. 2016, ApJ 818, 47 ⚠ | 6IXNFZJJ | 10.3847/0004-637X/818/1/47 (arXiv 1601.00974) |
| Graham 2026 (MNRAS submitted) | VX56BGZR | arXiv 2606.05808 |
| Dullo, Gil de Paz & Knapen 2020, ApJ | 3VFSNNB9 | arXiv 2012.04471 |
| Sahu, Graham & Davis 2019, ApJ 876, 155 | EKMRYCCG | 10.3847/1538-4357/ab0f32 (arXiv 1903.04738) |
| Hartmann et al. 2014, MNRAS | 8BQZQ5IA | arXiv 1309.2634 |
| King 2026, MNRAS 547, stag295 ⚠ | B2QKSH47 | 10.1093/mnras/stag295 |
| Shen et al. 2015, ApJ 805, 96 | BZ4GU5GU | 10.1088/0004-637X/805/2/96 |
| Farrah et al. 2025, ApJ | ITQYBMLC | 10.3847/1538-4357/adb0c7 |
| Pacucci & Loeb 2024 | MEFG9T7M | arXiv 2401.04159 |
| Bhattacharyya & Mangalam 2018 | 7CZKFVMU | arXiv 1808.04536 |
| Reines & Volonteri 2015, ApJ 813, 82 | EIVT6ZVS | 10.1088/0004-637X/813/2/82 (arXiv 1508.06274) |
| Greene, Strader & Ho 2020, ARA&A 58, 257 (IMBH) | S7XU8UGY | 10.1146/annurev-astro-032620-021835 (arXiv 1911.09678) |
| GRAVITY+ Project 2022 | NXHXS7QP | arXiv 2301.08071 |
| Golubchik et al. 2024 | F7RHI3SH | arXiv 2408.00073 |
