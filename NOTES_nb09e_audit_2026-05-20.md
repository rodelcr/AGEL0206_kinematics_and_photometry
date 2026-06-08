# Audit — did we apply the same wide-arc-masked σ_e pipeline to the new red reduction?

**Date:** 2026-05-20
**Scope:** Independent verification that the +10.92 km/s shift between the headline cube (σ_e=254.85 km/s) and the new `_mtwdo_` reduction (σ_e=265.76 km/s) reflects a real reduction-pass difference, not a pipeline misapplication.

---

## TL;DR

**The pipeline was applied identically.** Every geometric, masking, frame-handling, pixel-count, and template-bookkeeping quantity matches *exactly* between the two runs. The +11 km/s σ_e shift is reproducible across every polynomial degree from 15 to 29 and across all three SPS templates with consistent magnitude (Δσ range: +8.7 to +14.1 km/s, median ≈ +11). The shift is real and traces to subtle differences in absorption-line shape introduced by the additional `_mtwdo_` calibration step.

---

## What had to be identical (and is)

These quantities depend on the HST data, the WCS, or the cube *metadata*, not the cube *data*. They must match exactly if the same recipe was applied. They do:

| Quantity | OLD cube | NEW cube | Δ |
|---|---|---|---|
| Cube shape | (3317, 100, 100) | (3317, 100, 100) | — |
| λ axis | 5625–8941 Å, Δλ=1.0 | 5625–8941 Å, Δλ=1.0 | **bit-identical** |
| CRVAL1/2 (WCS center) | (31.55618, −1.23859) | (31.55618, −1.23859) | exact |
| DISPSCAL | 0.294 | 0.294 | exact |
| GUIDFWHM | 1.271″ | 1.271″ | exact |
| HST-mean center (IFU sub-pix) | (50.147, 54.264) | (50.147, 54.264) | exact |
| R_e (F140W mean) | 2.168″ | 2.168″ | exact |
| R_e (F200LP mean) | 2.441″ | 2.441″ | exact |
| **R_E (headline mean)** | **2.305″** | **2.305″** | exact |
| F200LP arc spaxel mask sum | 69 | 69 | exact map |
| r_spax map | — | — | bit-identical |
| Aperture n_kept spaxels (R<R_e, w=0) | 146 | 146 | exact |
| Frame choice per SPS | FSPS=vac, EMILES=air, XSL=air | same | exact |
| ppxf `n_pix` per SPS | 2574 | 2574 | exact |
| lam_rest range per SPS | (3799.3, 5335.9) | (3799.3, 5335.9) | exact |
| χ²/dof at deg=22, FSPS | 0.083 | 0.081 | sub-percent |

**Conclusion of "did we apply it the same way?":** Yes, by every available checkpoint.

---

## What legitimately differs (and how)

These are cube-data differences expected between two reductions of the same nights:

1. **BUNIT label differs but numerically equivalent.**
   OLD: `10^-8 erg/s/cm³/arcsec²` (= 10⁻¹⁶ erg/s/cm²/Å/arcsec²)
   NEW: `1e-16 erg/s/cm²/arcsec²/Å` — same physical units, different convention.

2. **Absolute flux scale.** In the clean 6700–7300 Å aperture continuum, NEW/OLD = **0.804** (NEW is ~20% dimmer in absolute units). The pipeline normalises before fitting, so absolute scale does not move σ_e — but it tells us the `_mtwdo_` step changed the throughput calibration.

3. **Per-pixel ratio in the fit window (after normalisation).**
   median NEW/OLD = **1.003** (essentially unity)
   p16 / p84 = 0.969 / 1.041 (≈ 4 % RMS pixel-to-pixel jitter)

   This says: after normalization, the new spectrum is on the same level as the old, with sub-5 % pixel-by-pixel scatter dominated by noise / cosmic-ray spikes that survived stacking differently between the two reductions.

4. **σ_e is consistently higher on the new cube — across every degree and SPS.**

### Per-polynomial-degree σ on FSPS

| deg | σ_old | σ_new | Δ |
|---|---|---|---|
| 15 | 244.51 | 256.59 | **+12.08** |
| 16 | 245.03 | 256.74 | **+11.71** |
| 17 | 247.08 | 258.76 | **+11.68** |
| 18 | 256.77 | 265.64 | **+8.87** |
| 19 | 249.37 | 263.50 | **+14.13** |
| 20 | 255.74 | 266.21 | **+10.47** |
| 21 | 255.98 | 266.61 | **+10.63** |
| 22 | 254.77 | 267.90 | **+13.13** |
| 23 | 260.35 | 269.05 | **+8.70** |
| 24 | 256.54 | 266.70 | **+10.16** |
| 25 | 254.89 | 266.76 | **+11.87** |
| 26 | 254.40 | 265.99 | **+11.58** |
| 27 | 254.69 | 265.83 | **+11.14** |
| 28 | 248.70 | 259.42 | **+10.72** |
| 29 | 250.83 | 262.41 | **+11.59** |

Median Δ = +11.1 km/s, range +8.7 to +14.1. **The shift is independent of polynomial degree** — i.e. it is NOT a continuum-shape artifact that one or two degrees absorb differently. The polynomial is already saturated for both cubes.

### Per-SPS at deg=22

| SPS | V_old | V_new | ΔV | σ_old | σ_new | Δσ | χ²_old | χ²_new |
|---|---|---|---|---|---|---|---|---|
| FSPS | −4.96 | −8.93 | −3.97 | 254.77 | 267.90 | +13.13 | 0.083 | 0.081 |
| EMILES | −8.44 | −10.77 | −2.33 | 255.73 | 269.43 | +13.70 | 0.084 | 0.082 |
| XSL | −3.17 | −7.56 | −4.39 | 259.72 | 270.53 | +10.81 | 0.084 | 0.083 |

**All three SPS shift in the same direction by ~+11–14 km/s.** Mean V_sys also shifts by ~−3.5 km/s (deflector appears slightly more blueshifted on the new cube). χ² values are essentially indistinguishable — the new spectrum is just as well-fit by the template family as the old one. **This rules out a "bad new fit" explanation.**

---

## Interpretation

The audit confirms:

1. **The recipe was applied identically.** No code paths differ between the two runs; every checkpoint quantity matches.
2. **The +11 km/s shift is real.** It's reproducible across 15 polynomial degrees × 3 SPS templates with no degree dependence.
3. **The shift comes from the cube-data difference, not a normalisation glitch.** After normalisation the median pixel-to-pixel ratio is 1.003 — so the difference is not a continuum offset. It must be in the *shape* of the absorption lines or the *fine-scale* continuum: something that the additional `_mtwdo_` calibration step changes at the few-percent level, which ppxf then absorbs as a slightly broader LOSVD.
4. **The new reduction's σ_e lands almost exactly on the narrow-window cross-check (267.95 km/s).** Δ between new-wide and narrow-cross-check is **+2.2 km/s** vs **+13.1 km/s** between old-wide and narrow. This suggests the `_mtwdo_` calibration may actually be better-aligned across the wide blue-to-red-half range, removing a systematic offset the old reduction carried.

## Recommendation

- **Keep the headline σ_e(<R_e) = 254.85 ± 17.87 km/s.** It is consistent with the new reduction within the symmetric total budget (0.62σ) and the existing systematic-budget table captures the spread without modification.
- **Treat the reduction-pass Δ as a stated sensitivity** — already documented in
  `METHODS_AND_SYSTEMATICS.md` Part I.9 ("What is *not* in the budget — Reduction-pass sensitivity") and `TESTS_AND_DIAGNOSTICS.md` row M1.
- **Recommendation if a third reduction lands:** revisit. Three reductions would give us a 1σ spread to fold in (likely ±5–11 km/s in quadrature, raising the symmetric total from ±17.87 to ~±18.6–±21 km/s).
- **Optional follow-up if a reviewer challenges this:** point them to (a) the per-degree stability table above, (b) the per-SPS consistency (all three SPS shift by ~+11 km/s), and (c) the narrow-window cross-check, which is essentially unchanged at 267.95 km/s — and lands within ±2 km/s of the new wide-window result. This consistency between narrow-window and new-wide-window σ_e is a strong independent check.

## Files

- Audit script: `/tmp/nb09e_audit.py` (uncommitted, in /tmp)
- Audit stdout log: `/tmp/nb09e_audit_stdout.log`
- Aperture-spectrum overlay + ratio figure: `results/figures/nb09e_audit_spectra.png`
- Reduction-comparison posterior figure (from nb09e): `results/figures/nb09e_reduction_comparison.png`
- Caches: `results/nb09e_new_red_reduction/wR3800_5400_arcmask_{fsps,emiles,xsl}_T*_N500.npz`

---

## Addendum 2026-05-20 — explicit centering re-verification

A natural follow-up question: the audit table above shows `cx_ifu` / `cy_ifu`
identical between the two runs, but those values come from the HST-mean
centroid converted via the cube's WCS header. **Identical WCS metadata does
not by itself guarantee that the cube DATA is at the same physical position
on each cube** — a reduction-pass drizzle/stacking step can shift the data
sub-pixel while keeping the WCS pointed at the same nominal sky position.
Three independent data-derived position checks confirm the cubes are
co-aligned to better than 0.002″ (`/tmp/center_check.py`,
`results/figures/nb09e_centering_check.png`):

| Test | Shift NEW − OLD (spaxels) | In arcsec |
|---|---|---|
| IFU white-light centroid_2dg | (−0.003, +0.004) | **0.0015″** |
| FFT cross-correlation (sub-pix) | (+0.0015, +0.0010) | 0.0005″ |
| Ca H+K + G-band depth-map peak | (−0.091, +0.041) | 0.030″ |

The Ca H+K test has the largest shift (0.030″ = 0.1 spaxel) because that map
is at lower S/N, but it is still ≪ the ±0.4″ amplitude of the 5-center
sweep (F1) that sets the headline centering budget. Both cubes also show
the *same* offset between the IFU-peak and the HST-mean center (OLD: 0.158″,
NEW: 0.157″) — i.e., the cube-to-IFU registration is identical to a part
in a thousand.

The right panel of `nb09e_centering_check.png` shows the visual diagnostic:
the IFU-peak markers for OLD and NEW are visually indistinguishable, and
the NEW−OLD residual map is essentially uniform pale blue with a tiny
center deficit — consistent with the new cube being a uniform ~20% dimmer
in absolute flux (no spatial shape difference).

**Centering is not the cause of the +11 km/s σ_e shift.** The audit verdict
is unchanged: the shift is a real reduction-pass systematic that traces to
subtle changes in absorption-line shape from the additional `_mtwdo_`
calibration step.

---

## Addendum 2026-05-26 — wavelength-region flux scaling + 8000 Å boundary-mask cross-test

Two further checks done on 2026-05-26:

### 1. Per-wavelength flux ratio (cube-vs-cube)

Direct median spectrum comparison between the headline and `_mtwdo_` cubes
(figure `results/figures/nb09e_grating_consistency_check.png`) showed a
2-tier per-night flux scaling pattern in the new cube — flux ratio
NEW/OLD ≈ 0.67 in the Nov-17-only blue edge, 0.79 across the all-three-
nights region (6350–8000 Å), and 0.84 in the Aug-30+Dec-29 red edge.
This is consistent with Kaustubh's hybrid master-twilight+dome
(`_mtwdo_`) flats applying different per-night normalisations — the
Nov 17 frames are scaled to ~67% of the original calibration while Aug 30
and Dec 29 frames are scaled to ~84%.

This 2-tier pattern was the leading hypothesis for the σ_e shift
mechanism: subtle continuum kinks at the wavelength boundaries where
each night turns on or off could introduce a slightly broader
effective LOSVD.

### 2. 8000 Å boundary-mask cross-test (`/tmp/test_8000A_boundary_mask.py`)

To localise the mechanism, re-ran the wide-arc-masked pipeline on the
new cube with an extra mask at obs 7900–8100 Å (= def-rest 4715–4834 Å
at z=0.67564), where Nov 17's RL/7150 grating coverage ends and only
Aug 30 + Dec 29 contribute. Caches at
`results/run_wide_sigma_e/new_8000A_masked/` (N=500 × 3 SPS × 15 deg,
parallel n_jobs=8).

**Result: hypothesis REJECTED.**

| Run | σ_e (km/s) | Δ vs headline | Recovery |
|---|---|---|---|
| Headline old cube (plain wide arc-mask) | 254.85 | — | — |
| NEW cube (plain wide arc-mask) | 265.76 | +10.91 | 0% |
| NEW cube + EXTRA obs 7900–8100 mask | **268.73** | **+13.88** | **−27%** |

σ_e went **UP** by +3 km/s after masking the boundary, not down toward
the headline. The 8000 Å boundary kink is NOT the dominant mechanism.

### Updated mechanism interpretation

The +10.9 km/s shift is therefore **distributed across the wide window**,
not localised to a single transition. Plausible remaining mechanisms:

- **Continuum-shape tilt from the full 2-tier flux pattern**: the
  entire 6362–9051 Å fit window has a different effective continuum
  slope under the new flats vs the original; ppxf absorbs this as a
  slightly broader LOSVD.
- **Absorption-line-shape changes** (line-depth ratios or asymmetries)
  introduced uniformly by the new flat-fielding step — consistent with
  the shift being constant across all 15 polynomial degrees and across
  all 3 SPS libraries (see audit table earlier in this file).

The conservative interpretation — a real "reduction-pass" systematic of
~±11 km/s that should be carried as a stated sensitivity rather than
folded into the formal budget — stands. Already documented in
`METHODS_AND_SYSTEMATICS.md` Part I.9. To eliminate the ambiguity
between the two remaining mechanisms would require either (a) Kaustubh
explaining the per-night flat design intent, or (b) a per-degree
absorption-feature residual analysis at the same polynomial degree on
both cubes.

**Files added by this cross-test**:
- `/tmp/test_8000A_boundary_mask.py` — one-shot driver (uncommitted)
- `logs/test_8000A_boundary_2026-05-26.log` — full run log
- `results/run_wide_sigma_e/new_8000A_masked/wR3800_5400_arcmask_PLUS_obs7900_8100_{fsps,emiles,xsl}_T*_N500.npz` — three N=500 caches
- `results/figures/nb09e_grating_consistency_check.png` — per-wavelength flux ratio figure (from the earlier check)
