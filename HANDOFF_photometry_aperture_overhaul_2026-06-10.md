# HANDOFF — Photometry / deflector-model / aperture overhaul (2026-06-10)

TL;DR of a long iterative arc that started from "are we really using a 2-component Sérsic?"
and ended with a matched-aperture, multi-estimator total stellar mass. **σ_e is untouched
(269.62 ± 13.27 km/s).** The M★ pipeline was substantially reworked.

## The arc of changes (in order, each prompted by a user catch)

1. **2-component → single-component Sérsic deflector model.** The deflector is an elliptical,
   so a single Sérsic is the physically-motivated light model. The old 2-comp (bulge+disk) was
   only a numerical kludge to tame residual-only mask growth. (`scripts/photometry_systematics.py`)

2. **Single-Sérsic fit was BROKEN.** `astropy.Sersic2D.amplitude` is I(r_eff), NOT the central
   peak. Initialising at the peak (n=4) over-predicts the centre ~2000× → LevMar collapses the
   amplitude to ≈0 → the whole galaxy reads as residual → catastrophic over-masking. Caught via
   the BLACK model panels in `photsys_*.png`. The intermediate "M★=11.34 ±0.013" was an ARTIFACT.
   Fix: init ~peak·exp(−b_n), bound amplitude to (peak·1e-3, peak·1.5).

3. **Color gate (HST) / morphology guard (deep JWST) on the IR mask growth.**
   - HST (F200LP, F140W): grow only into pixels **bluer** than the deflector core by ≥0.5 mag in
     F200LP−F140W (blue z=1.302 source vs red passive deflector).
   - Deep JWST (F150W2, F322W2): the color gate is HST-limited where JWST is deepest and was
     MISSING diffuse arc emission (user catch via deep-stretch). Replaced with a **morphology
     guard**: grow into pixels where the source excess beats the model (`resid > MORPH_FRAC·model`,
     i.e. OUTSIDE the Sérsic-dominated region), k=2. Captures the diffuse arc; the deflector body
     (model-dominated) is protected because the single Sérsic holds 66–87% of the flux.

4. **WCS cross-correlation registration.** Raw reprojection left a 0.09–0.18″ HST↔JWST/IR offset
   (JWST i2d GWCS vs FITS-header WCS). `register_shift()` cross-correlates the reprojected F200LP
   arc + color layers to each band → residual 0.00″. Fixed both the arc mis-coverage and a
   companion slipping past the mask edge. (Do NOT anchor on the F200LP centroid — arc-biased.)

5. **Aperture was too small + F200LP mismatched (user catch).** The hand-drawn per-band apertures
   were inconsistent: F200LP a=1.19″ (~0.5 R_e, 17% of light) vs IR a=2.66″ (~1 R_e, 41–48%). The
   tiny F200LP aperture under-measured the blue band → SED artificially red → M★ biased HIGH.
   Provenance: apertures hand-drawn per band in the interactive GUI (`example_outputs/*_params.json`).

6. **Matched 2 R_e aperture + companion masking.** ONE elliptical aperture (a=4.61″, deflector
   PA/axis-ratio, HST-mean centre, ALL bands the same physical region). A 2 R_e aperture encloses
   field companions the old tight aperture excluded → added **companion masking** (DAOStarFinder
   peaks beyond 1.25 R_e, round, not arc, masked with fixed 0.7″ circles; deflector core preserved).
   Default to the **global** arc mask (user preference). 2 R_e chosen over 2.5 R_e: 2.5 captures
   ~8% more light but TRIPLES the companions (2→6) for only 0.036 dex less aperture correction.

## Final numbers (4+1 estimators, matched 2 R_e aperture, 10% flux floor)

| estimator | what it is | log M★ (2 R_e) |
|-----------|-----------|----------------|
| **raw** | empirical aperture, contaminants masked (not filled) | 11.22 |
| **raw + apcorr** | most-empirical total: raw + model WINGS only (masked interior NOT filled) | 11.35 |
| **filled** | + masked pixels filled with single-Sérsic model (within aperture) | 11.36 |
| **total (aperture-corrected)** | filled + model wings beyond aperture — **HEADLINE** | 11.47 |
| **Sérsic-total** | pure model, integrate fitted Sérsic to ∞ | 11.41 |

**Sérsic-only (full-model) systematic budget** (`scripts/sersic_total_systematic.py`):
log M★ ≈ 11.41 ± 0.12 (stat) ± 0.12 (sys) [±0.171 total]. Components (dex): mask 0.102
(dominant), model-form 0.057, fit-n 0.027, flux-floor 0.023, apcorr-recon 0.010.

- raw converges across 2↔2.5 R_e (~11.22); total converges across 1/2/2.5 R_e (~11.45–11.49) —
  validates the aperture correction. Old mismatched-aperture headline 11.36 was biased ~+0.14 high
  by the F200LP mismatch.
- **σ_e unchanged: 269.62 ± 13.27 km/s.** σ_e mask-approach systematic re-derived under the new
  masks = ±4.58 km/s (was ±5.85), still subdominant to the F200-mask w-sweep (±6.65).

## Flux uncertainties (user catch: the 10% floor is far above the real error)

Propagated statistical errors: F200LP 1.58%, F140W 0.41%, F150W2 0.05%, F322W2 0.03% — i.e.
**6–362× below the adopted 10% floor.** The 10% (now reported as such, not as the measurement
error) is a *systematic* floor (zeropoint, aperture, drizzle correlated noise, SED-model mismatch);
~5% is better-motivated. Decision on 5 vs 10% still open.

## Method citations (verified)
- Aperture-flux → Sérsic aperture-correction-to-total feeding M★: **Taylor et al. 2011, MNRAS 418,
  1587** (GAMA `fluxscale`). Our hybrid (empirical aperture + model wings) is the additive variant.
- Strong-lens deflector photometry (Sérsic + mask arc & interlopers → M★): **Sonnenfeld et al. 2013,
  ApJ 777, 97** (SL2S III).
- Sérsic total-magnitude formalism: **Graham & Driver 2005, PASA 22, 118**. SDSS `cmodel` analog.

## Scripts / outputs (deterministic)
- `scripts/photometry_systematics.py` — single-Sérsic + color/morph gate + WCS registration; arc masks.
- `scripts/aperture_2re_companions.py` — matched aperture + companion masks → `aperture_2re_masks.npz`.
- `scripts/aperture_matched_photometry.py` — 4+1 estimators + Bagpipes → `aperture_matched_photometry.npz`.
- `scripts/sersic_total_systematic.py` — Sérsic-only M★ systematic budget → `sersic_total_systematic.npz`.
- `scripts/mask_attempts_comparison.py`, `scripts/mstar_masking_budget.py`.
- Figures: `results/figures/{final_aperture_masks, aperture_2re_companions, mask_attempts_comparison,
  photsys_*, nb08_sersic_vs_aperture}.png`.

## Open items
- Headline M★ value to commit (recommend aperture-corrected total ~11.47, report all estimators).
- Flux-error floor 5 vs 10%.
- Sérsic-only systematic budget (`sersic_total_systematic.py`) — building.
- Companion mask aggressive in deep JWST (harmless for deflector phot; hand-check if paper-bound).
- AXIS_RATIO=0.75, COMP_RGROW=0.7″ are un-tuned sensible defaults.
- **Paper Figs 2 (SED) & 4 (M•–M★) — DONE 2026-06-10.** Regenerated in `AGEL_0206_ApJL_Figures/figures.ipynb`
  (backup `figures_backup_2026-06-10.ipynb`): Fig 4 hardcoded M★ updated to 10^11.47 (+0.09/−0.15);
  Fig 2 SED+inset repointed to `results/bagpipes_sed_results_aperture.npz` (headline total photometry,
  via `scripts/_headline_sed_samples.py`). Verified: Fig 4 point at M★≈3e11, Fig 2 inset 11.47 +0.09/−0.15.
  Figs 1, 3 unaffected.
- **Flux-floor 5% sensitivity — DONE** (`results/aperture_floor5_check.npz`, `scripts/_floor5_check.py`):
  5%↔10% shifts central ≤+0.04 dex (total 11.47→11.50); 10% kept for headline (robust).
- **R_e–masking sensitivity — CHECKED (`scripts/re_mask_sensitivity.py`):** the new color/morph+companion
  masks pull the F140W+F200LP CoG R_e **2.305″→~2.08″** (−0.21″); the Sérsic r_eff confirms the shift
  (real, not a CoG-annulus artifact — the masks remove outer diffuse-arc + interloper light). **Headline
  2.305″ kept** (expert mask, validated, embedded; the 2.08″ has ~30% masked annuli). Impact MINIMAL:
  σ_e ≤~3 km/s, already inside the D7 R_e-source systematic (±6.13, spans 2.168–2.902″); M⋆ *total*
  estimators R_e-robust. Documented in DRAFTING §2.1.2 + CLAUDE.md as the R_e analog of the σ_e/M⋆
  masking systematics.
