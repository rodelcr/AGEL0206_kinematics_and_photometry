<!-- pv-skip-file: dated snapshot — records historical (now-superseded) numbers; current headline in CLAUDE.md / results/PAPER_VALUES.json -->

# Handoff — clean production headline (bad-pix mask + Balmer-unmasked)

**Date:** 2026-05-27
**Sessions covered:** 2026-05-20 through 2026-05-27
**Status at handoff:** Production headline reset on cleaned-reduction pipeline; figures + docs all in sync.

---

## TL;DR — new paper headline

**σ_e(<R_e) = 268.98 ± 18.08 km/s** (symmetric)
**σ_e(<R_e) = 268.98 −18.19 / +17.98 km/s** (asymmetric, preserves stat-side bootstrap skew)

on the wide arc-masked window of the NEW `_mtwdo_` reduction of the
KCWI red-arm cube, with the bad-pixel mask (52 local-MAD-flagged
cosmic-ray residuals) baked into the production pipeline AND the ppxf
default Balmer-absorption mask removed (Hδ, Hγ, Hβ are absorption in
the passive deflector, not emission).

**Source:** `scripts/run_wide_sigma_e.py --cube new_clean --n_bootstrap 500 --n_jobs 8`
**Caches:** `results/run_wide_sigma_e/new_clean/wR3800_5400_arcmask_{fsps,emiles,xsl}_T*_N500.npz`
**Figures regenerated:** `AGEL0206_sigma_e_SED_final_wide.pdf`, `Mbh_sigma_relation.pdf`

## Budget (post-clean, refined 2026-05-26)

```
σ_e(<R_e) = 268.98 ± 18.08 km/s

stat:        ±5.1  (asym -5.47/+4.73)  ← tighter than legacy ±6.1 (+120 Balmer pixels)
I-shape:     ±1.5  (carried from old-cube sweep)
F200 mask:   ±3.8  (carried)
frame:       ±5.0  (carried)
centering:   ±4.0  (carried)
fit-window:  ±15.0 (carried, dominant)
reduction:   ±4.27 ← REFINED post-clean (was ±5.5 pre-clean)
TOTAL sys:   ±17.35
TOTAL:       ±18.08 sym  /  -18.19 / +17.98 asym
```

The reduction-pass ±4.27 is the half-Δ between the two cleaned cubes
(NEW_CLEAN 268.98 − HEADLINE_CLEAN 260.44 = 8.54 → /2 = 4.27).
**Cleaning shrinks the inter-reduction gap** from +10.91 km/s
(legacy 265.76 vs 254.85) → +8.54 km/s (clean 268.98 vs 260.44).

## What's in the fit window

- **Obs 6363–9052 Å** (= rest 3800–5400 Å at z=0.67564)
- 2574 cube pixels, **~2225 used** (≈86%, post-all-masks)
- Excluded (~350 pix = ~14%):
  - **Arc mask (4 bands, 140 rest-Å)**: Mg II (3835–3855), O₂ telluric (4525–4545), [O II] (5115–5135), Mg b+[Ne III] (5260–5340) — `ARC_MASKS_REST`
  - **Bad-pixel mask (26 ranges, 46 rest-Å)**: cosmic-ray residuals at rest 4010–5237 Å (biggest is 6-pix cluster at 5232–5237 Å, 26σ CR spike) — `BAD_PIXELS_REST`
  - **Forbidden lines** ([O III] 4959/5007, ±13 Å each at ±800 km/s) — ppxf default, kept
  - **Template edge padding** (±900 km/s) — ppxf default
- **Kept in fit (Balmer absorption + Lick indices)**: Hδ 4102, Hγ 4340, Hβ 4861, Ca H+K, G-band, Mg b, Fe4383, Fe4531, Fe5270, Fe5335

## Cross-checks at all four (cube × clean) corners

| Cube | clean? | σ_e (km/s) | preset |
|---|---|---|---|
| NEW `_mtwdo_` | clean (headline) | **268.98 −18.19/+17.98** | `--cube new_clean` |
| OLD original | clean (2nd reduction) | 260.44 −18.36/+17.97 | `--cube headline_clean` |
| NEW `_mtwdo_` | legacy (pre-clean) | 265.76 −18.20/+17.59 | `--cube new` |
| OLD original | legacy (pre-clean) | 254.85 −18.20/+17.59 | `--cube headline` |

Narrow Ca H+K + G window cross-check (old cube, nb09d): **267.95 ± 30.10 km/s**.
The cleaned new-cube headline (268.98) lands within 1 km/s of this narrow-window
value — strong consistency across both reductions, both pipelines, and both
window choices.

## What this session arc accomplished

### 1. Discovered the new `_mtwdo_` reduction (2026-05-20)

Kaustubh Gupta provided a new red-arm reduction at
`raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`
using "Hybrid Master-Twilight + Dome" flats (`_mtwdo_`). Same input
frame set as the original Kaustubh stack (verified via Drive folder
structure + co-aligned to <0.002" via 3 independent shift tests).

### 2. Single-pipeline reproducibility script (2026-05-20)

Built `scripts/run_wide_sigma_e.py` — one entry point that runs the
exact same wide-arc-masked σ_e pipeline (load_setup → extract_aperture
→ frame-aware ppxf × 3 SPS × 15 polynomial degrees × N parallel
bootstrap at wR3800_5400_arcmask) on whichever cube via `--cube
<alias>`. Four presets:

- `headline`: legacy OLD cube, arc-mask only
- `new`: legacy NEW _mtwdo_ cube, arc-mask only
- `headline_clean`: OLD cube + bad-pix mask + no-Balmer
- `new_clean`: NEW cube + bad-pix mask + no-Balmer ← CURRENT HEADLINE

Caches at `results/run_wide_sigma_e/{headline,new,headline_clean,new_clean}/`.
Parallel n_jobs=8 with BLAS pinned to 1 thread per worker via
`bootstrap_ppxf_parallel.run_bootstrap_single_degree_parallel`.

### 3. Provenance verification via Google Drive observing logs (2026-05-26)

Authenticated to Google Drive and pulled the Keck PILogin night logs.
Verified definitively that AGEL0206 was observed on **three nights**
(not just one):

| Night (UT) | Program | PI | DESJ0206 frames | RED min | BLUE min |
|---|---|---|---|---|---|
| 2025-08-30 | K409 | TBD | 12 RED `kr250830_00090-00101` + 4 BLUE `kb250830_00052-00055` | 60 | 66 |
| 2025-11-17 | U002 | Jones | 20 RED `kr251117_00129-00148` + 5 BLUE `kb251117_00087-00091` | 100 | 110 |
| 2024-12-29 | U204 | TBD | 4 RED `kr241229_00092-00095` + 1 BLUE `kb241229_00085` | 20 | 22 |
| **TOTAL** | | | **36 RED + 10 BLUE** | **180** | **198** |

(BLUE arm stack actually uses 8 of 10 — the Dec 29 BLUE was excluded
by Kaustubh.) Stale "Sept 29 K409" mention in earlier docs was wrong
— that night had zero DESJ0206 entries. Full per-night frame
inventories + observer lists in `reference_kcwi_data_properties.md`
"Paper-ready boilerplate" section.

### 4. Reduction-pass shift discovery + mechanism narrowing

Discovered the new `_mtwdo_` cube gives σ_e = 265.76 vs original
cube's 254.85 — **+10.91 km/s shift** in σ_e. Ruled out four
hypotheses:

| Hypothesis | Test | Verdict |
|---|---|---|
| Centering shift | 3 independent shift tests | ❌ <0.002″ |
| Localized 8000 Å boundary kink | Mask obs 7900-8100 Å | ❌ σ went UP |
| Bad-pixel outliers (ppxf clean=True) | clean=True at noise×1 | ❌ 0 rejected (noise overestimated) |
| Bad-pixel outliers (local-MAD) | 3σ rolling MAD | ⚠ 52 pix flagged but only ±1.4 km/s on σ_e |

Cleaning both cubes by the same recipe (bad-pix + no-Balmer) brings
their σ_e closer together (gap shrinks 10.91 → 8.54 km/s). The
remaining ~8.5 km/s gap is smooth continuum-shape-tilt from the
per-night flux scaling differences between flat-field calibrations.

### 5. ppxf Balmer-masking investigation (2026-05-26)

User noticed Hβ was in a gray band on Figure 2. Traced to
`ppxf.ppxf_util.determine_goodpixels` which masks all common emission
lines (Hδ, Hγ, Hβ, Hα, plus forbidden lines) at ±800 km/s by default.
For a passive elliptical deflector, the Balmer lines are
**absorption-dominated** — masking them throws away ~120 stellar
absorption pixels in the wide window.

Solution: added `_determine_goodpixels_no_balmer()` to
`scripts/bootstrap_ppxf.py` that keeps only the truly forbidden lines
masked; activated via `mask_balmer=False` parameter on
`setup_ppxf_inputs_from_spectrum`. The `_clean` presets in
`run_wide_sigma_e.py` use this automatically.

Effect on σ_e:
- Balmer-only N=100 test: +1.76 km/s
- Combined with bad-pix at N=500 production: net effect of both = +3.22 km/s

### 6. Bad-pixel mask baked into production

Added `BAD_PIXELS_REST` constant (26 ranges, 52 pixels) to
`scripts/run_window_sweep.py` and matching `_apply_bad_pixels_mask`
helper. Identified via local-MAD outlier detection (rolling 75-pixel
window, |residual|/local_MAD > 3σ) on the canonical residuals.
Replicated on OLD cube (M6): same pixels flagged → bad pixels are
intrinsic to the data, not reduction-specific.

### 7. Tools added this arc

| Script | Purpose |
|---|---|
| `scripts/run_wide_sigma_e.py` | Single-pipeline driver (4 presets, --compare mode) |
| `scripts/bootstrap_ppxf.py` (extended) | `mask_balmer` parameter + `_determine_goodpixels_no_balmer` |
| `scripts/run_window_sweep.py` (extended) | `BAD_PIXELS_REST` + `_apply_bad_pixels_mask` |
| `notebooks/09e_new_red_reduction.ipynb` | Initial nb09e parallel-cube smoke notebook |
| (uncommitted /tmp/ scripts kept for reproducibility) | All cross-tests: 8000Å boundary, ppxf clean=True, local-MAD clip, OLD-cube replication, no-Balmer test, noise-scaling scan |

## Files modified this arc

### Code
- `scripts/run_window_sweep.py` — added BAD_PIXELS_REST (26 ranges, 52 pix) + `_apply_bad_pixels_mask`
- `scripts/bootstrap_ppxf.py` — added `_determine_goodpixels_no_balmer` + `mask_balmer` parameter on both `_prep_spectrum_for_ppxf` and `setup_ppxf_inputs_from_spectrum`
- `scripts/run_wide_sigma_e.py` — 4 presets (`headline`, `new`, `headline_clean`, `new_clean`); `_clean` auto-applies bad-pix + no-Balmer; SYS_QUAD updated to 17.35

### Docs
- `METHODS_AND_SYSTEMATICS.md` — §0 headline table (cleaned-production values), Part I.9 budget table, Part I.9 "What is not in the budget" expanded with M1/M5/M6/M7, Part III.1 final tabulated systematics, Part III.5 added item #0 (Hδ TODO)
- `TESTS_AND_DIAGNOSTICS.md` — §0 headline values, M5 (local-MAD), M6 (OLD-cube replication), M7 (noise-scaling scan)
- `CLAUDE.md` (project) — Key Results section rewritten with cleaned-production headline + 7-component budget
- `reference_kcwi_data_properties.md` (memory) — Full 3-night provenance table + Drive ID inventory + paper-ready boilerplate
- `project_nb09e_reduction_systematic.md` (memory) — Promoted to headline + addenda (b) 8000Å test rejected, (c) ppxf clean=True null result

### Figures regenerated
- `AGEL_0206_ApJL_Figures/figures.ipynb` cells: `fig2_wide` (loads from new_clean caches via fig2_data.npz), `d55da320` (Fig 3 M_BH-σ, AGEL0206 point at σ=268.98)
- `AGEL0206_sigma_e_SED_final_wide.pdf` — title shows σ_⋆(<R_e) = 269⁺⁵₋₅(stat) ± 17(sys) km/s; Lick Fe labels added (Fe4383, Fe4531, plus Fe5270, Fe5335 now visible inside masked bands via fallback y)
- `Mbh_sigma_relation.pdf` — AGEL0206 point at σ=268.98, asym errs -18.19/+17.98

### Caches added (all under `results/run_wide_sigma_e/`)
- `headline/` and `new/` — legacy un-cleaned (existed before this arc)
- `headline_clean/` and `new_clean/` — NEW post-clean production (3 SPS × N=500 each + figure2_data.npz in new_clean/)
- `new_local_mad_clip_N100/`, `old_local_mad_clip_N100/` — local-MAD cross-test caches
- `new_sigma_clip_N100/` — ppxf clean=True null-result cache
- `new_noise_scaled_clean_N100/` — noise-scaling cross-test cache
- `new_8000A_masked/` — 8000Å boundary cross-test cache
- `new_no_balmer_mask_N100/` — Balmer-only N=100 test cache

## Outstanding work / open items

(see `METHODS_AND_SYSTEMATICS.md` Part III.5 for the full list)

1. **Hδ targeted treatment (TODO 2026-05-26)** — Hδ may still need
   narrow/partial masking; re-run local-MAD on no-Balmer-mask fit
   residuals to verify. Inline comment + Part III.5 item #0 flag it.

2. **K409 PI** for the Aug 30 night needs confirmation (paper
   acknowledgements / observing-program citation).

3. **Per-night DIMM seeing** for the K409 Aug 30 night needs
   extraction from observing-log header (airmass 1.07–1.08 verified;
   GUIDFWHM=1.27" from Nov 17 used as conservative estimate).

4. **Optional: Balmer-only N=500 production isolated** — current N=500
   `new_clean` runs both bad-pix + no-Balmer. An N=500 Balmer-only-no-
   bad-pix run would isolate the Balmer effect at production stats
   (currently only N=100 isolated test: +1.76 km/s).

5. **Optional: re-run I-shape, mask-weight, fit-window sweeps on the
   new cube** — currently budget components I-shape ±1.5, mask ±3.8,
   window ±15 are carried from the old-cube sweeps. Strict re-derivation
   would cost ~50 min × 3 sweeps. Flagged but not blocking.

6. **Manuscript kinematics section** — still needs rewriting around
   the new cleaned headline with the new reduction-pass narrative.

## Key insights for future-you returning to this work

1. **The wavelength-window systematic (±15 km/s) is now the dominant
   residual** in the budget. Reduction-pass (±4.27), frame (±5),
   centering (±4), F200 mask (±3.8), I-shape (±1.5) all quadrature-
   sum to ~9 km/s — combined are smaller than the single window
   component.

2. **Cleaning works on both cubes the same way** (M6 replication):
   the 52 bad pixels are intrinsic to the data, not reduction-specific.
   So if any future re-reduction lands, applying the same clean
   recipe will give a self-consistent comparison.

3. **The new reduction is genuinely better** at consistency across
   methodology choices — cleaned σ_e = 268.98 lands within 1 km/s
   of the narrow-window cross-check (267.95). Old cube was ~13 km/s
   discrepant before clean and ~7.5 km/s discrepant after clean.

4. **ppxf's default `clean=True` is useless on this data** because
   the noise array is overestimated (median residual / noise << 1).
   Local-MAD-based outlier detection on rolling residuals is the
   correct approach for cosmic-ray rejection here. Caveat for future
   AGEL targets that share the kcwiRedux pipeline.

5. **For PASSIVE elliptical deflectors specifically, do NOT use
   `ppxf.ppxf_util.determine_goodpixels` as-is** — it masks Balmer
   absorption as if it were emission. Use `mask_balmer=False` on
   `setup_ppxf_inputs_from_spectrum` (or pass a custom goodpixels
   list to ppxf directly).

6. **The KCWI observing logs (Keck PILogin) ARE machine-readable** —
   the autolog PDFs in the Drive's "Logs by date" tree have all the
   per-frame metadata (Object, RA/Dec, EXPTIME, airmass) extractable
   via Drive's text-content read. Used for the 3-night provenance
   verification 2026-05-26.

## Reproduction recipe (one command)

```bash
conda activate ISMgas

# Single-command headline reproduction (~13 min at n_jobs=8):
python scripts/run_wide_sigma_e.py --cube new_clean --n_bootstrap 500 --n_jobs 8

# Full --compare mode (runs both cleaned cubes, ~26 min):
python scripts/run_wide_sigma_e.py --compare --n_bootstrap 500 --n_jobs 8

# Regenerate figure2_data.npz after a re-run:
python /tmp/prep_fig2_new_clean_data.py

# Re-render Fig 2 + Fig 3 (cells 0, 19, 22, 23, 24, 28, 29):
python /tmp/run_fig2_fig3.py
```
