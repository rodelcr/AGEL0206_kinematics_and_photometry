# Handoff — §6cum I(r) sweep × {masked, no-mask}

**Goal:** Run §6cum cumulative I-weighted ppxf at R<R_e for **6 I-source choices × 2 masking states (F200LP-masked, no-mask) = 12 σ_e measurements**, all at N=500 bootstrap × 3 SPS × 15 polynomial degrees. Compare to determine the relative magnitude of the I(r) systematic vs the masking systematic on the headline σ_e.

**Why this matters:** We've measured the masking systematic (Δσ_e = −16.4 km/s for mask on/off) and the centering systematic (3.7 km/s spread) directly, but the I(r) systematic for §6cum has not been measured. Until it is, we cannot say definitively whether the headline 267 ± 24 km/s is robust against the I(r) choice in the cumulative path. nb07a measured this for §7 (annular) and got ≤ 3–8% spread, but §6cum was never tested.

---

## Current state of the σ_e analysis (paper headline numbers, as of 2026-04-28)

- **σ_e(<R_e) = 267.32 ± 24 km/s** — F200LP-masked, §6cum cumulative I-weighted ppxf, IFU 6500–7500 Å I-weight (current paper headline)
- **σ_e(<R_e) = 250.96 ± 23 km/s** — no arc mask, same I-weight (sensitivity test, Δ = −16.4 km/s)
- σ_e robust to ±0.4″ center shifts (3.7 km/s spread across 5 candidate centers)
- Per-SPS at R<R_e: FSPS 253, EMILES 268, XSL 280 km/s (combined 267 ± 24)
- §7 cross-checks (255, 271, 386 unfiltered) are settled — see `reference_cumulative_vs_annular_sigma_e.md`

**The mask on/off shift (−16.4 km/s) is the largest systematic identified so far on §6cum.** It's about 68% of the SPS budget. The centering systematic is much smaller (15% of SPS budget). The I(r) systematic is the missing piece.

---

## What's failed so far (lessons learned — DO NOT repeat)

I tried three from-scratch implementations of an I(r) sweep, all crashed:

1. **First attempt (6 I-maps with 2D Sersic):** Crashed because my Sersic2D fit had a shape mismatch — I used `np.zeros_like(img140_full, dtype=float)` to receive both F140W and F200LP Sersic models, but the two HST cutouts have *different* pixel scales (0.0298″ vs 0.0186″) and different array shapes. The Sersic-to-IFU reprojection then tried to broadcast into a 0×0 sub-slice.

2. **Second attempt (4 I-maps without Sersic):** Skipped 2D Sersic but the 1D-CoG-based F140W and F200LP I-maps had boundary issues. My naive `interp1d(rs, I, fill_value=(I[0], 0.0))` produced near-zero weights at the aperture edges, which propagated into the I-weighted aperture sum and gave NaN/overflow inside ppxf.

3. **Pre-existing nb07a code works.** nb07a fit Sersic2D *separately* on F140W and F200LP cutouts (different `sub_shape`, `x0`, `y0` per band), with `weights = (~mask).astype(float)`. nb06 has a `build_I_weight_map(source, mask_strategy)` function that already handles all 6 I-sources × 3 mask strategies. **THE RIGHT MOVE IS TO CALL nb06's tested function, not reimplement.**

---

## Plan for the re-implementation

### Step 1: Reuse nb06's `build_I_weight_map` function (cell 31)

The function in `notebooks/06_final_sigma_Re_apertures.ipynb` cell 31 takes `(source, mask_strategy)` and returns a (ny, nx) intensity map on the IFU spaxel grid. It supports:

| source | mask_strategy options |
|---|---|
| `IFU_wl` (full white-light) | `unmasked`, `15pct_psf` (PSF-aware contamination map) |
| `IFU_band` (6500-7500 Å) | `unmasked`, `15pct_psf` |
| `F140W` | `unmasked`, `arc_only_ifu`, `hst_mask_excl` |
| `F200LP` | `unmasked`, `arc_only_ifu`, `hst_mask_excl` |

The HEADLINE choice in nb06 is `('IFU_band', 'unmasked')`. The 8 viz combos used in nb06 cell 32 are documented as the standard comparison set:
```python
[('IFU_band','unmasked'), ('IFU_band','15pct_psf'),
 ('IFU_wl',  'unmasked'), ('IFU_wl',  '15pct_psf'),
 ('F140W',   'unmasked'), ('F140W',   'arc_only_ifu'),
 ('F200LP',  'unmasked'), ('F200LP',  'arc_only_ifu')]
```

For the §6cum I(r) sweep, **use these same 8 I-source/mask combos** — that gives a direct 6-vs-7 reproducibility check (apples to apples with nb06's PowerBin Gültekin sum).

### Step 2: Wire `build_I_weight_map` into a §6cum pipeline

The §6cum recipe (from nb07c cell 31):
```python
sel = (r_spax < R_e) & ~arc_spax_mask  # r_spax computed from HST_mean center
w = I_map[sel]              # I_map from build_I_weight_map(source, mask_strategy)
w_norm = w / w.sum()
flux_agg = np.sum(cube[:, sel] * w_norm[None, :], axis=1)
# → ppxf bootstrap on flux_agg
```

For the **no-mask** half of the test, drop `& ~arc_spax_mask` from `sel`.

But note: when `mask_strategy = 'arc_only_ifu'` or `'15pct_psf'`, the I_map ALREADY zeroes the arc region. So the mask × I-map combinations to actually run are:

| I-map | mask strategy on sel | Actually meaningful? |
|---|---|---|
| `IFU_band, unmasked` | `~arc_spax_mask` | YES — current headline |
| `IFU_band, unmasked` | no mask on sel | YES — no-mask sensitivity (already done) |
| `IFU_band, 15pct_psf` | `~arc_spax_mask` | redundant (mask both ways) |
| `IFU_band, 15pct_psf` | no mask on sel | NEW — PSF-aware mask via I-weight only |
| `F140W, unmasked` | `~arc_spax_mask` | YES |
| `F140W, unmasked` | no mask on sel | YES |
| `F140W, arc_only_ifu` | no mask on sel | YES — same arc removal but via I=0 |
| `F200LP, unmasked` | `~arc_spax_mask` | YES |
| `F200LP, unmasked` | no mask on sel | YES |
| `F200LP, arc_only_ifu` | no mask on sel | YES |

So the actual minimum-non-redundant matrix is ~10–12 distinct measurements.

### Step 3: Cache structure

Use `results/annular_bootstrap_07c_isource/` (already created, currently has corrupt files from failed runs — wipe first). One npz per `(I_source, mask_strategy, sel_mask, sps)` tuple. Filename convention:
```
results/annular_bootstrap_07c_isource/{Isrc}_{maskstrat}_{selstrat}_{sps}.npz
```
where `selstrat ∈ {arcmask, nomask}`.

Each npz contains: `V_orig, sig_orig, V_boot (15, 500), sig_boot (15, 500), degrees (15,), r_max, sps, isource, maskstrat, selstrat`.

### Step 4: Output

Final summary table:
```
I_source         mask_strategy   sel_strategy   σ_e p50   −1σ   +1σ   N_kept   Δ vs headline
IFU_band         unmasked        arcmask        267.3     24.4  24.1  147     0.00 (headline)
IFU_band         unmasked        nomask         251.0     23.0  23.4  185    -16.4 (mask test)
IFU_band         15pct_psf       nomask         ?         ?     ?     ?      ? (PSF-aware mask via I)
F140W            unmasked        arcmask        ?
F140W            unmasked        nomask         ?
F140W            arc_only_ifu    nomask         ?
F200LP           unmasked        arcmask        ?
F200LP           unmasked        nomask         ?
F200LP           arc_only_ifu    nomask         ?
... etc
```

Plot: 4-panel figure
- (A) Bar chart of σ_e per (I_source, mask_strategy, sel_strategy) with combined-SPS error bars; group by mask state to highlight on-vs-off effect per I-source
- (B) σ_e distribution histogram overlay for ~3–4 representative cases
- (C) The I-maps themselves on the IFU grid (small multiples) for visual sanity
- (D) The aperture spectra for ~3 representative cases overlaid

### Step 5: Decision criteria

Once the matrix is filled in:

| Outcome | Implication |
|---|---|
| All 6 I-sources give σ_e within 5 km/s of each other | I(r) is sub-dominant; quote single headline |
| Spread is 5–10 km/s | I(r) is comparable to centering; add to systematic budget in quadrature |
| Spread is > 10 km/s | I(r) is comparable to or larger than the SPS systematic; need a justified single choice for headline |
| `15pct_psf` (PSF-aware) gives different σ_e from `arc_only_ifu` (raw F200) by > 5 km/s | PSF-broadening of the arc mask matters; should use nb06's PSF-aware mask, not nb07c's raw F200 reprojection |

---

## Implementation skeleton (use this, don't rewrite)

```python
"""§6cum I(r) sweep × mask sensitivity at R<R_e (N=500).

Uses nb06's build_I_weight_map (extracted to script form) and nb07c's §6cum
flux extraction recipe. Saves per-config npz to results/annular_bootstrap_07c_isource/.
"""
import os, sys, time
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates, gaussian_filter
from photutils.centroids import centroid_2dg

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

os.chdir('/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry')
sys.path.insert(0, '.')
from ppxf.ppxf import ppxf
from scripts.bootstrap_ppxf import setup_ppxf_inputs_from_spectrum, NOISE_SLICE
from scripts.bootstrap_ppxf_parallel import run_bootstrap_single_degree_parallel

# 1. Load IFU + center it (use HST_mean from F140W+F200LP centroid_2dg, same as nb07c)
# 2. Build PSF-aware contam map (nb06 cell 5):
#    spaxel_contam = gaussian_filter(F200_mask, sigma=KCWI_seeing_FWHM/HST_pix/2.355)  reproject to IFU
#    spaxel_masked = spaxel_contam > 0.15
# 3. Reproject F200LP raw mask (order=0): arc_spax_mask  (nb07c convention)
# 4. EXTRACT build_I_weight_map function from nb06 cell 31 verbatim
#    Modify only: change the 'F140W' source to read from your local cutouts,
#    point HST_CUTOUT/HST_MASK to the same paths nb07c uses
# 5. For each (Isrc, maskstrat) combo + each selstrat ∈ {arcmask, nomask}:
#       I_map = build_I_weight_map(Isrc, maskstrat)
#       sel = (r_spax < R_e) & (~arc_spax_mask if selstrat=='arcmask' else True)
#       w = I_map[sel]; w_norm = w / w.sum()
#       flux_agg = np.sum(cube[:, sel] * w_norm[None, :], axis=1)
#       For sps in [fsps, emiles, xsl]:
#          inp = setup_ppxf_inputs_from_spectrum(flux_agg, noise_sky, hdr, sps, z=0.67564)
#          For each deg in arange(15, 30):
#             pp = ppxf(...); rb = run_bootstrap_single_degree_parallel(..., n_jobs=8)
#          save to results/annular_bootstrap_07c_isource/{Isrc}_{maskstrat}_{selstrat}_{sps}.npz
# 6. Build summary table + plot
```

**Estimated wall time:** 12 configurations × 3 SPS × 15 deg × N=500 bootstrap with `n_jobs=8` parallelization. Each (config, SPS) takes ~2–3 min based on §6cum-nomask runtime. Total: ~75–110 min wall.

**Pre-flight checks before the long run:**
- Sanity-check `build_I_weight_map(Isrc, maskstrat)` output for each combo: print N_nonzero, peak location, sum-within-R_e — should match nb06 cell 32's printed values
- Confirm the I-weighted aperture spectra at R<R_e have S/N_band ≥ 5 for all configurations (the failed earlier runs had S/N=3.8 which is borderline; if any config drops below S/N=3, that config will likely give NaN ppxf)
- Run a smoke test at N=10 first (~2 min) to catch any bugs before the 1.5h N=500 run

---

## Other open items from the session (not blocking the I(r) sweep but worth flagging)

### Centering investigation — RESOLVED
- All 5 candidate centers (HST_mean, IFU_WL_peak, F140W_only, Sersic2D_F140W, IFU_peak_arcmask) agree at < 0.4″ and give σ_e within 3.7 km/s spread (264.7–268.4)
- Headline 267 ± 24 stands; centering is sub-dominant
- See `NOTES_centering_investigation_2026-04-27.md` for the full investigation
- Cache: `results/annular_bootstrap_07c_centers/` (5 centers × 3 SPS × N=500)

### F200LP mask sensitivity — RESOLVED
- Mask on: 267.32 ± 24 km/s; mask off: 250.96 ± 23 km/s; **Δ = −16.36 km/s** (~68% of SPS budget)
- The mask is real and quantifiable, not negligible
- Should be quoted as a sensitivity in the paper, with 267 ± 24 as fiducial
- See `NOTES_methodology_2026-04-27.md` and the diagnostic figure `results/figures/nb07c_s6cum_nomask_diagnostic.png`

### KCWI cube provenance — RESOLVED
- The cube used (`Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`, 253 MB, 3317×100×100) IS the final reduced product. Only its FITS headers are mislabeled (DATE-OBS=2025-11-17, PROGNAME=U002, PROGPI=Jones describes only the last stacking pass)
- Multi-night data: Aug 29 2025 K409 (canonical) + Sept 29 K409 + Dec 29 2024 + Nov 17 U002
- For the paper, cite the multi-night provenance, not the header
- See `HANDOVER.md` and `~/.claude/.../memory/reference_kcwi_data_properties.md`

### PSF-aware mask comparison — STILL TODO (separate from I(r) sweep)
- nb06 used `spaxel_contam = gaussian_filter(HST_mask, sigma=PSF_HST_pix) > 0.15` — a KCWI-seeing-broadened mask
- nb07c uses raw HST F200LP reprojected with order=0 nearest-neighbor — narrower than the seeing-convolved one
- The I(r) sweep above includes both `15pct_psf` and `arc_only_ifu` strategies, so this gets answered as a side effect
- If the I(r) sweep finds 15pct_psf vs arc_only_ifu σ_e differs by > 5 km/s, the PSF-broadening matters and we should adopt nb06's mask going forward

### Gültekin formula attribution — RESOLVED
- Original Gültekin 2009 eq. 1 is 1D longslit form; our 2πr is the 2D IFU extension matching SAURON+ATLAS3D+MaNGA convention. Memory file fixed.

### Method choice (§6cum vs §7) — RESOLVED
- §6cum is best (7 specific reasons; no binning to defend, single LOSVD fit, etc.)
- §7 with equal-N is the cross-check; never quote §7 as headline
- Documented across `reference_cumulative_vs_annular_sigma_e.md`, `project_sigma_e_gultekin.md`, `CLAUDE.md`

---

## Recent commits this session

- `0cce2cd` — Centering investigation: σ_e robust to ±0.4″ center shifts
- `6fde707` — Add masked-vs-nomask diagnostic + dual-quote σ_e in figures.ipynb
- `6283d98` — Add §6cum-nomask sensitivity test: F200LP mask is −16 km/s effect
- `fd19e27` — HANDOVER.md (compaction state)

---

## Resources

### Code paths
- `notebooks/06_final_sigma_Re_apertures.ipynb` cell 31 — `build_I_weight_map(source, mask_strategy)`
- `notebooks/07c_sigma_e_equalN.ipynb` cells 31–32 — §6cum aperture extraction + bootstrap (current headline pipeline)
- `scripts/bootstrap_ppxf.py` — `setup_ppxf_inputs_from_spectrum`, `NOISE_SLICE`
- `scripts/bootstrap_ppxf_parallel.py` — `run_bootstrap_single_degree_parallel(..., n_jobs=8)`

### Data
- IFU cube: `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` (final reduction, mislabeled headers)
- HST F140W: `../velocity_dispersion_from_IFU/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits` + `_mask.fits`
- HST F200LP: `..._F200LP_..._cutout_L3.fits` + `_mask.fits`

### Caches
- `results/annular_bootstrap_07c/` — masked headline (used)
- `results/annular_bootstrap_07c_nomask/` — no-mask sensitivity (used)
- `results/annular_bootstrap_07c_centers/` — center comparison (used)
- `results/annular_bootstrap_07c_isource/` — **EXISTS BUT CORRUPT — wipe before re-running**

### Figures
- `results/figures/nb07c_s6cum_nomask_sensitivity.png` — mask comparison plot
- `results/figures/nb07c_s6cum_nomask_diagnostic.png` — 4-panel diagnostic (spaxel maps + spectra + posterior)
- `results/figures/center_check.png` and `center_candidates.png` — centering diagnostics

### Methodology notes
- `NOTES_methodology_2026-04-27.md` — full session methodology summary
- `NOTES_centering_investigation_2026-04-27.md` — centering investigation
- `HANDOVER.md` — compaction-state document
- `~/.claude/.../memory/reference_cumulative_vs_annular_sigma_e.md` — §6cum is best + binning judgment
- `~/.claude/.../memory/reference_gultekin_implementation.md` — exact discrete formula + 2D/1D attribution
- `~/.claude/.../memory/reference_sps_systematic.md` — SPS pooling math
- `~/.claude/.../memory/project_sigma_e_gultekin.md` — full methodology narrative
- `~/.claude/.../memory/project_nb07e_arc_subtraction.md` — nb07e residual-arc test (only tests *residual*, not bulk mask sensitivity)

---

## Bottom line for the next session

**The job:** Run the 8–12 configuration §6cum I(r)-mask sweep using nb06's tested `build_I_weight_map`. Wall time ~1.5h N=500 parallelized. Result fills the last gap in the §6cum systematic budget for the paper.

**The trap:** Do NOT reimplement Sersic2D, CoG, or any I-map machinery from scratch. nb06 already has them tested. Import or copy verbatim.

**What success looks like:** A 8–12-row table of σ_e(<R_e) values + 1σ widths showing the I(r) spread relative to the −16 km/s mask shift. If the I-source spread is < 10 km/s, headline 267 ± 24 is fully validated. If > 10 km/s, the systematic budget needs revising.
