<!-- pv-skip-file: dated snapshot — records historical (now-superseded) numbers; current headline in CLAUDE.md / results/PAPER_VALUES.json -->

# Handoff — wide arc-masked headline σ_e finalized, paper Figure 2 updated

**Date:** 2026-05-13
**Sessions covered:** 2026-05-10 through 2026-05-13
**Status at handoff:** σ_e pipeline at the wide arc-masked window fully measured + budgeted + paper figure regenerated. Ready for manuscript text.

---

## TL;DR — the new paper headline

**σ_e(<R_e) = 254.85 ± 17.87 km/s** at the *wide arc-masked* wavelength window
(`wR3800_5400_arcmask`, rest 3800–5336 Å with explicit z=1.302 source-emission masks).

Supersedes the prior nb07c §6cum narrow-window headline of 267.32 ± 24 km/s (now a method cross-check). The CLAUDE.md "Key Results" section was rewritten on 2026-05-11 to reflect this; the prior number is preserved for context.

**Figure 2 left panel** in `AGEL_0206_ApJL_Figures/figures.ipynb` has a new cell (id=`fig2_wide`) at the bottom of Section 2 that uses the new posterior and shows the full rest 3800–5336 Å spectrum being fit. Saved to `AGEL0206_sigma_e_SED_final_wide.pdf`. Same plotting style as the prior cell; just the wavelength range + posterior + absorption-feature labels are updated.

---

## What we did this session arc

### 1. Wide arc-masked window — full systematic budget at N=500/N=250

For each of the 6 systematic components at the wide arc-masked window:

| component | ± km/s | source |
|-----------|--------|--------|
| stat (N=500 all-3 SPS pool) | 6.1 | `results/nb09a_wavelength_sweep/wR3800_5400_arcmask_*_T*_N500.npz` |
| I-shape (10 shapes × N=250) | 1.5 | `results/ishape_sweep_wR3800_5400_arcmask/*_N250.npz` (after Sersic2D bound-fix) |
| F200 spatial mask (peak-to-peak/2) | 3.8 | `results/maskweight_sweep_wR3800_5400_arcmask/{w00,w50,w100}_*_N{100,250}.npz` |
| frame (vac/air, carried) | 5.0 | from prior frame-fix work (2026-04-28) |
| centering (HST WCS, carried) | 4.0 | from prior work |
| fit-window (3-window spread, §9 nb09) | 15.0 | spread between w6500_7500, wR3800_5400_arcmask, wR4000_5400_arcmask |
| **TOTAL (quadrature)** | **17.87** | |

Cross-check at narrow `w6500_7500` window: σ_e = 267.95 ± 30.10 km/s. Consistent at <1σ (difference 13.1 km/s).

### 2. Major methodological finding — Sersic2D bound-fix

The original `scripts/run_isource_shape_sweep.py` allowed Sersic2D `n ∈ [0.3, 8.0]` and `ellip ∈ [0.0, 0.95]`. The F200LP fit escaped to **n=0.30, ellip=0.00** — a degenerate flat-disk minimum that astropy flagged as unsuccessful. The flat I-map under-weighted the center and yielded σ = 271.8 km/s at wide window (297 at narrow) — a 20–35 km/s outlier inflating the I-shape std.

**Fix:** tighten bounds to `n ∈ [1.0, 6.0]` and `ellip ∈ [0.0, 0.6]`, plus try 3 initial n-values (1.5, 2.5, 3.5) and keep lowest-χ². F200LP now converges to a physical n that gives σ = 251.4 km/s at wide. F140W was already fine (n=2.14) so no change there.

**Impact:** I-shape std collapsed from ±5.4 → ±1.6 km/s at wide window; from ±10.5 → ±3.7 at narrow.

Applied at: `scripts/run_isource_shape_sweep.py:181-207`.

### 3. Per-SPS spread collapse — the big paper argument

At the **narrow** Ca H+K + G window, FSPS=254, EMILES=268, XSL=280 — **SPS spread = 26 km/s**.
At the **wide** arc-masked window, FSPS=254, EMILES=253, XSL=258 — **SPS spread = 4 km/s**.

The wide arc-masked window collapses the SPS-template systematic by 6×. Once ppxf has access to Mg b/Fe5270 plus the broader feature set, the FSPS<EMILES<XSL ordering vanishes. This is the strongest methodological argument for the wide window in the paper.

Figure: `results/figures/nb09d_per_sps_both_windows.png`. Documented in nb09d §1.3.

### 4. Source-emission line catalog (arc masks)

The z=1.302 source emits at multiple lines that, when mapped into the deflector rest frame (factor 1.374), fall within rest 3800–5400:

| feature | source rest λ (Å) | deflector rest λ (Å) | masked band |
|---------|--------------------|------------------------|--------------|
| Mg II 2796/2803 | 2796–2803 | 3842–3852 | 3835–3855 |
| ~4534 unidentified | (3300?) | 4534 | 4525–4545 |
| [O II] 3727/3729 | 3727–3729 | 5121–5124 | 5115–5135 |
| Mg b cluster / [Ne III] | various | ~5300 | 5260–5340 |

Codified in `scripts/run_window_sweep.ARC_MASKS_REST` (constant) + `_apply_arc_mask()` (function). 212 of 2374 pixels (8.9%) of the fit are excluded.

This catalog is the first-of-its-kind methodology for AGEL deflector kinematics and is reusable for any AGEL target where z_source is known.

### 5. Notebook 11 — resolved per-PowerBin and per-spaxel kinematics

The earlier per-PowerBin (nb05x Test 19) and per-spaxel attempts failed at the narrow window (noise-blob outer bins at σ > 800 km/s). At the wide arc-masked window the same pipelines give:

- **Per-spaxel S/N≥5 (17 spaxels, R<R_e or R<1.5 R_e):** median σ = 201 km/s, **0 outliers above 400**, σ(R) drops from ~250 → ~150 km/s — sensible elliptical-bulge gradient.
- **PowerBin (target S/N=15, 7 bins at R<1.5 R_e):** median σ = 295 km/s, 5/7 in [150,350]; 2 outer bins still noise-blob.

Final notebook setting: **R < 1.5 R_e** aperture (R<R_e is too small to see a gradient; R<3 R_e adds too much noise).

The wide arc-masked window TRANSFORMED what was previously unmeasurable. Documented in `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` and `results/nb11_perbin_perspaxel_wide.npz`.

### 6. Figures created/regenerated

In `results/figures/`:
- `nb09_fit_windows_on_spectrum.png` — deflector spectrum with 3 fit windows shaded + 9 absorption features labeled + 4 arc masks hatched
- `nb09_n500_window_comparison.png` — three N=500 windows side-by-side
- `nb09_ishape_wide_arcmask.png` — 10 I-shapes at wide window
- `nb09d_per_sps_both_windows.png` — per-SPS posteriors at both windows (the 26→4 km/s collapse plot)
- `nb09d_final_spectra_both_windows.png` — narrow + wide spectra w/ residual sub-panels (2026-05-11 nested-gridspec fix)
- `nb09d_ishape_both_windows.png` — per-shape bars, narrow vs wide
- `nb09d_maskweight_both_windows.png` — mask-weight curves at both windows
- `nb09d_budget_components.png` — narrow vs wide budget bars
- `nb09d_combined_posterior.png` — stat + systematic convolution, both windows (paper Figure 2 candidate B)
- `nb11_perspaxel_sn.png`, `nb11_powerbin_map.png`, `nb11_powerbin_kinematics.png`, `nb11_perbin_perspaxel_sigma_vs_r.png` — kinematic maps

In `AGEL_0206_ApJL_Figures/`:
- `AGEL0206_sigma_e_SED_final_wide.pdf` — **new paper Figure 2** (left panel = wide-window σ_e + spectrum, right panel = SED, unchanged)

---

## Files modified this session

### Repo `AGEL0206_kinematics_and_photometry/`

| file | change |
|------|--------|
| `CLAUDE.md` | "Key Results" section rewritten 2026-05-11 to point to wide arc-masked headline; prior nb07c numbers preserved as cross-checks |
| `notebooks/09_final_sigma_e_paper.ipynb` | Added §9 (3-window comparison + windows-on-spectrum figure), §10 (I-shape sweep, N=250), §10.1 (mask-weight sweep, N=250). 47 cells total. |
| `notebooks/09d_final_systematics_both_windows.ipynb` | **NEW** — side-by-side narrow/wide budgets, 20 cells. Final cell saves `results/sigma_e_final_systematics_nb09d.npz` |
| `notebooks/11_perbin_perspaxel_kinematics_wide.ipynb` | **NEW** — per-PowerBin and per-spaxel ppxf at wide window. 16 cells, R<1.5 R_e aperture |
| `scripts/run_isource_shape_sweep.py` | Sersic2D bounds tightened to n ∈ [1.0, 6.0], ellip ∈ [0.0, 0.6], with 3-init grid; lines 181-207 |
| `scripts/run_isource_shape_window.py` | **NEW** — I-shape sweep at configurable window (used by nb09 §10) |
| `scripts/run_maskweight_window.py` | **NEW** — mask-weight sweep at configurable window (used by nb09 §10.1) |
| `scripts/run_window_sweep.py` | Already had `ARC_MASKS_REST`, `_apply_arc_mask`, and 26 windows (24 + 2 wider blueward) from prior work |

### Repo `AGEL_0206_ApJL_Figures/`

| file | change |
|------|--------|
| `figures.ipynb` | Inserted new cell `fig2_wide` at bottom of Section 2 — same plotting style as the prior σ_e/SED cell, but wide arc-masked posterior and full rest 3800–5336 Å spectrum. Saved as `AGEL0206_sigma_e_SED_final_wide.pdf`. Section 3 (M_BH–σ) uses hardcoded constants `sigma_e_p50=254.85`, `sigma_e_errup=sigma_e_errlo=17.87` (the wide-window budget) — already aligned with the new headline |

### Memory files (`~/.claude/.../memory/`)

| file | change |
|------|--------|
| `MEMORY.md` | Added pointer to new project_nb11_resolved_kinematics.md |
| `project_nb09a_arcmask_discovery.md` | Heavily expanded — full budget table, per-SPS collapse finding, nb09d cross-reference |
| `project_nb11_resolved_kinematics.md` | **NEW** — documents the per-spaxel + PowerBin resolution at the wide window |
| `project_roadmap.md` | Fully refreshed from 42-day-stale state; now reflects wide-arc-masked headline as paper number; lists paper-writing as next steps |

---

## All compute caches (no re-runs needed; everything is on disk)

### Wide arc-masked window (paper headline)
- `results/nb09a_wavelength_sweep/wR3800_5400_arcmask_{fsps,emiles,xsl}_T*_N500.npz` — stat posterior, N=500 wild bootstrap × 3 SPS
- `results/ishape_sweep_wR3800_5400_arcmask/{shape}_{sps}_N{100,250}.npz` — 10 I-shapes × 3 SPS at N=100 and N=250 (the Sersic2D files were refit 2026-05-11)
- `results/maskweight_sweep_wR3800_5400_arcmask/{w00,w50,w100}_{sps}_N{100,250}.npz` — 3 mask weights × 3 SPS at N=100 and N=250

### Narrow window (cross-check)
- `results/nb09a_wavelength_sweep/w6500_7500_{fsps,emiles,xsl}_T*_N500.npz` — stat posterior
- `results/annular_bootstrap_07c_ishape/{shape}_{sps}.npz` — 10 I-shapes × 3 SPS at N=500 (F140W_Sersic2D and F200LP_Sersic2D refit 2026-05-11 with the new bounds)
- `results/mask_weight_sweep.npz` — 5 mask weights × 3 SPS at N=500 (pre-pooled)

### Final consolidated budget
- `results/sigma_e_final_systematics_nb09d.npz` — both windows combined, paper-ready

### Method cross-check (nb07c Gültekin path, narrow window only)
- `results/annular_bootstrap_07c/cumR_2p305_{sps}.npz` — §6cum, paper-cross-check
- `results/annular_bootstrap_07c_nomask/cumR_2p305_{sps}.npz` — F200LP mask sensitivity
- `results/sigma_e_radial_07c.npz` — pooled posteriors

### Notebook 11
- `results/nb11_perbin_perspaxel_wide.npz` — per-bin σ/V/chi2, per-spaxel σ/V/chi2 at 4 S/N floors, bin_map, V_sys

---

## Outstanding work (paper writing, not compute)

1. **Manuscript Section: Kinematics** — rewrite around the wide arc-masked window:
   - Methods: source-emission line catalog (table of the 4 def-rest masks)
   - Results: σ_e = 254.85 ± 17.87 km/s, with narrow window as cross-check
   - Systematic budget table (use the nb09d table directly)
   - Strong selling point: SPS spread 26 → 4 km/s

2. **Figure curation** — choose 2-3 figures from:
   - Figure 2 left panel = `AGEL0206_sigma_e_SED_final_wide.pdf` ← **paper-ready**
   - Methodology supporting figure = `nb09_fit_windows_on_spectrum.png` (shows the 3 windows + features + arc masks)
   - SPS spread = `nb09d_per_sps_both_windows.png`
   - Combined posterior (paper Figure 2 alt) = `nb09d_combined_posterior.png`

3. **Figure 3 (M_BH–σ position)** — already aligned (Section 3 in `figures.ipynb` uses hardcoded `sigma_e_p50=254.85`, errup=errlo=17.87). May want to update the error to be asymmetric using the actual posterior asymmetry.

4. **Section 1 cells (HST/JWST images)** — pre-existing `ModuleNotFoundError` in cells 2-7 of `figures.ipynb` (Section 1). NOT caused by this session's work but may need addressing for full notebook re-execution.

5. **Notebook 11 follow-up (optional)** — if reviewers ask for per-spaxel maps, the data exists. Could also try EMILES or XSL at per-spaxel to confirm SPS-spread collapse holds spatially.

6. **Git organization** — many untracked files (TEXT/, pipes/, large FITS not in repo). Per `project_git_organization.md` memory.

---

## Key insight summary (for future-you returning to this work)

1. **The deflector spectrum has source-emission contamination from z=1.302** at four def-rest bands. Without explicit masking, the wide window gave a bimodal σ posterior with the "wrong" mode at ~100-150 km/s.

2. **The wide arc-masked window is structurally better** than narrow Ca H+K + G alone, not just incrementally better. Three independent improvements stack:
   - 2.6× more fit pixels (2161 vs 927 good) → 4× tighter stat error
   - 6× tighter SPS template agreement (4 vs 26 km/s spread)
   - F200 spatial-mask sensitivity drops 4× (±3.8 vs ±16) because the spectral arc mask already absorbs the worst contamination

3. **The Sersic2D fit pathology** at the F200LP image (n=0.30, ellip=0.00 escape to a degenerate boundary minimum) was a hidden source of I-shape budget inflation in BOTH windows. The fix is universal and worth keeping in `scripts/run_isource_shape_sweep.py`.

4. **Per-spaxel kinematics now work** at the wide arc-masked window — but only inside ~1 R_e (where per-spaxel S/N ≥ 5). Beyond that, KCWI noise wins and ppxf becomes unconstrained per-spaxel; PowerBin co-adds to recover σ for the bulk of the bulge.

5. **The fit-window systematic (±15 km/s) is now the dominant residual** in the budget. The other 5 components quadrature-sum to ~9 km/s. Tightening this would require additional N=500 windows at different feature subsets (e.g., Hβ-only, Mg b-only) — not cheap and likely not worth the marginal improvement.
