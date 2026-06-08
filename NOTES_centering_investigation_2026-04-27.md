<!-- pv-skip-file: dated snapshot — records historical (now-superseded) numbers; current headline in CLAUDE.md / results/PAPER_VALUES.json -->

# Centering investigation — 2026-04-27

Triggered by user noticing the cyan + (deflector center) appeared off-center from the IFU brightness in the §6cum-nomask diagnostic figure. This note captures the investigation, findings, candidate centers tested, and σ_e sensitivity to center choice.

---

## 1. The puzzle

The §6cum-nomask diagnostic figure showed the deflector center (cyan +) sitting at (0, 0) but the bright region of the IFU 6500–7500 Å image visually peaking ~0.5–0.6″ to the upper-right.

Two possible explanations:
- (a) Visual artifact — arc adds bright UV continuum on one side of the deflector, biasing the eye
- (b) Real disagreement between HST-derived center (used by nb07c) and IFU brightness peak

**Both are partly true** — see §3 below.

---

## 2. What each notebook uses for the center

| Notebook | Center method | Code |
|---|---|---|
| **nb06** | **Smoothed IFU white-light peak** | `cy, cx = unravel_index(argmax(gaussian_filter(sum(cube, axis=0), sigma=2)), shape)` |
| **nb07** | HST F140W + F200LP centroid_2dg with their masks → mean | `find_center_centroid(...)` |
| **nb07a** | Same as nb07 (`find_center` function, identical logic) | uses RA_DEFL, DEC_DEFL = 31.55611, -1.23817 as init |
| **nb07c** (HEADLINE) | Same as nb07 / nb07a | identical `find_center` function |

**Switching from nb06 (IFU-WL) to nb07 family (HST centroid) was a methodology change that introduced a center disagreement.** nb06 had this aspect "solved" by using the IFU brightness directly.

The `find_center` function has a soft check `print(f'F140W vs F200LP centroid |Δ|={d_tot:.3f}"  → {"PASS" if d_tot<0.2 else "FAIL"}')` — but does NOT abort on FAIL.

With the current data, F140W vs F200LP centroid disagreement is **0.473″** — fails the 0.2″ check.

---

## 3. Five candidate centers compared

Computed with photutils.centroid_2dg (HST methods) and gaussian-smoothed argmax (IFU methods):

| Center name | RA (deg) | Dec (deg) | IFU pix | Δ vs HST_mean |
|---|---|---|---|---|
| **HST_mean** (nb07c current) | 31.556154 | −1.238184 | (49.80, 54.31) | 0.000″ |
| **IFU_WL_peak** (nb06) | 31.556137 | −1.238127 | (50.00, 55.00) | 0.214″ |
| **F140W_only** | 31.556097 | −1.238149 | (50.49, 54.73) | 0.240″ |
| **Sersic2D_F140W** | 31.556132 | −1.238103 | (50.07, 55.29) | 0.302″ |
| **IFU_peak_arcmask** | 31.556137 | −1.238210 | (50.00, 54.00) | 0.111″ |

All 5 agree at < 0.4″. Pairwise max separation: 0.386″ (Sersic2D ↔ IFU_peak). Sub-IFU-spaxel level (0.30″/spaxel).

So contrary to the original visual impression, the centers actually agree well. The "off-center" appearance in the diagnostic figure was a combination of:
- The arc adding bright UV continuum on one side of the deflector
- Asymmetric bright IFU region extending from arc-mask region into the lower-right (likely PSF-broadened arc light not caught by the HST-tight F200LP mask)
- My diagnostic plot using the slightly stale CLAUDE.md hardcoded RA/Dec (0.087″ off from HST mean)

---

## 4. The PSF-broadening issue (separate finding)

Looking at the IFU 6500–7500 Å image with 5 candidate centers overlaid, the bright IFU region extends from the F200 mask region (upper-right) into the lower-right of R_e. The lower-right bright lobe is NOT in the F200 mask region.

Two interpretations:
- **(i) Real deflector asymmetry** — unlikely; AGEL0206 is described as an elliptical bulge with mild ellipticity (e ~ 0.13 from Sersic fit)
- **(ii) PSF-broadened arc light** — the F200LP mask was built on HST (PSF ~ 0.15″); the IFU has KCWI seeing FWHM 1.27″, so arc light gets smeared into spaxels OUTSIDE the F200-mask boundary

**nb06 addressed this with its PSF-aware contamination map**: 
```python
psf_sigma_hst_pix = (KCWI_SEEING_FWHM / HST_PIXSCALE) / 2.355
mask_convolved = gaussian_filter(mask_hst.astype(float), sigma=psf_sigma_hst_pix)
spaxel_contam = reproject_to_IFU(mask_convolved)
spaxel_masked = spaxel_contam > 0.15
```

This produces a wider, KCWI-seeing-convolved mask. nb06 used this as the I-weight mask. **nb07c uses the raw HST-resolution F200LP reprojected with order=0 nearest-neighbor — narrower than the seeing-convolved one.**

Implication: nb07c's "F200LP-masked" aperture may admit some seeing-broadened arc light that nb06's PSF-aware mask would catch. This is an additional systematic on top of the F200-mask-on/off sensitivity test (which gave Δσ_e = −16 km/s).

Cross-check: re-run §6cum at R<R_e using nb06's PSF-aware mask (15% threshold) and compare to the F200 hard mask. Expected effect: ≤ 5 km/s on σ_e (since most of the arc IS in the F200 mask; only the seeing wings differ).

---

## 5. σ_e-per-center comparison RESULTS (N=500, completed 2026-04-27 22:33)

Computed σ_e(<R_e) for each of the 5 candidate centers using same F200LP arc mask + IFU 6500–7500 Å I-weight + N=500 bootstrap × 3 SPS × 15 polynomial degrees. Wall time 28.7 min. Cache: `results/annular_bootstrap_07c_centers/`.

| Center | σ_e(<R_e) p50 | −1σ | +1σ | N_kept | Δ vs HST_mean |
|---|---|---|---|---|---|
| **HST_mean (nb07c headline)** | **266.41** | 23.14 | 22.90 | 147 | 0.00 km/s |
| IFU_WL_peak (nb06) | 268.42 | 24.76 | 23.57 | 148 | +2.01 km/s |
| F140W_only | 268.11 | 24.71 | 23.40 | 148 | +1.70 km/s |
| Sersic2D_F140W | 264.70 | 25.27 | 23.41 | 149 | −1.71 km/s |
| IFU_peak_arcmask | 266.40 | 24.19 | 22.31 | 147 | −0.01 km/s |

**Total spread: 3.72 km/s (264.70 → 268.42)** — ~15% of the SPS-systematic budget (±24). **Centering at the 0.3″ level is NOT a dominant systematic** on σ_e.

### Verdict
- The headline σ_e(<R_e) = 267 ± 24 km/s stands — robust to ±0.4″ center shifts
- nb07c's HST_mean choice is fine. Switching to IFU_WL_peak (nb06's choice) would shift σ_e by ~+2 km/s, well below 1σ
- HST_mean and IFU_peak_arcmask are essentially equivalent (Δ = 0.01 km/s) — both centers round to adjacent spaxels but the I-weighted aperture is dominated by the bright deflector core regardless
- Sub-km/s difference between this run's HST_mean (266.41) and the earlier nb07c headline (267.32) is bootstrap noise — different seed offsets produce different draws of the same underlying posterior

### Recommendation
Continue using HST_mean (nb07c's existing choice) for the headline. **Do NOT re-run the σ_e analysis with a different center** — the answer doesn't change at the level that matters.

---

## 6. Sersic2D fit sanity check (nb07a-style)

Re-ran the Sersic2D fit with weights = (~F140W_mask).astype(float) on a 6″ box around the F140W centroid. Result:
- n = 2.39 (lower than typical ellipticals; nb07a's earlier value not directly compared yet)
- r_eff = 1.39″ (but this is the F140W subimage Sersic, not the full curve-of-growth R_e)
- ellip = 0.13
- x_0, y_0 → projected to (50.07, 55.29) in IFU pixels

The Sersic fit's centroid (50.07, 55.29) is 0.302″ from HST_mean. This is the center the deflector light profile prefers under a single-component Sersic model.

Note: Sersic n=2.39 < 4 (de Vaucouleurs) suggests either the deflector is not a pure de Vaucouleurs profile (possible — AGEL0206 may have disk component) OR the masked deflector core biased the n estimate. nb07a's full 2D Sersic on F140W and F200LP separately should be re-checked.

---

## 7. Open items

1. **σ_e-per-center comparison results** (running) — will determine if centering at the 0.3″ level matters
2. **PSF-aware mask comparison** — re-run §6cum with nb06's seeing-convolved mask vs nb07c's hard F200 reprojection. Quantify the seeing-wing arc-light contribution
3. **Decide on canonical center for the paper** — given that all 5 centers agree at < 0.4″ and σ_e is unlikely to shift by more than the SPS systematic, the choice is partly conventional. Recommend either:
   - Stick with HST_mean (matches nb07c, closest to the IFU peak after arc-masking)
   - Or switch to IFU_WL_peak / IFU_peak_arcmask for self-consistency with the kinematic data
4. **Re-examine the asymmetric bright IFU lobe** — is it really seeing-broadened arc, or is there a separate bright source / deflector substructure?

---

## Files

- `results/figures/center_check.png` — original 5-center diagnostic plot (full IFU FoV)
- `results/figures/center_candidates.png` — zoomed visualization with R_e circle
- `results/annular_bootstrap_07c_centers/` — per-center σ_e cache (in progress)
- `/tmp/center_compare.py` — comparison script
- `/tmp/center_compare.log` — running log
- This file — methodology investigation summary
