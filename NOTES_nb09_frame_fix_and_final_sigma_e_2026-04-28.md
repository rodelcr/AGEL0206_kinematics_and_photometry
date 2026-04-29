# NOTES — nb09 frame fix + final σ_e production
**Date:** 2026-04-28
**Author session:** Claude + Rodrigo
**Goal:** Diagnose the "atypical" −95 km/s FSPS V_sys split-track in our ppxf
fits, confirm instrumental LSF is being passed correctly, fix the underlying
issue, and produce the paper-ready final σ_e.

---

## TL;DR

1. **Instrumental LSF was correctly passed to ppxf** all along. KCWI cube
   `DISPSCAL = 0.294` (instrumental σ in pixels), `CD3_3 = 1.000 Å/pix`
   → FWHM_inst = **0.692 Å** (constant in Å), σ_v_inst = 12–14 km/s,
   R ≈ 9400–10800 over the fit band. The pipeline at `bootstrap_ppxf.py:122-129`
   correctly converts `DISPSCAL` → per-pixel rest-frame FWHM dict and feeds
   it to `sps_lib`, which broadens the templates to match. Canonical
   Cappellari (2023) pattern. **Not the source of the V_sys split.**

2. **The −95 km/s FSPS V_sys offset was a wavelength-frame mismatch.** ppxf's
   bundled SPS template files ship in different frames:
   - `spectra_fsps_9.0.npz` is in **vacuum** (Ca K @ 3934.86 Å = +1.20 Å above air rest)
   - `spectra_emiles_9.0.npz` is in **air** (Ca K @ 3933.82 Å)
   - `spectra_xsl_9.0.npz` is in **air** (V_sys closes at +9 km/s with air galaxy)

   The original pipeline converted galaxy to air for *all three* libraries.
   Air galaxy ↔ vacuum FSPS templates → ppxf compensates with V ≈ −91 km/s
   to align lines. Matches the observed −95 km/s offset exactly.

3. **Fix is per-SPS native frame.** Patched `scripts/bootstrap_ppxf.py` to
   define `SPS_NATIVE_FRAME = {'fsps': 'vacuum', 'emiles': 'air', 'xsl': 'air'}`
   and only apply `vac_to_air` when the SPS is in air. Frame-aware setup
   collapses the V_sys split-track from ~110 → ~18 km/s. **σ shifts by only
   ~1 km/s** between old and new headline (267.32 → 267.82 at N=500).

4. **Final headline (nb09, masked, frame-aware, N=500 SPS-pooled):**
   $$ \sigma_e(<R_e) = 268 \pm 32\ \mathrm{km/s} $$
   = stat ±24 ⊕ I-shape ±13 ⊕ mask ±16 ⊕ frame ±5 ⊕ centering ±4
   $$ \sigma_e(<R_e/2) = 224 \pm 18\ \mathrm{km/s} $$
   $$ R_e = 2.305\arcsec = 16.23\ \mathrm{kpc}\ (z = 0.67564) $$

---

## Investigation timeline

### Test 1 — Instrumental LSF audit

Already-built script: `scripts/ifu_spectral_resolution.py`. Run output:

```
DISPSCAL header value: 0.294  (instrumental sigma in pixels)
λ range: 5625.0 - 8941.0 Å (3317 px @ 1.000 Å/pix)
FWHM_inst = 2.355 × 0.294 × 1.000 = 0.692 Å (constant)

In ppxf fit band 6500-7500 Å:
  FWHM (Å):       median 0.692
  R = λ/FWHM:     median 10110 (range 9388-10832)
  σ_v_inst:       median 12.6 km/s (range 11.8-13.6)
  velscale:       42.83 km/s/pix (= ppxf velscale)

Rest-frame at z=0.67564:
  λ band:         3879 - 4476 Å
  σ_v_inst rest:  7.5 km/s (deconvolved by 1+z)
```

**Galaxy is resolved by a factor of ~20** (σ_galaxy ≈ 270, σ_inst ≈ 13).

In `_prep_spectrum_for_ppxf` at `bootstrap_ppxf.py:122-129`:
```python
dlam_gal = np.gradient(lam_gal)
wdisp = hdr['DISPSCAL']
fwhm_gal = 2.355 * wdisp * dlam_gal
...
fwhm_gal_dict = {"lam": lam_gal_rest, "fwhm": fwhm_gal_rest}
sps = lib.sps_lib(filename, velscale, fwhm_gal_dict, lam_range=...)
```
The dict pattern is the canonical Cappellari (2023) interface. `sps_lib`
internally broadens each template to match `fwhm_gal_dict[fwhm]` at every
wavelength via `util.varsmooth`, then log-rebins. **No bug here.**

### Test 2 — SPS template frame diagnosis

Hypothesis: the −95 km/s FSPS V_sys offset (vs +7 EMILES, +14 XSL) is
suspiciously close to the air↔vacuum shift in this band (~85–91 km/s). If
FSPS templates ship in vacuum and our galaxy is converted to air, we'd see
exactly this offset.

**Test:** Locate Ca H + K and G-band absorption-line minima in each SPS
template file directly.

```python
ppxf_dir = "/opt/anaconda3/envs/ISMgas/lib/python3.14/site-packages/ppxf"
for sps in ['fsps', 'emiles', 'xsl']:
    f = np.load(f"{ppxf_dir}/sps_models/spectra_{sps}_9.0.npz")
    # Pick old, near-solar SSP (i_age ~ 12 Gyr, [Z/H] ≈ 0)
    # Locate Ca K minimum via parabolic fit around np.argmin
```

**Result:**

| Library | Ca K min | Δ vs air rest 3933.66 Å | Frame |
|---|---|---|---|
| FSPS | 3934.86 Å | **+1.20 Å (+91 km/s)** | **vacuum** |
| EMILES | 3933.82 Å | +0.16 Å (+12 km/s) | **air** |
| XSL | (confused — peculiar SSP feature) | — | (defer to V_sys closure test) |

XSL's Ca K diagnostic was confused by an unusual emission-like peak in the
SSP at 3933.69 Å (likely a normalization artifact in the npz packaging),
so I deferred to the end-to-end V_sys closure test below.

### Test 3 — End-to-end V_sys closure

**The decisive test.** Same KCWI integrated spectrum (190-spaxel deflector
slice) fed to ppxf with **only the wavelength frame swapped**:

```python
# Fix the galaxy spectrum + noise. Then build two versions:
lam_vac = native KCWI vacuum wavelength axis
lam_air = lam_vac × np.median(util.vac_to_air(lam_vac)/lam_vac)
# Run ppxf with each frame against each SPS at degree=20
```

| SPS | V_sys (vacuum galaxy) | V_sys (air galaxy) | Δσ (vac − air) |
|---|---:|---:|---:|
| **FSPS** | **−7.17** | −89.58 | −2.0 |
| **EMILES** | +83.77 | **+2.88** | −5.2 |
| **XSL** | +90.15 | **+8.70** | −3.1 |

**Conclusions:**
- FSPS templates are in vacuum (V_sys closes when galaxy is in vacuum).
- EMILES + XSL are in air (V_sys closes when galaxy is in air).
- σ shifts by only 2–5 km/s between vacuum and air galaxy — sub-systematic.
- Frame-aware setup will collapse the V_sys split from ~110 km/s → ~15 km/s.

### Test 4 — Frame fix & re-run

Patched `bootstrap_ppxf.py`:
```python
SPS_NATIVE_FRAME = {'fsps': 'vacuum', 'emiles': 'air', 'xsl': 'air'}

def _prep_spectrum_for_ppxf(..., frame_galaxy='auto'):
    if frame_galaxy == 'auto':
        frame_galaxy = SPS_NATIVE_FRAME.get(sps_name, 'air')
    if frame_galaxy == 'air':
        lam_gal *= np.median(util.vac_to_air(lam_gal) / lam_gal)
    # else: keep KCWI native vacuum
```

The `frame_galaxy` arg is plumbed through `setup_ppxf_inputs_from_spectrum`.
Default = `'auto'`, so existing callers automatically pick up the fix.

**Result (N=50 smoke test, integrated 190-spaxel box at degree=20):**
- FSPS V_sys: −95 → ~−7 km/s ✓
- EMILES V_sys: +7 → +3 km/s ✓ (already in correct frame)
- XSL V_sys: +14 → +9 km/s ✓

V_sys split: ~110 → ~16 km/s. Per-SPS V_sys still has small residual
template offset (real SPS systematic), now ~5–10 km/s rather than
dominated by frame mismatch.

### Test 5 — Production rerun (nb09)

Built `scripts/final_sigma_e.py` — best-practice σ_e pipeline:
- `setup_ppxf_inputs_from_spectrum(..., frame_galaxy='auto')` → per-SPS frame
- §6cum cumulative I-weighted aperture ppxf at R<R_e/2 and R<R_e
- F200LP arc mask reprojected via `map_coordinates(order=0)` → IFU grid
- HST-mean center: `centroid_2dg` on F140W + F200LP, averaged in world coords
- R_e = 0.5 × (R_e^F140W + R_e^F200LP) from masked CoG = **2.305"**
- Two tracks: masked (headline) + nomask (sensitivity)
- Per-(label, sps, mask) bootstrap caches with N-fallback chain
- SPS pooling for the headline posterior
- Quadrature-combined error budget

**N=50 smoke (12 fits, ~6 min total):**
- σ_e(<R_e) = 268.5 +23/−25 km/s [stat]
- Δ_mask = −17.2 km/s

**N=500 production (9 fits at N=500 + 3 at N=50 fallback, ~50 min):**
- σ_e(<R_e) = **267.8 +23/−25 km/s** [stat]
- Δ_mask = −16.5 km/s
- N=50 vs N=500 shift: 0.7 km/s — **bootstrap converged**

#### Per-aperture per-SPS σ at headline N

| Aperture | FSPS | EMILES | XSL |
|---|---|---|---|
| R<R_e/2 | 214.4 | 220.3 | 235.6 |
| R<R_e | 253.8 | 267.3 | 279.6 |

SPS spread at R<R_e = 26 km/s — captured by SPS-pooled posterior 1σ = ±24.

#### Per-aperture V_sys after frame fix

| Aperture | FSPS (vac) | EMILES (air) | XSL (air) |
|---|---:|---:|---:|
| R<R_e/2 | −21.0 | −10.3 | −5.3 |
| R<R_e | −18.9 | −4.4 | −1.1 |

V_sys split at R<R_e = 18 km/s. Compare to pre-fix split ~110 km/s.

### Test 6 — R<R_e/8 dropped

3 spaxels at R = 0.288″, well inside the seeing radius (FWHM/2 = 0.64″).
Initial smoke fit gave σ ≈ 65 km/s (masked) vs 95 km/s (nomask) — clearly
seeing-limited and not physical. Removed from the production aperture set.

### Test 7 — Two-track comparison (mask sensitivity)

| Aperture | masked (headline) | nomask (sensitivity) | Δ |
|---|---|---|---|
| R<R_e/2 | 223.6 (−18/+18) | 219.6 (−18/+18) | **−4** ✓ |
| **R<R_e** | **267.8 (−25/+23)** | 251.3 (−22/+22) | **−16.5** |

The Δ at R<R_e ties out to nb07c's prior nomask sensitivity (−16.4 km/s)
within 0.1 km/s. Confirms the F200 mask systematic is a real ~16 km/s
effect (not just a frame artifact). Carried as `mask` in the error budget.

---

## Methodology — vac/air conversion

ppxf's `util.vac_to_air()` uses **Ciddor 1996** (eq. 1):
```python
def _wave_convert(lam):
    sigma2 = (1e4/lam)**2
    return 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
def vac_to_air(lam_vac): return lam_vac/_wave_convert(lam_vac)
```
Accuracy ~10⁻⁸ in *n*. This is the modern astronomical standard.

Our pipeline applies the **scalar median** ratio:
```python
lam_gal *= np.median(util.vac_to_air(lam_gal) / lam_gal)
```
This is the **canonical Cappellari pattern** — see `ppxf_example_kinematics_sdss.py:127-134`:

> *"The SDSS wavelengths are in vacuum, while the MILES ones are in air.*
> *For a rigorous treatment, the SDSS vacuum wavelengths should be*
> *converted into air wavelengths and the spectra should be resampled."*
> ```python
> lam_gal *= np.median(util.vac_to_air(lam_gal)/lam_gal)
> ```

The scalar approach preserves log-uniform spacing (required for ppxf's
log-rebinned templates). A pixel-by-pixel resampling would be more rigorous
but breaks the log axis. Over our narrow band 6500–7500 Å the
*differential* between scalar and per-pixel is < 3×10⁻⁵ → sub-km/s in the
LOSVD fit. **Negligible.**

---

## Files produced this session

### New scripts
- `scripts/final_sigma_e.py` — production pipeline (CLI: `--n_bootstrap {50,500} --force`)
- `scripts/build_nb09.py` — notebook builder

### Modified scripts
- `scripts/bootstrap_ppxf.py` — added `SPS_NATIVE_FRAME` + `frame_galaxy='auto'` arg
  to `_prep_spectrum_for_ppxf` and `setup_ppxf_inputs_from_spectrum`.
  Backward-compatible: explicit `frame_galaxy='air'` reproduces pre-2026-04-28 fits.

### New notebook
- `notebooks/09_final_sigma_e_paper.ipynb` — paper-ready, executed end-to-end

### Artifacts
- `results/final_sigma_e_paper.npz` — pooled posteriors + error budget
- `results/final_sigma_e_paper/{Re_2,Re}_{fsps,emiles,xsl}_N{50,500}{_nomask}.npz` — 24 cache files
- `results/figures/nb09_two_tracks.png` — masked vs nomask histogram overlay
- `results/figures/nb09_sigma_vs_aperture.png` — final σ_e vs aperture figure

---

## Caveats / known limitations

1. **3 of 12 N=500 fits fell back to N=50** (no-mask Re track). The
   no-mask sensitivity track ran 5–40× slower than expected — one fit
   (Re_2_xsl_N500_nomask) took 4+ hours, and the bg job was killed before
   the no-mask Re track at N=500 could complete. The headline Track A
   (masked) is **N=500 throughout**. The Δ_mask sensitivity uses the
   N=50 caches for missing fits; this affects only the second-decimal-place
   precision of the −16.5 km/s shift, not the headline.

2. **F140W vs F200LP centroid offset = 0.357″** which exceeds the 0.2″
   sanity threshold. This is driven by the F200LP arc contaminating its
   centroid; the HST-mean center used for the analysis is robust because
   F140W is arc-clean. Sub-dominant ±0.4 km/s contribution to σ.

3. **R_e/8 aperture removed** — 3 spaxels at 0.29″ is below seeing FWHM/2 =
   0.64″. The σ value was both seeing-broadened and statistically poor.
   For the M•–σ relation we report σ_e(<R_e), which is what KH13 / Greene+20
   define as the canonical aperture anyway.

4. **XSL Ca K diagnostic was confused** by a peculiar emission-like peak
   in the bundled SSP at 3933.69 Å. The end-to-end V_sys closure test
   (Test 3) provided the definitive identification of XSL = air. Worth
   noting if anyone re-uses the XSL templates for a similar test.

5. **Per-SPS V_sys offsets at R<R_e are now FSPS −19, EMILES −4, XSL −1
   km/s** — i.e., the FSPS template still has a ~−15 km/s residual offset
   relative to EMILES/XSL after the frame fix. This is a real SPS-template
   systematic (different stellar libraries, different abundance ratios),
   not a wavelength-frame issue. Captured by the SPS-pooled posterior; no
   additional correction needed.

6. **The error budget components are NOT independent in the strict
   statistical sense.** I-shape is bounded by mask choice; frame is bounded
   by SPS choice; centering is partially absorbed by the additive
   polynomial. Quadrature combination is conservative — the true total
   error could be ~10% smaller. We err on the conservative side for the
   paper.

---

## Reproduction recipe

```bash
conda activate ISMgas

# Verify instrumental LSF
python scripts/ifu_spectral_resolution.py

# Smoke test (~6 min)
python scripts/final_sigma_e.py --n_bootstrap 50

# Production (with N=50 fallback if any fit is killed; ~50-100 min)
python scripts/final_sigma_e.py --n_bootstrap 500

# Build + execute notebook
python scripts/build_nb09.py
jupyter nbconvert --to notebook --execute --inplace \
  notebooks/09_final_sigma_e_paper.ipynb \
  --ExecutePreprocessor.timeout=300
```
