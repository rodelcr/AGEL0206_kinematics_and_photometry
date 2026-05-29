# NOTES — Arc-light mask verification (F200LP/F140W) — 2026-05-29

## TL;DR

The hand-painted F200LP arc mask is **reproducible from objective criteria**, and the
photometry headline (mag, R_e, M★) is **invariant to the masking method**:

| quantity | no mask | expert (hand) | Sersic-residual (auto) | color (auto) |
|----------|---------|---------------|------------------------|--------------|
| F200LP AB mag | 22.192 | **22.613** | 22.629 (Δ **+0.016**) | 22.393 (Δ −0.220) |
| F200LP R_e ('') | 3.165 | 3.053 | 2.965 (Δ −0.09) | 3.038 (Δ −0.02) |
| F140W AB mag | 18.765 | **19.140** | 19.343 (Δ +0.21, over-masks) | 19.004 (Δ −0.13) |
| log M★ | — | **11.332** (−0.091/+0.073) | — | **11.366** (auto F200LP+F140W) Δ **+0.034** |

**Headline:** masking *matters* (F200LP no-mask is 0.42 mag = 47% brighter — real arc
contamination), but *which objective method* you use to build the mask does not: the
**Sersic-residual method reproduces the expert F200LP photometry to 0.016 mag** (≪ the
adopted 10% / 0.1-mag flux error), R_e to 3%, and M★ to +0.034 dex (≪ ±0.08 posterior).
→ The Sersic-residual method can be **promoted to the documented standard for F200LP**.

## What was done

New backing script `scripts/arc_mask_verification.py` (importable + `__main__` CLI), display
to be wired into `notebooks/12_arc_mask_verification.ipynb`. Two independent, objective
arc-pixel selectors, an S/N-regime sweep, and the full photometry-invariance chain.

### (A) Color (F200LP − F140W)
Reproject F140W surface brightness onto the F200LP grid (`scipy.ndimage.map_coordinates` + WCS
— `reproject` is **not** installed in ISMgas), build an AB color map in a **fixed blue−red
convention** (m_F200LP − m_F140W), so the blue z=1.302 arc is always *more negative* than the
red deflector core on both grids. Flag pixels bluer than the core by `DCOLOR_THRESH=0.5` mag,
with S/N > threshold and outside a 0.4″ core. Recovers up to ~73% of the expert mask on F200LP.

### (B) Sersic-residual
Fit a 2D Sersic to the deflector (local multi-init fit with **free position angle**, seeded
from the aperture-geometry ellip/θ in `example_outputs/*_params.json` — the shared
`sersic_total_photometry.fit_sersic2d` and `run_isource_shape_sweep` versions freeze θ=0 and
rail ellip→0), render on the full grid, subtract, flag positive excess > `k·σ_sky`. Diagnostic
plots (`results/figures/arcmask_{band}_subtraction.png`) confirm: model tracks the azimuthal
profile, residual reveals the arc as coherent positive excess, and the arc-free residual is
≈ N(0,1) (a positive tail = the arc).

### S/N-regime sweep (the "how much masking vs noise" question)
`SNR_REGIMES=[2,3,5,8,10,15,20]` (color), `KSIGMA_REGIMES=[2,3,5,8]` (residual). For each:
masked-px count, IoU vs expert, expert-recovery, **masked-flux fraction from the Sersic model**
(over-masking indicator), and Δmag. Key behaviour:
- **F200LP**: model-flux-fraction stays low (0.03–0.10) across the whole sweep → the automated
  masks remove arc light, not deflector light. Δmag from the Sersic mask saturates near the
  expert value at k≳3. This is the natural stopping point.
- **F140W**: model-flux-fraction is high (0.18–0.45) and rises as the threshold loosens → in
  the IR the arc sits on **bright deflector light**, so aggressive masking eats the galaxy.
  The "right" amount is ~k>5 (Δmag +0.05); k>3 over-masks (+0.21). This is the over-masking
  onset made quantitative (cf. nb07d: dilating the F200 mask inflates σ — noise-driven).

## Caveats / honest findings

1. **F140W Sersic-residual over-masks.** A single Sersic leaves a coherent negative ring in the
   F140W residual (model mismatch to the bright IR galaxy; visible in
   `arcmask_F140W_subtraction.png`). The **color method is more robust on F140W**; and the
   F140W arc is intrinsically faint (blue source, red IR deflector) and its expert mask was
   imported from F200LP anyway. Recommendation: Sersic-residual standard for F200LP, color
   cross-check for F140W.
2. **F200LP Sersic fits ellip≈0** (genuinely rounder in the blue / arc-filled); F140W fits
   ellip 0.10, θ≈165°. F200LP result is excellent regardless because the blue deflector is
   faint, so model imperfection barely shifts the aperture flux.
3. **`measure_Re.hst_Re` hardcodes pixscale=0.08″/pix for BOTH bands** (`HST_FILES['F200LP']`),
   but the F200LP cutout is actually **0.05″/pix**. ΔR_e *between masks* is unaffected (cancels),
   but the **absolute** F200LP R_e from `measure_Re.py` is biased — flag for a separate fix; may
   affect the headline R_e=2.305″ if that band entered the CoG mean at 0.08. **Do not silently
   fix here** — needs its own check against the F200LP CoG provenance.

## Files

- New: `scripts/arc_mask_verification.py`, `results/arc_mask_verification.npz`,
  `results/arc_mask_bagpipes.npz`, `results/figures/arcmask_*.png`.
- Edited (backward-compatible): `scripts/measure_Re.py` (`hst_Re` gains `mask_override=`),
  `scripts/bagpipes_sersic_refit.py` (factored `run_bagpipes_for_mags(...)`).
- Reused: `sersic_total_photometry.load_hst/find_center_2dg`, `measure_Re.hst_Re`.

## Outstanding (later TODO)

- **Spectroscopic invariance** (explicit later TODO): reproject the automated F200LP mask to the
  IFU grid and re-run `scripts/run_wide_sigma_e.py --cube new_clean_hei` to confirm σ_e is
  unchanged. NOT done in this arc.
- Fix the `measure_Re` F200LP pixscale (item 3) under its own provenance check.
- Notebook `12_arc_mask_verification.ipynb` display driver + TESTS_AND_DIAGNOSTICS rows.
