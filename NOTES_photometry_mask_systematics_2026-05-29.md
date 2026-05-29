# NOTES — Photometry arc-mask systematics (4-band, principled) — 2026-05-29

Successor to `NOTES_arc_mask_verification_2026-05-29.md` (HST-only verification). This
arc extends the principled arc-masking to **all four photometry bands** (F200LP, F140W,
F150W2, F322W2), quantifies the **deflector-light-under-the-arc** systematic via Sersic
fill-in, and propagates everything to a **M★ budget (stat ⊕ sys)**.

## Decisions taken (user, 2026-05-29)

1. **Headline basis** = reprojected-HST masks + **raw** photutils photometry; quote the
   Sersic **fill-in as a systematic** (not the central value).
2. **Fill model** = upgrade single-Sersic → **2-component (bulge+disk) Sersic** (PSF
   convolution unavailable in ISMgas — `webbpsf` not installed, no PSF FITS).
3. **Mask locator** = **F200LP** (most accurate for the blue z=1.302 source light);
   the deeper/sharper **IR bands reveal the source's further extent** → grow the F200LP
   footprint into IR-significant source pixels (region-growing).
4. Compare **per-band masks vs a global union mask** as a systematic.
5. Flux errors: **10% and 20%** (20% to be conservative about flux removed by masking).
6. Quote **statistical (posterior) ⊕ systematic** on the final M★.

## Method

- **Mask locator**: F200LP hand mask (2512 px, the master arc footprint) reprojected to
  every band via WCS + `scipy.ndimage.map_coordinates` (no `reproject` pkg in env).
- **NOT the union of F200LP∪F140W**: verified F140W's mask covers the deflector CORE
  (17 px at r<0.4″, 50 px at r<0.8″ when reprojected to F200LP) — cf.
  `feedback_masking_strategy`. Using F200LP-only avoids core contamination.
- **Deflector model**: 2-component Sersic2D (bulge n∈[2,6] + disk n∈[0.5,2], shared
  center), masked pixels excluded from the fit ("photometry re-informs the Sersic").
  **Validated**: on JWST the single Sersic had median fractional residual **+1.000**
  (model ≈ 0 vs the galaxy — total under-fit, the cause of the catastrophic over-masking),
  while 2-comp cuts galaxy-body RMS **3.3–3.5×** (F150W2 3.25×, F322W2 3.46×) and the
  pedestal to −0.25. F140W single was already adequate (1.05×).
- **IR source extension** (`grow_extension`): threshold the 2-comp residual at
  **k=3·σ_sky** (positive excess only), region-grow from the F200LP arc seed
  (`scipy.ndimage.label`, keep components touching the seed), within 0.4″<r<5.0″. This
  CANNOT re-grab the symmetric core pedestal or unrelated field companions (not contiguous
  with the arc).
- **per-band mask** = F200LP-reproj ∪ (band's IR extension).
- **global mask** = union of ALL per-band masks reprojected onto each band.
- **Photometry** (`flux_mag`): elliptical aperture from `example_outputs/*_params.json`,
  canonical ZPs. **raw** = photutils masked sum (under-counts under-arc light); **filled**
  = masked pixels replaced by the 2-comp model, then summed (recovers under-arc light).

## Findings (photometry side)

### Audit / correctness
My aperture photometry reproduces the official headline AB mags under the expert mask to
Δ ≤ 0.007 (F200LP/F150W2/F322W2 = 0.000, F140W +0.007) — deltas below are trustworthy.

### Independent per-band Sersic-residual mask OVER-MASKS on F140W/JWST
(from `scripts/principled_mask_photometry.py`, single-Sersic) — removes 0.21–0.36 mag too
much; mdlfrac low (0.02–0.03) was a FALSE reassurance (low because the under-fit model is
faint, not because masked pixels are clean). **Do not use independent per-band Sersic as
the headline** — only clean on F200LP. This motivated the 2-comp upgrade + IR-extension.

### IR reveals real extra source extent (per-band masks)
| band | F200-reproj px | +IR-extension px | per-band px |
|------|---------------|------------------|-------------|
| F200LP | 2512 | +272 | 2784 |
| F140W | 988 | +1121 | 2109 |
| F150W2 | 6649 | +9119 | 15768 |
| F322W2 | 1582 | +2676 | 4258 |
The arc is much more extended in the deep IR (F150W2 mask ~6× the HST footprint), as
expected — JWST sees source flux HST misses.

### raw is mask-sensitive; FILLED is mask-robust (key result)
| band | per-band raw | global raw | per-band filled | global filled |
|------|--------------|-----------|-----------------|---------------|
| F200LP | 22.7009 | 23.4761 | 22.4578 | 22.5175 |
| F140W | 19.4329 | 19.7285 | 18.8948 | 18.9024 |
| F150W2 | 19.2747 | 19.5145 | 18.6851 | 18.6828 |
| F322W2 | 18.9862 | 19.1452 | 18.3142 | 18.3134 |
- **raw** swings up to **0.78 mag** between per-band and global (global imposes the big IR
  extent onto the blue F200LP where the source is compact → over-masks). raw also depends
  on the IR extension depth.
- **filled** agrees per-band vs global to **≤0.06 mag** (F150W2/F322W2 to 0.001) — the
  model refills whatever is masked, so it is nearly mask-definition-independent.
- **Sersic fill-in correction (filled − raw)** = −0.24 to −0.96 mag depending on mask size:
  the deflector light hidden under the arc that raw photometry discards.

**Interpretation:** the filled photometry is the robust quantity; raw's sensitivity to the
mask definition is itself a systematic. (Earlier decision was raw-as-headline; the data now
argue the fill-in/mask-definition systematic on raw is sizable — flag for the M★ budget.)

## M★ budget (results/photometry_systematics_Mstar.npz)

8 Bagpipes fits: {perband,global}×{raw,filled}×{10%,20%}. nb02 priors. Headline (expert,raw,10%)=11.332.

| vector | 10% | 20% |
|--------|-----|-----|
| perband_raw    | 11.155 (±0.077) | 11.036 (±0.138) |
| perband_filled | 11.463 (±0.077) | 11.359 (±0.127) |
| global_raw     | 11.152 (±0.076) | 11.070 (±0.133) |
| global_filled  | 11.471 (±0.078) | 11.361 (±0.127) |

**KEY RESULT — the picture inverted vs the HST-only analysis:**
- The IR-extended masks are large (F150W2 +9119px), so **raw photometry is biased LOW and
  mask-size-dependent** (11.16) — it discards the deflector light under the now-large arc mask.
- **filled is robust**: per-band 11.463 vs global 11.471 (mask-def-independent to 0.01 dex),
  and matches the HST-only filled value (11.48). The 2-comp model refills whatever is masked.
- → The coherent method is **IR-extended mask + Sersic FILL** (exclude source-contaminated
  pixels AND restore the deflector light under them). raw+large-mask is incoherent (throws
  away deflector flux). This **inverts the earlier raw-central/fill-systematic decision**
  (made when HST masks were small and raw≈filled).
- **fill-in shift** (filled−raw) ≈ +0.31 dex (huge, because masks are big).
- **per-band vs global**: negligible for filled (0.01 dex), ~0.003 for raw.
- **10% vs 20% error**: shifts the CENTRAL by ~0.10 dex (looser errors → prior pulls mass
  down) AND widens stat (±0.08→±0.13). Not just a width effect.

**ADOPTED HEADLINE (user decision 2026-05-29): empirical raw-central, one-sided-UP systematic, quote both errors.**
  log M★ = **11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors [**11.04 ± 0.14 +0.32** at 20%].
  - central = perband_raw (the empirical masked measurement with the principled IR-extended mask).
  - stat = Bagpipes posterior half-width at the adopted flux error.
  - **systematic is one-sided UP** = the Sérsic fill-in (deflector light under the arc that raw
    discards) → reaches log M★ = 11.46 (10%) / 11.36 (20%). mask-definition term ±0.003–0.034 (folded in).
  - Rationale (user): "stay empirical and show the uneven systematic that indicates it could well be
    higher." Raw masked photometry is effectively a lower bound; the fill-in marks the upper reach.
  - Supersedes the older expert-aperture 11.33 (smaller masks; sits mid-bracket).
  - Figure: `results/figures/Mstar_budget.png`. Propagated to DRAFTING_FACTS / CLAUDE.md / memory.

## Explicit masking-approach systematic on M★ (2026-05-29, user-requested)

**±0.16 dex** (10%: ±0.160, 20%: ±0.162) — peak-to-peak/2 of the M★ medians across ALL masking
approaches {expert-aperture 11.33; per-band/global × raw/filled spanning 11.15–11.47}.
`results/Mstar_masking_systematic.npz`, fig `Mstar_masking_systematic.png`. Decomposition (10%):
- **under-arc light (raw↔filled): ±0.15** — dominant; = the adopted one-sided +0.31.
- **mask-definition (per-band↔global): ±0.004** — negligible (filled is mask-def-independent).
- **mask-extent (expert↔IR-extended-raw): 0.18** — the IR-extension effect on raw.
±0.16 is the symmetric statement of the same effect the headline quotes one-sided (+0.31).

## R_e pixscale fix (2026-05-29) — measure_Re.hst_Re

`measure_Re.hst_Re` hard-coded 0.08″/pix for BOTH bands; the F200LP cutout is **0.05″/pix**, so the
diagnostic F200LP R_e was biased high by 0.08/0.05 = 1.6×. **Fixed**: `hst_Re` now reads the pixel
scale from the WCS (`proj_plane_pixel_scales`); HST_FILES F200LP literal corrected 0.08→0.05.
- F200LP proper-mask R_e: 3.05″ → **1.91″** (corrected); F140W unchanged (2.21″; 0.08 was right).
- **Headline R_e = 2.305″ CONFIRMED UNAFFECTED**: it comes from `scripts/final_sigma_e.py`
  (separate path that already reads pix scale from WCS — F200LP=2.52″ at 0.05″, F140W=2.16″,
  mean ≈ 2.34″ ≈ 2.305″). Only the measure_Re diagnostic npz was wrong.
- **Flag (pre-existing, NOT the pixscale bug):** the two curve-of-growth implementations disagree
  on F200LP even at 0.05″ — measure_Re 1.91″ (1-pix annuli to 80px) vs final_sigma_e 2.52″
  (r_step=0.08″ to 6″). Headline uses final_sigma_e. Reconciling the two CoG algorithms is a
  separate TODO.

## TODO — masking systematic on σ_e (analogue of the M★ one)

Reproject each masking approach (F200LP-located, IR-extended per-band, global) to the IFU grid and
re-run `scripts/run_wide_sigma_e.py --cube new_clean_hei` under each → quote an explicit
masking-approach systematic on σ_e, same as the ±0.16 dex M★ one. NOT done in this arc.

## Files

- Scripts: `scripts/principled_mask_photometry.py` (single-Sersic audit, 4-band),
  `scripts/mask_method_comparison.py` (expert/HST-reproj/Sersic × raw/filled/total),
  `scripts/photometry_systematics.py` (2-comp, per-band+IR-extension vs global, raw/filled).
- Edits: `bagpipes_sersic_refit.run_bagpipes_for_mags(..., err_frac=)`,
  `measure_Re.hst_Re(mask_override=)`.
- Caches: `results/{principled_mask_photometry,mask_method_comparison,photometry_systematics}.npz`,
  `results/{arc_mask_bagpipes,mask_method_Mstar,photometry_systematics_Mstar}.npz`.
- Figures: `results/figures/{principled_,maskcompare_,photsys_}*.png`.

## Caveats

- **PSF effect on the fill — QUANTIFIED (2026-05-29), no longer just flagged.**
  `scripts/psf_fill_model.py` forward-models a PSF-convolved 2-component Sersic (Gaussian PSF at
  the per-band instrument FWHM: F200LP 0.08″, F140W 0.14″, F150W2 0.06″, F322W2 0.13″; env lacks
  webbpsf so parametric) and recomputes the filled mag. **ΔPSF(filled) ≤ 0.004 mag in every band**
  (max 0.004, F200LP) → ΔM★ ≪ 0.01 dex, negligible vs stat ±0.08 and masking ±0.16. The fill
  region is the arc at r~1–4″, outside the PSF core, so the fill is PSF-robust.
  `results/psf_fill_model.npz`. (Empirical-PSF refinement possible but unnecessary given ≤0.004 mag.)
- **(historical)** 2-comp Sersic without PSF was used for the headline fill; the line above shows that
  light + analytic totals carry PSF-related uncertainty; the arc-region fill (what we use)
  is less PSF-sensitive.
- `measure_Re.hst_Re` hardcodes 0.08″/pix for F200LP (actual 0.05″) — ΔR_e between masks
  cancels; absolute F200LP R_e from that script biased. Flagged, NOT fixed.
- 2-comp `medfrac ≈ −0.25` on JWST = slight over-subtraction in the body → the IR extension
  is CONSERVATIVE (won't flag the over-subtracted body as source) — good for over-masking
  safety, may slightly under-grow the true extent.
