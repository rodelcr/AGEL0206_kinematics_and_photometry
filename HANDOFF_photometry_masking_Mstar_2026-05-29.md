# HANDOFF — Principled arc masking + M★ re-derivation — 2026-05-29

## TL;DR

Built and validated an **objective, reproducible arc-masking recipe** for the AGEL0206
deflector photometry, extended it to all four bands (incl. JWST), and **re-derived the stellar
mass** with a full systematic budget.

- **NEW headline M★ = log(M⋆/M☉) = 11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors
  [**11.04 ± 0.14 +0.32** at 20%]. One-sided +sys = deflector light under the arc (Sérsic fill-in),
  reaching **11.46 / 11.36**. **Supersedes the old expert-aperture 11.33** (smaller masks).
- **Explicit masking-approach systematic on M★ = ±0.16 dex** (peak-to-peak/2 across all masking
  approaches; symmetric statement of the one-sided +0.31).
- σ_e headline **unchanged** (269.62 ± 11.77 km/s) — independent of this work.
- R_e headline **unchanged** (2.305″) — confirmed; only a diagnostic script was buggy (now fixed).

## What we did, in order

1. **HST arc-mask verification** (`scripts/arc_mask_verification.py`, `notebooks/12` Part I).
   Two objective selectors reproduce the hand-painted F200LP mask: (i) color m_F200LP−m_F140W
   (arc bluer than the red deflector), (ii) Sérsic-residual (subtract a 2D model, flag positive
   excess >3σ). **Sérsic-residual at k=3 reproduces the expert F200LP photometry to 0.016 mag**,
   R_e to 3%. S/N-regime sweep (snr {2..20}, k {2..8}) shows k=3 is the clean saturation point
   on F200LP. Subtraction diagnostics confirm the arc is isolated (arc-free residual ≈ N(0,1)).

2. **4-band audit** (`scripts/principled_mask_photometry.py`). My aperture photometry reproduces
   the official headline mags under the expert mask to ≤0.007 mag (audit passes). **Independent
   per-band single-Sérsic masks OVER-MASK F140W/JWST** by 0.21–0.36 mag (single Sérsic under-fits
   the bright extended IR galaxy; median residual +1.0) — do NOT use them there.

3. **Mask-method comparison** (`scripts/mask_method_comparison.py`). **Reprojecting the F200LP
   footprint reproduces the expert JWST mags to 0.01–0.02 mag** = the transferable standard.
   Quantified the raw-vs-Sérsic-fill-in and analytic-total fluxes per band/mask.

4. **2-component model + IR extension** (`scripts/photometry_systematics.py`). Upgraded to a
   **2-component (bulge+disk) Sérsic** (cuts JWST galaxy-body RMS 3.3–3.5×). F200LP locates the
   arc; the deeper IR bands **extend** it (region-grow into contiguous 2-comp-residual source
   pixels). Compared **per-band vs global** (union) masks. Key result: **raw photometry is
   mask-size-dependent** (up to 0.78 mag), **the Sérsic fill-in is mask-definition-independent**
   (per-band vs global agree to 0.01 dex).

5. **M★ budget** — 8 Bagpipes fits {per-band, global} × {raw, filled} × {10%, 20%}
   (`results/photometry_systematics_Mstar.npz`, `Mstar_budget.png`). Adopted (user decision):
   **empirical raw-central, one-sided-up systematic, both errors quoted**. The IR-extended masks
   are large, so raw is biased low (discards under-arc light) and the fill-in marks the upper reach.

6. **Explicit masking systematic on M★** = ±0.16 dex (`results/Mstar_masking_systematic.npz`,
   `Mstar_masking_systematic.png`). Decomp: under-arc(raw↔filled) ±0.15, mask-def(per-band↔global)
   ±0.004, mask-extent(expert↔IR-extended) 0.18.

7. **PSF-convolved fill model** (`scripts/psf_fill_model.py`, `results/psf_fill_model.npz`).
   PSF-convolved 2-comp Sérsic (Gaussian at instrument FWHM; env lacks webbpsf) shifts the filled
   mag by **≤0.004 mag** all bands → ΔM★ ≪0.01 dex. **The fill leg is PSF-robust** (arc outside
   the PSF core). Firmed up, no longer a flagged unknown.

8. **R_e pixscale fix** (`scripts/measure_Re.py`). `hst_Re` hard-coded 0.08″/pix for both bands;
   the F200LP cutout is **0.05″** → diagnostic F200LP R_e was 1.6× too high. Fixed (reads WCS).
   F200LP diagnostic R_e 3.05→1.91″. **Headline R_e=2.305″ UNAFFECTED** — it comes from
   `final_sigma_e.py`, which already reads the WCS scale (F200LP=2.52″, F140W=2.16″, mean≈2.34″).

9. **ApJL paper figures updated** (`../AGEL_0206_ApJL_Figures/figures.ipynb`, cells 19/35/36;
   backup `figures.ipynb.bak.Mstar_principled_2026-05-29`). **Figure 2** SED inset → new posterior
   (median 11.16) + asymmetric title 11.16 +0.32/−0.08 + shaded fill-in reach to 11.46.
   **Figure 4** M•–M★ point → 10^11.16 with asymmetric xerr −0.08/+0.32. PDFs regenerated.

## Files

**New scripts:** `arc_mask_verification.py`, `principled_mask_photometry.py`,
`mask_method_comparison.py`, `photometry_systematics.py`, `psf_fill_model.py`.
**Edited scripts (backward-compatible):** `bagpipes_sersic_refit.py` (`run_bagpipes_for_mags(...,err_frac=)`),
`measure_Re.py` (`hst_Re(mask_override=)` + WCS pixscale).
**New notebook:** `12_arc_mask_verification.ipynb` (32 cells, Part I + Part II, executed).
**New caches:** `results/{arc_mask_verification, principled_mask_photometry, mask_method_comparison,
photometry_systematics, photometry_systematics_Mstar, mask_method_Mstar, arc_mask_bagpipes,
Mstar_masking_systematic, psf_fill_model}.npz`.
**New figures:** `results/figures/{arcmask_*, principled_*, maskcompare_*, photsys_*, Mstar_budget,
Mstar_masking_systematic}.png`.
**Docs:** `DRAFTING_FACTS_paper_2026-05-29.{md,tex,pdf}` (headline M★, §2.1.1b, §2.1.4, §3.2),
`TESTS_AND_DIAGNOSTICS.{md,tex,pdf}` (N1–N9, L4, A3b), `CLAUDE.md`,
`NOTES_arc_mask_verification_2026-05-29.md`, `NOTES_photometry_mask_systematics_2026-05-29.md`,
memory `project_arc_mask_verification.md`.
**iCloud (for phone access):** `~/Library/Mobile Documents/com~apple~CloudDocs/AGEL0206_Draft/`
(md, pdf, ApJL figure PDFs, key figures, self-contained `12_arc_mask_verification.html`).

## Outstanding TODOs

1. **σ_e masking systematic** (analogue of the ±0.16 dex M★ one): reproject each masking approach
   (F200LP-located, IR-extended per-band, global) to the IFU grid, re-run
   `scripts/run_wide_sigma_e.py --cube new_clean_hei` under each, quote the spread. NOT done.
2. **Curve-of-growth reconciliation:** `measure_Re.hst_Re` (1.91″) vs `final_sigma_e.curve_of_growth`
   (2.52″) for F200LP disagree by ~0.6″ even at the same pixscale — different annulus binning.
   Headline uses final_sigma_e; reconcile the two algorithms.
3. **Decision to revisit:** the headline M★ dropped 11.33→11.16 because the IR-extended masks remove
   more flux (raw biased low, +0.31 one-sided up). Confirm this empirical framing reads well in the
   draft (DRAFTING_FACTS §3.2, ApJL Fig 2/Fig 4) before it propagates into prose.

## Key decisions (user, 2026-05-29)

- Headline = **reprojected-HST masks** as the standard; **raw-central, fill-in as one-sided +sys**.
- Mask **locator = F200LP**; IR bands **extend** (not independent per-band Sérsic — it over-masks).
- Fill model = **2-component Sérsic** (PSF firmed up at ≤0.004 mag).
- Quote **both 10% and 20%** flux errors; quote the **explicit masking systematic** (±0.16 dex).
