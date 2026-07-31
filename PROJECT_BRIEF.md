# PROJECT BRIEF — AGEL0206 deflector kinematics & photometry

**Audience.** Designer / agent producing a conference poster or talk
slides. This brief is self-contained — you don't need to read the
repo. For the technical catalogue see `TESTS_AND_DIAGNOSTICS.pdf`.

**Last updated.** 2026-06-14

---

## 1. One-line headline

We measure the central stellar velocity dispersion **σ_e = 267 ± 12 km/s**
inside one effective radius and the stellar mass **M★ = 3.2 × 10¹¹ M☉**
of the foreground deflector galaxy in the strong-lens system
**AGEL J020613-011417** at z = 0.676, anchoring it on the M•–σ
relation.

## 2. Three-bullet story

- **System.** AGEL0206 (DES J0206−01) is a galaxy-scale strong-lens
  system from the AGEL DR2 survey. The lensed source is a spiral
  galaxy at z = 1.302; the foreground deflector is a passive
  elliptical at z_def = 0.6756 (line-fit redshift; Δz = +95 km/s vs
  the AGEL DR2 catalogue value).
- **Measurement.** Keck/KCWI integral-field spectroscopy provides
  a 2D stellar-kinematic map. We pool three SPS template libraries
  (FSPS, EMILES, XSL) and run frame-aware ppxf with an arc-mask to
  extract σ_e inside R<R_e via cumulative-aperture I-weighted fitting
  — the literature-standard definition (Cappellari+ 2006 eq. 1, used
  in KH13, Greene+ 2020, SAURON/ATLAS3D/MaNGA).
- **Headline result.** σ_e(<R_e) = 267.31 ± 11.77 km/s; combined with
  M★ ≈ 3.0 × 10¹¹ M☉ (log M★/M☉ = 11.46; Bagpipes SED on HST WFC3 + JWST NIRCam
  photometry, validated-fit apcorr + quiescent prior) this places AGEL0206 cleanly on the M•–σ and M•–M★
  scaling relations and provides an independent strong-lens-side
  σ-aperture constraint for cosmographic uses.

## 3. Headline numbers (paper-ready)

```
σ_e(<R_e)            = 267.31 ± 11.77 km/s  (paper headline; asym −11.98/+11.58)
σ_e(<R_e/2)          = 224 ± 18 km/s     (gradient diagnostic)

R_e                  = 2.097″ = 15.26 kpc   (HST F140W+F200LP best-mask CoG)
log(M★/M☉)           = 11.5 ± 0.1 (stat) ± 0.2 (sys)  [2dp 11.46]  (validated-fit apcorr + quiescent SED prior, 2 R_e, 10%; was 11.50)
log(M★/M☉) [20% err] = 11.04 ± 0.14       (fill-in reach 11.46; supersedes 11.33)

z_deflector          = 0.67564          (line fit — Ca H+K, [O II], etc.)
z_source             = 1.302
```

### σ_e error budget

The ±32 km/s is a quadrature sum of five independent, propagated
systematics:

| Component | km/s | Origin |
|---|---|---|
| Statistical (bootstrap pooled 1σ) | ±24 | wild-bootstrap, N=500, 3-SPS pool |
| I(r) shape (10-shape sweep) | ±13 | F140W/F200LP/IFU/Sersic2D weight choices |
| Mask on/off Δ | ±16 | F200LP arc-mask vs no-mask |
| Frame fix (vacuum/air per SPS) | ±5 | SPS-template native frame correction |
| Centering (5-position sweep) | ±4 | HST F140W/F200LP centroid + IFU peak |
| **Quadrature total** | **±32** | |

## 4. Figure inventory (poster-ready, in repo)

All paths relative to repo root.

| File | Use it for | What it shows |
|---|---|---|
| `results/figures/nb09_sigma_vs_aperture.png` | **Hero figure 1.** | σ_e increases with aperture from R<R_e/2 → R<R_e. Masked headline (blue) vs no-mask sensitivity (red). ±total-budget band at R_e. |
| `results/figures/nb09_sed_and_Mstar.png` | **Hero figure 2.** | Two-panel: observed-vs-modelled SED at four filter pivots + log M★ posteriors (aperture vs Sersic-total). |
| `results/figures/nb09_crosschecks.png` | Method-validation panel. | Three independent estimators (§6cum / §7 / §7b) agreeing within ±total budget. |
| `results/figures/nb09_re_source_systematic.png` | R_e-choice systematic. | σ_e(<R_e) for four R_e definitions (mean / F140W / F200LP / Ca H+K). |
| `results/figures/nb09_mask_weight_sweep.png` | Mask-sensitivity supp. fig. | σ_e(w) curve for w∈{0, 0.25, 0.5, 0.75, 1}. Threshold-dominated, super-linear. |
| `results/figures/nb09_two_tracks.png` | Tracks A vs B comparison. | Posterior distributions for masked / no-mask at R<R_e. |
| `results/figures/nb09_sigma_vs_degree.png` | Polynomial-saturation diagnostic. | σ vs degree 15-29; flat → polynomial saturated. |
| `results/figures/nb07f_recovery.png` | Dilution-mechanism panel. | Track A/B/D bar chart with 145% recovery fraction. |
| `results/figures/nb07f_posteriors.png` | Dilution-mechanism histograms. | σ_e posteriors for masked, no-mask, no-mask+arc-sky. |
| `figures/Re_comparison.png` | Effective-radius slide. | F140W/F200LP/IFU CoG curves with R_e markers. |
| `figures/ifu_spectral_resolution.png` | Instrument-LSF slide. | KCWI σ_inst ≈ 12-14 km/s (galaxy σ resolved by ~20×). |
| `results/AGEL0206_spectra_SED_fit.pdf` | Full Bagpipes corner+SED. | Bagpipes-native posterior+SED plot. |

## 5. Methodology in plain language

Fits stellar absorption lines (Ca H+K, G-band, rest 3879-4476 Å) in a
single I-weighted aperture spectrum to recover the line-of-sight
velocity dispersion σ. Three things that make this measurement
non-trivial here and how we handle each:

1. **The lensed arc contaminates the spectrum.** The arc is at z=1.302
   and contributes featureless UV continuum at our deflector rest-frame
   fit window. We hard-mask the arc using a F200LP-tuned segmentation;
   sensitivity to that mask choice (±16 km/s) is propagated. We
   independently verified that the no-mask shift is fully explained
   by continuum dilution by re-fitting with the arc spectrum as a
   free-amplitude `sky` template (recovery fraction 145%).
2. **Different SPS template libraries disagree on V_sys.** FSPS ships
   in vacuum wavelengths; EMILES and XSL ship in air. Treating them
   uniformly (the bug we found and fixed) caused a ~99 km/s split in
   the recovered systemic velocity. Frame-aware ppxf collapses that
   split to ~8 km/s and halves the cross-check error bars.
3. **Aperture choice matters.** We use cumulative I-weighted ppxf
   inside R<R_e — the literature-standard σ_e definition that
   matches what KH13 / Greene+20 / SAURON / ATLAS3D / MaNGA use.
   Validated against three architecturally-independent annular paths
   (§7 discrete Gültekin sum, §7b flat-σ extrapolation, arc-spectrum
   subtraction) and four R_e definitions — all consistent within ±1σ.

## 6. Why it matters / "so what" framing

- **First σ_e for a strong-lens deflector at this redshift with this
  precision and systematic propagation.** Most published lens-deflector
  σ measurements use single SPS, single aperture, no mask sensitivity,
  no frame audit — and quote ~10% uncertainties that are likely
  optimistic.
- **Anchors the M•–σ relation at the strong-lens-deflector regime.**
  Combined with the lens-model derived enclosed mass, gives a 4D
  constraint (σ, M_enclosed, M★, z) on the central black hole + bulge
  scaling relations.
- **Templates a reusable pipeline.** Frame-aware multi-SPS ppxf,
  arc-mask sensitivity propagation, R_e-source systematic, all
  catalogued at `TESTS_AND_DIAGNOSTICS.pdf` (513 lines, 9 categories
  A-I, ~40 individual tests). Future AGEL deflectors can be
  pipelined through this with N≈1 night of human time.

## 7. Caveats to be transparent about

- σ_e at R<R_e/8 was dropped because that aperture (=0.288″) sits
  inside half the seeing FWHM (=0.64″). We do NOT report a central σ.
- The Bagpipes SED has a 8-11% residual in the NIR bands (F140W,
  F150W2). A single-τ exponential SFH is our prior; a more flexible
  SFH would tighten that fit. Effect on M★ is at the ±0.07-0.15 dex
  level which we already quote.
- The KCWI cube header is mislabeled (`DATE-OBS = 2025-11-17`,
  `PROGNAME = U002`, `PROGPI = Jones`) — that describes only the last
  stacking pass, not the observing provenance. The actual data
  combines Aug/Sept/Dec 2025 nights. Cite the multi-night combine,
  not the header. (See `NOTES_methodology_2026-04-27.md`.)

## 8. Reusable assets

| Asset | Path | Notes |
|---|---|---|
| Headline notebook (paper-driving) | `notebooks/09_final_sigma_e_paper.ipynb` | All headline + cross-check + systematic figures inline |
| Master test catalogue | `TESTS_AND_DIAGNOSTICS.{md,pdf}` | 9 categories A-I, ~40 tests, paths to caches |
| Production driver | `scripts/final_sigma_e.py` | One-shot N=500 production at all 3 mask tracks |
| Audit script | `scripts/audit_ppxf_methodology.py` | Four-orthogonal correctness audits (LSF, frame, z×air-vac, fwhm-dict) |
| Cosmology | FlatLambdaCDM(H0=70, Om0=0.3) | 1″ ≈ 7.04 kpc at z=0.67564 |

## 9. References to cite

- Cappellari & Emsellem (2004) — ppxf
- Cappellari (2017, 2023) — ppxf revisions
- Cappellari et al. (2006) — Cappellari+2006 eq. 1 (σ_e definition)
- Kormendy & Ho (2013) — M•–σ scaling for ellipticals
- Greene et al. (2020 ARA&A) — σ_e aperture definition + M•–σ updates
- Graham & Driver (2005) — Sersic 2D total-flux formula
- AGEL DR2 (Carleton et al., paper in prep / arXiv:24XX.YYYY)

## 10. Suggested poster layout (3 columns × 4 rows)

```
┌──────────────────┬──────────────────┬──────────────────┐
│ Title + authors  │   System sketch  │   Headline box:  │
│   + AGEL logo    │  (HST color cut, │  σ_e = 267 ± 12  │
│                  │   arc + lens)    │  log M★ = 11.33  │
├──────────────────┼──────────────────┼──────────────────┤
│ Motivation       │   Method (ppxf+  │   Hero figure 1: │
│ • M•–σ at lens   │   cumulative I-  │   σ_e vs aperture│
│   deflectors     │   weighted)      │  (nb09_sigma_vs_ │
│                  │                  │   aperture.png)  │
├──────────────────┼──────────────────┼──────────────────┤
│ Key plot:        │   Hero figure 2: │   Cross-checks   │
│ KCWI cube +      │   SED + M★ panels│  (nb09_cross-    │
│ deflector spec   │  (nb09_sed_and_  │   checks.png)    │
│                  │   Mstar.png)     │                  │
├──────────────────┼──────────────────┼──────────────────┤
│ Error budget     │   R_e + mask     │   Significance + │
│ (bar chart of    │   systematics    │   refs           │
│  5 components)   │  (two small      │                  │
│                  │   panels)        │                  │
└──────────────────┴──────────────────┴──────────────────┘
```

For talk slides: ~12 slides, one per section above (1, 2, 3, 4, 5,
hero-1, hero-2, 7-row error budget table, two systematics, caveats,
significance, refs).
