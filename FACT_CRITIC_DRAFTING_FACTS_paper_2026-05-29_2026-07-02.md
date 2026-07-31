# Fact-critic report — DRAFTING_FACTS §3.4 "Structural nature of the deflector" (2026-07-02)

**Scope:** §3.4 only (structural classification: elliptical / bulge-vs-pseudobulge / cored-vs-coreless).
**Active adapters:** Zotero (citation, priority) + CrossRef-web (citation) + PAPER_VALUES/npz (number); web ON.
Overlay `local/agel-sources.md` active. Numbers routed to `results/PAPER_VALUES.json` + `results/*.npz`;
citations to Zotero then CrossRef. **No claim PASSed from memory.**

## Numbers (routed to registry / npz)

| id | claim | verdict | source |
|----|-------|---------|--------|
| N1 | σ_e = 267 km/s | ✅ PASS | `PAPER_VALUES sigma_e.central` = 267.31 |
| N2 | log M⋆ ≈ 11.46 | ✅ PASS | `PAPER_VALUES logMstar.central_10pct` = 11.463 |
| N3 | mass-weighted age 3.8 Gyr | ✅ PASS | `PAPER_VALUES sed.mass_weighted_age` = 3.78 |
| N4 | A_V ≈ 0.2 | ✅ PASS | `sed.dust_Av` = 0.211 / `mstar_headline_quiescent.npz prop_dust_Av` |
| N5 | sSFR ≈ −10.5 | ✅ PASS | `PAPER_VALUES sed.ssfr` = −10.53 |
| N6 | b/a ≈ 0.80–0.86 | ✅ PASS | `sersic_parameter_table.npz` F140W 0.86 / F150W2 0.80 / F322W2 0.85 |
| N7 | PA ≈ 4–11° | ✅ PASS | `sersic_parameter_table.npz` 4.3 / 10.7 / 5.5 |
| N8 | V ≈ 40 km/s (within R_e) | ✅ PASS | `radial_sigma_combined_posterior.npz V_p50` = 40.1 |
| N9 | V/σ ≈ 0.15 | ✅ PASS | 40.1 / 267.31 = 0.150 (derived from N1,N8) |
| N10 | single-Sérsic n ≈ 1.2–1.6 | ✅ PASS | `sersic_parameter_table.npz` 1.22–1.59 (F200LP 1.24) |
| N11 | F140W n 0.9→2.0 over box 4→10″ | ❌→**AUTO-FIXED** to **1.2→2.0** | nb18 Test A output: "n: 1.17→1.98 over box 4→10″" (0.9 was the box-3 value) |
| N12 | r_b < 0.024″ = 0.17 kpc | ✅ PASS | `core_sersic_test.npz` F150W2 HWHM 0.024″ × 7.2764 = 0.175 kpc |
| N13 | PSF F150W2 FWHM 0.049″ (HWHM 0.024″) | ✅ PASS | `core_sersic_test.npz F150W2_fwhm_as` = 0.049 |
| N14 | ΔBIC = 386 / 43 | ✅ PASS | `core_sersic_test.npz` dBIC 385.8 (F150W2) / 43.1 (F140W) |
| N15 | r_b,fit = 2.5–3.6 kpc | ✅ PASS | F140W 0.343″×7.2764=2.50 / F150W2 0.497″×7.2764=3.62 |
| N16 | "0.02–0.1 kpc = a few–40 mas" | ❌→**AUTO-FIXED** to **≈3–14 mas** | geometry: 0.02–0.1 kpc ÷ 7.2764 = 2.7–13.7 mas (the "40 mas" upper end was ~3× too high) |
| N17 | expected r_b ~ 0.02–0.1 kpc | ⚠ SUSPECT (provenance) | a *literature-expectation* range, not an in-repo number; should be anchored to a cited M•–r_b scaling (Lauer+2007 / Rusli+2013), not asserted. See X1. |

**Auto-fixes applied to the draft:** N11, N16 (both trace to a single canonical source for the
unambiguously-identified AGEL0206 deflector; conclusions unchanged — the core stays ≪ PSF).

## Citations (routed to Zotero → CrossRef)

| id | claim | verdict | source / action |
|----|-------|---------|-----------------|
| C1 | Kormendy & Ho 2013 (ARA&A 51, 511) — classical/pseudo criteria + pseudobulges off M•–σ | ✅ PASS | in Zotero (verified this session); the canonical review — support holds |
| C2 | Greene+2020 E+S0 calibration | ✅ PASS | Zotero; Suppl. Table 5 "Early" verified 2026-06-17 (`NOTES_relation_offset_provenance`) |
| C3 | Faber+1997 (core/coreless dichotomy) | ⚠ SUSPECT→**verified real, ADD** | NOT in Zotero. CrossRef: Faber et al. 1997, **AJ 114, 1771, DOI 10.1086/118606**. Add to library/.bib. |
| C4 | Lauer+2007 (core/coreless dichotomy) | ❌ **MIS-ATTRIBUTION RISK** | The Zotero "Lauer 2007" (item 3TUYW9FH) is the *M•–σ selection-bias* paper — **not** the cusp/core paper. The dichotomy paper is Lauer et al. 2007, **ApJ 664, 226, DOI 10.1086/519229** ("Bimodal Central Surface Brightness Profiles"). Cite/add THAT; ensure the §3.4 key does not resolve to the selection-bias entry. |
| C5 | Kormendy+2009 (core/coreless dichotomy) | ⚠ SUSPECT→**verified real, ADD** | NOT in Zotero (the "Kormendy 2009" hit is Gültekin+2009). CrossRef: Kormendy, Fisher, Cornell & Bender 2009, **ApJS 182, 216, DOI 10.1088/0067-0049/182/1/216**. Add. |
| C6 | pseudobulge criteria / mass & σ / V/σ classification (currently → KH13 only) | ⚠ SUSPECT (needs primary ref)→**ADD** | KH13 covers it, but the *primary* refs are **Fisher & Drory 2008, AJ 136, 773, DOI 10.1088/0004-6256/136/2/773** (links pseudobulges to Sérsic n<2; low-n bulges → higher V/σ, lower σ — directly supports our n/V-σ argument) and **Fisher & Drory 2016, ASSL 418, 41, DOI 10.1007/978-3-319-19378-6_3** (observational guide). Add + cite for the criteria. |

(Per policy, citations are flag-only — none auto-edited.)

## Consistency

- **X1** — expected-core range disagrees across docs: §3.4 says r_b ~ 0.02–0.1 kpc (=3–14 mas, post-fix)
  while `REPORT_cored_vs_coreless_2026-06-17.md` says "r_b ~ 0.005–0.05″" (=0.036–0.36 kpc). ~3× mismatch.
  Both are below the PSF (conclusion — unresolvable — unaffected), but the number should be reconciled to
  ONE cited M•–r_b scaling. **Flag** (not auto-fixed: >1 differing candidate, and it's a literature value).

## Verdict

§3.4's **in-repo numbers are sound** (15/17 PASS; 2 auto-fixed typos, conclusions unchanged). The
**citation layer is the weak point**: C1/C2 solid, but **C3 (Faber97), C4 (Lauer07 — wrong paper in
library!), C5 (Kormendy09) and C6 (Fisher & Drory primary refs) must be added/corrected** before §3.4 is
drafted. All four are CrossRef-verified real papers with DOIs above — handed to the Zotero curation pass.
One provenance flag (N17/X1): pin the expected-core radius to a cited scaling rather than asserting it.
