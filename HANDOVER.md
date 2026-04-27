# HANDOVER — AGEL0206 ApJL σ_e analysis (last updated 2026-04-27)

This document captures the state of the AGEL0206 stellar velocity dispersion analysis at the point of context compaction. Read this first when resuming.

---

## 1. Headline numbers (paper-ready, N=500)

**σ_e(<R_e) = 267.32 ± 24 km/s** — for the M•–σ relation (Kormendy & Ho 2013, Greene+2020)

Extracted via `notebooks/07c_sigma_e_equalN.ipynb`, §6cum (cumulative I-weighted ppxf), combined SPS pool of FSPS + EMILES + XSL, N_BOOTSTRAP = 500.

Cross-checks (agree at <1σ):
- §7 discrete Gültekin annular sum (arc-filtered to R<R_safe=1.72″): **254.99 −24.2/+28.4 km/s**
- §7b flat-σ extrapolation into outer annulus: **271 −33/+35 km/s**
- nb07e arc-spectrum subtraction sibling: identical to nb07c within 0.1 km/s — arc dilution is **NOT** a detectable systematic at production statistics

Other apertures (cumulative):
- σ_e(<R_e/8 ≈ 0.76″) = 209.18 ± 20 km/s
- σ_e(<R_e/2 ≈ 1.10″) = 225.78 ± 18 km/s
- σ_e(<R_safe = 1.72″) = 238.78 ± 20 km/s

R_e = **2.305″ = 16.23 kpc** (mean of HST F140W + F200LP masked CoG).

Per-SPS spread at R_e: FSPS 253 / EMILES 268 / XSL 280 (spread 27 km/s, captured by combined pool ±24).

---

## 2. Cube provenance — header metadata is mislabeled (NOT a data issue)

**The cube we used IS the final, most-reduced data product** — `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` (253 MB, shape 3317×100×100, symlinked from `../velocity_dispersion_from_IFU/`). The σ_e numbers above are correct and do NOT need re-running.

**However, the FITS header metadata is mislabeled.** The headers say `DATE-OBS=2025-11-17`, `PROGNAME=U002`, `PROGPI=Jones`, suggesting a single Nov 17 frame, but the actual data combines multiple nights — at minimum:

- **Aug 29 2025 HST = Aug 30 2025 UTC** (Program K409, frames `kr250830_00091`–`00100+` at PAs 0°/45°/90°) — the canonical observing night per user
- Sept 29 2025 (Program K409)
- Dec 29 2024 (frames `kr241229_00092`–`00095`)
- Nov 17 2025 (Program U002 — only the most recent stacking pass, hence the header date)

**For the paper, do NOT cite Nov 17 / U002 / Jones from the header** — that's a stacking-pass artifact, not the observation provenance. Cite the multi-night combine including Aug 29 K409 as the canonical night.

A separate, larger combine (`DESJ0206-red_medium_combined.fits`, 469 MB, shape 4269×120×120, 5384–9652 Å) was discovered in the shared Google Drive and copied locally to `raw_KCWI/`. This may be a more recent reduction with wider FoV/wavelength coverage, but **we did not re-run on it** because the cube we have is already the validated final product. The local copy + provenance files are kept as a reference:

```
raw_KCWI/
├── red/DESJ0206-red_medium_combined.fits     (469 MB — alternate reduction, NOT used)
├── blue/DESJ0206_medium_combined.fits        (285 MB — alternate blue, NOT used)
└── provenance/
    ├── DESJ0206-red_medium_shifts.list
    ├── DESJ0206_medium_shifts.list
    ├── DESJ0206_medium_hdrReference.txt
    ├── desj0206-red.py / desj0206.py        (kcwiRedux combine scripts)
    ├── dec29_2024_rl8000.list / .par / .shift.list
    ├── K409_2025-08-30_UTC.html              (Aug 29 night observing log)
    └── K409_2025-09-29_UTC.html
```

The raw `.fits` frames (`kr250830_*.fits`, `kb250830_*.fits`) are NOT in the shared Drive — only night logs. They presumably live on a collaborator's local machine. For paper writing this doesn't matter; the reduced cube is what we cite.

### Open paper-writing questions

1. **PI of K409** (NOT Jones — that was U002 on Nov 17 only). Look up Keck program archive.
2. **Total on-source exposure** — enumerate from each night's log: ≥ 10 frames × 300 s on Aug 29 alone, plus Sept 29 + Dec 29 + Nov 17.
3. **Aug 29 night seeing** — extract from `provenance/K409_2025-08-30_UTC.html` header.
4. **Whether to also use the BLUE arm** — Mg b 5170 Å falls at observed λ ≈ 8660 Å (rest 5170 × 1.6756 = 8666 Å). BLUE cube ends at 5872 Å so Mg b is NOT in BLUE for z=0.676 — it's in RED. So BLUE is mostly redundant for σ; might still be useful for [O II] check on z.
5. **Should we re-run on the 469 MB alternate reduction** for sanity? Quick check (re-point symlink, re-run nb07c at N=500, ~1h) — but optional, not blocking.

---

## 3. Notebooks family — current state

All in `notebooks/`:

| Notebook | Purpose | Status |
|---|---|---|
| `06_final_sigma_Re_apertures.ipynb` | σ at <R_e/8, R_e/2, R_e via earlier combine-posterior path | Run, superseded by nb07 family |
| `07_sigma_e_radial_clean.ipynb` | First clean Gültekin pipeline (5 equal-width annuli) | N=500, σ_e = 257 ± 22 |
| `07a_sigma_e_sersic_I.ipynb` | 8 annuli + F200 mask + 2D Sersic I-fits + 6-way I-source | N=50 |
| `07b_sigma_e_5bins.ipynb` | 5-annulus version of 07a | N=50 |
| `07c_sigma_e_equalN.ipynb` | **5 equal-N inner bins + 1 outer flagged** (HEADLINE) | **N=500, σ_e = 267 ± 24** |
| `07d_sigma_e_forceful_mask.ipynb` | Dilated F200 arc mask experiment — NEGATIVE finding (over-masks) | Settled, do not revisit |
| `07e_sigma_e_arc_subtract.ipynb` | nb07c + per-annulus α_arc residual subtraction | **N=500, identical to nb07c within 0.1 km/s** |

Three architecturally-independent σ_e estimators in nb07c:
- §6cum: I-weighted aperture spectrum → ppxf (one fit per R_max). Headline = 267 ± 24.
- §7: per-annulus ppxf → discrete F_j × (V²+σ²) sum. Filtered = 255 ± 26.
- §7b: §7 + ann5 inherits (V, σ) from ann4 (bulge-physics extrapolation). = 271 ± 34.

§6cum is the recommended headline — see memory `reference_cumulative_vs_annular_sigma_e.md` for full pros/cons discussion.

Bootstrap caches at `results/annular_bootstrap_07{c,e}/` (production N=500). N=50 caches preserved at `..._N50_backup/`.

---

## 4. Memory files (in `~/.claude/projects/.../memory/`)

Index: `MEMORY.md`. Latest method-relevant references:

- `reference_kcwi_data_properties.md` — **CORRECTED 2026-04-27** with multi-night provenance, K409/Aug 29 details, observing log analysis
- `reference_cumulative_vs_annular_sigma_e.md` — when to use §6cum vs §7/§7b
- `reference_sps_systematic.md` — FSPS/EMILES/XSL pooling math, V_sys offsets
- `reference_gultekin_implementation.md` — exact discrete formula, code locations
- `reference_dilution_correction.md` — EW-based diagnostic (qualitative only, not used in paper)
- `feedback_sigma_e_estimator_independence.md` — don't cross-contaminate the three paths
- `feedback_deflector_morphology.md` — AGEL0206 is an elliptical bulge, σ(r) should DECREASE with R
- `feedback_masking_strategy.md` — F200LP mask for arc, F140W mask is for Sersic only
- `project_sigma_e_gultekin.md` — full chain narrative with N=500 numbers
- `project_nb07d_overmasking_finding.md` — negative finding, do not revisit
- `project_nb07e_arc_subtraction.md` — full method + N=500 result table

---

## 5. Recent git history (this session)

```
1a558de Correct N=500 σ_e numbers; document cumulative-vs-annular choice
cf97113 N=500 production: σ_e(<R_e) = 267 ± 24 km/s (nb07c+nb07e)
1c6f8c6 Bump nb07c and nb07e to N_BOOTSTRAP=500 for production run
f47c74f Run nb07e smoke (N=50): arc subtraction → σ_e essentially unchanged
e3937c2 Build σ_e Gültekin pipeline (nb06, nb07/a/b/c/d, nb07e-skeleton)
```

`raw_KCWI/` directory NOT yet committed (large binary FITS files, blocked by `.gitignore` *.fits rule). Provenance files in `raw_KCWI/provenance/` are smaller and could be committed if useful — most are .list/.par/.html which aren't in `.gitignore`.

---

## 6. Suggested next steps (priority order)

1. **Get K409 PI + total exposure** from Keck program archive or PI contact. Update CLAUDE.md and memory boilerplate so the paper cites the multi-night provenance, not the mislabeled Nov 17 / U002 / Jones header.
2. **Extract per-night seeing** from K409 logs (Aug 29 + Sept 29) and Nov 17 + Dec 29 frame headers.
3. **Paper draft** — use updated boilerplate from `reference_kcwi_data_properties.md` once §1 questions are resolved. Headline σ_e numbers stand.
4. (Optional) Re-run nb07c at N=500 on the 469 MB alternate reduction `raw_KCWI/red/DESJ0206-red_medium_combined.fits` as a sanity check. Re-point symlink + re-extract spectra (~1h). Not blocking — only do if the alternate reduction is judged to be the canonical paper product.
5. **Decide cube symlink convention** — keep current `Nov17_2025_*` filename (misleading but stable) or rename to `DESJ0206_RL_combined.fits`.

---

## 7. Tool/environment notes

- ISMgas conda env: `source /opt/anaconda3/etc/profile.d/conda.sh && conda activate ISMgas` (NOT miniconda3)
- Parallel bootstrap: `scripts/bootstrap_ppxf_parallel.py`, `n_jobs=8`, 12 perf cores available
- All notebooks executed in-place via `jupyter nbconvert --to notebook --execute --inplace`
- Run spreadsheets / observing logs live in Google Drive (CfA account) under `KCWI arcs / Run spreadsheets by date / 2025/` and `KCWI arcs / Observing logs + finders / Logs by date / 2025/`
- Raw KCWI frames are NOT shared on Google Drive — only logs and reduced/combined cubes
