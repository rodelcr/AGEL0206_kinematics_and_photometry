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

## 2. ⚠️ CRITICAL DATA-PROVENANCE ISSUE — must address before paper submission

**The σ_e numbers above were computed on the WRONG cube** — a single 300s frame from UT 2025-11-17 (`Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`, 253 MB, shape 3317×100×100, single Nov17 frame from program U002 PI Jones).

**The CORRECT paper cube is** the multi-night combined product `DESJ0206-red_medium_combined.fits` (469 MB, shape 4269×120×120, 5384–9652 Å) — combines frames from at least 4 nights:

- **Aug 29 2025 HST = Aug 30 2025 UTC** (Program K409, frames `kr250830_00091`–`00100+` at PAs 0°/45°/90°) — the canonical observing night per user
- Sept 29 2025 (Program K409)
- Dec 29 2024 (frames `kr241229_00092`–`00095`)
- Nov 17 2025 (Program U002, the frame that mistakenly became the headline)

Both cubes (red + blue) and provenance files (shifts.list, .py reduction scripts, observing logs) have been copied locally to:

```
raw_KCWI/
├── red/DESJ0206-red_medium_combined.fits     (469 MB — paper headline)
├── blue/DESJ0206_medium_combined.fits        (285 MB — Mg b/Hβ region in BLUE arm)
└── provenance/
    ├── DESJ0206-red_medium_shifts.list
    ├── DESJ0206_medium_shifts.list
    ├── DESJ0206_medium_hdrReference.txt
    ├── desj0206-red.py / desj0206.py        (kcwiRedux combine scripts)
    ├── dec29_2024_rl8000.list / .par / .shift.list
    ├── K409_2025-08-30_UTC.html              (Aug 29 night observing log)
    └── K409_2025-09-29_UTC.html
```

The actual `.fits` raw frames (`kr250830_*.fits`, `kb250830_*.fits`) are **NOT** in the Google Drive — only the night log is shared. They presumably live on a collaborator's local machine (e.g., Yuguang Chen's `~/obs/2025aug29/`). For paper writing, the combined cube is enough; for re-reduction, get raw frames from PI/reducer.

### What changes if we re-run on the multi-night cube

The σ_e pipeline (notebooks 06, 07, 07a–07e) reads the cube via the symlink at the repo root → `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`. To switch:

1. Re-point the symlink (or change `IFU_FILE` config) to `raw_KCWI/red/DESJ0206-red_medium_combined.fits`
2. Re-extract integrated/aperture/annular spectra (cube has different shape and wavelength range)
3. Re-run N=500 bootstrap fits in nb07c + nb07e → ~2h sequential
4. The headline σ_e value will likely shift slightly (more frames = higher S/N, possibly tighter ±) but the three-method agreement should hold

### Open paper-writing questions

1. **PI of K409** (NOT Jones — that was U002). Look up Keck program archive.
2. **Total on-source exposure** — enumerate from each night's log: ≥ 10 frames × 300 s on Aug 29 alone.
3. **Aug 29 night seeing** — extract from `provenance/K409_2025-08-30_UTC.html` header.
4. **Whether to also use the BLUE arm** — Mg b 5170 Å falls at observed λ ≈ 8660 Å (rest 5170 × 1.6756 = 8666 Å). BLUE cube ends at 5872 Å so Mg b is NOT in BLUE for z=0.676 — it's in RED. So BLUE is mostly redundant for σ; might still be useful for [O II] check on z.

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

1. **Re-point pipeline to multi-night combined cube** in `raw_KCWI/red/DESJ0206-red_medium_combined.fits` and re-run nb07c at N=500 (~1h). Compare headline σ_e shift.
2. **Get K409 PI + total exposure** from Keck program archive or PI contact. Update CLAUDE.md and memory boilerplate.
3. **Extract per-night seeing** from K409 logs (Aug 29 + Sept 29) and Nov 17 + Dec 29 frame headers.
4. **Decide on cube symlink convention** — should the canonical cube live at the repo root (rename current symlink target) or stay in `raw_KCWI/`?
5. **Paper draft** — use updated boilerplate from `reference_kcwi_data_properties.md` once questions above are resolved.

---

## 7. Tool/environment notes

- ISMgas conda env: `source /opt/anaconda3/etc/profile.d/conda.sh && conda activate ISMgas` (NOT miniconda3)
- Parallel bootstrap: `scripts/bootstrap_ppxf_parallel.py`, `n_jobs=8`, 12 perf cores available
- All notebooks executed in-place via `jupyter nbconvert --to notebook --execute --inplace`
- Run spreadsheets / observing logs live in Google Drive (CfA account) under `KCWI arcs / Run spreadsheets by date / 2025/` and `KCWI arcs / Observing logs + finders / Logs by date / 2025/`
- Raw KCWI frames are NOT shared on Google Drive — only logs and reduced/combined cubes
