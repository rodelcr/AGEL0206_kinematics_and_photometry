# HANDOFF — project stock-take before machine update (2026-07-30)

**Purpose:** the laptop is being updated. This doc is (a) the state-of-the-project snapshot to
pick up from afterward, and (b) a **migration risk register** — what is *not* protected by git
and would be lost if the disk is wiped.

Previous handoff: `HANDOFF_structure_relations_tutorials_figures_2026-07-06.md`.
Everything below was verified by inspection on 2026-07-30, not from memory.

---

## 0. TL;DR — where the science stands

All headline numbers come from `results/PAPER_VALUES.json` (registry, single source of truth).
**Drift linter is GREEN** as of today:

```
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md \
    TESTS_AND_DIAGNOSTICS.md METHODS_AND_SYSTEMATICS.md PROJECT_BRIEF.md
→ check: OK — all headline statements match the registry (5 live files, 0 exempt)
```

| quantity | value | note |
|---|---|---|
| σ_e(<R_e) | **267.31 ± 4.51 (stat) ± 10.87 (sys)** km/s | total ±11.77 sym (−11.98/+11.58); stat asym −5.04/+3.98 |
| R_e | **2.097″ = 15.26 kpc** | best-mask CoG mean; method sys ±0.100″ |
| log M★/M☉ | **11.463** (+0.069/−0.057 stat) ± 0.168 sys | validated photom + quiescent SED prior, 2 R_e, 10% floor |
| M_BH | **5.2 (+1.7/−1.4) × 10⁸ M☉**, log = 8.716 | ⚠ **PROVISIONAL** — Ferrami free-BH EPL envelope |
| θ_E / γ | 1.36″ / 1.31 ± 0.08 | ⚠ PROVISIONAL (Ferrami draft §3.3.2) |
| SED (passive) | age 3.78 Gyr, SFR 7.9, log M_formed 11.73 | quiescent prior |

The science is stable. **The risk is entirely on the storage side.**

---

## 1. ⚠ MIGRATION RISK REGISTER — read before wiping anything

Ordered most-severe first. Nothing here is done yet.

### 1.1 Unpushed git history — 83 + 7 commits

| repo | branch | state |
|---|---|---|
| `AGEL0206_kinematics_and_photometry` | `photometry-masking-mstar-2026-05-29` | **83 commits ahead of `origin/main`**, 0 behind |
| " | `main` | 71 commits ahead of `origin/main` |
| `AGEL_0206_ApJL_Figures` | `main` | **7 commits ahead of `origin/main`** |

Remote is `git@github.com:rodelcr/AGEL0206_kinematics_and_photometry.git` (SSH — the **new machine
needs the SSH key**, or push before wiping). Nothing has been pushed since the best-mask cascade.
**This is the single highest-value, lowest-effort save.**

### 1.2 Uncommitted work in the main repo

11 tracked files modified, **+3435 / −505 lines**, none staged:

```
DRAFTING_FACTS_paper_2026-05-29.tex   +2489 / −358   ← the big one
DRAFTING_FACTS_paper_2026-05-29.md     +558 /  −78
scripts/paper_values.py                +274 /  −26
TESTS_AND_DIAGNOSTICS.md                +76 /   −8
METHODS_AND_SYSTEMATICS.md, PROJECT_BRIEF.md, CLAUDE.md,
scripts/{bagpipes_sersic_refit,sersic_parameter_table}.py, 2 PDFs
```

Plus **~30 untracked files that have never been committed**, including 12 analysis scripts
(`aperture_correction_validated.py`, `bh_mass_combine.py`, `build_psf_models.py`,
`cog_total_light.py`, `core_sersic_test.py`, `greene_sample_redshift.py`, `make_latex_tables.py`,
`mstar_drop_f200lp.py`, `mstar_headline_quiescent.py`, `psf_star_census.py`,
`relation_offset_significance.py`, `sed_properties_final.py`, `sed_quiescence_check.py`),
6 HANDOFF/NOTES/REPORT docs, `notebooks/18_cored_vs_coreless_sersic.ipynb`, and
`cosmology_reference_AGEL0206.ipynb`. **`mstar_headline_quiescent.py` produces the M★ headline** —
it is untracked.

The ApJL repo likewise has 9 modified figures/notebooks + ~40 untracked files.

### 1.3 Precious data that git will NEVER protect

`.gitignore` excludes `*.fits *.npz *.npy *.h5 *.zip results/ pipes/ TEXT/ logs/
notebooks/tutorials/`. So **pushing does not back these up**:

| path | size | why it matters |
|---|---|---|
| `results/` | **601 MB, 1129 files** | every N=500 bootstrap cache. Re-running is ~12 min *per SPS per config* — days of compute |
| `raw_KCWI/` | **1.0 GB** | `New_red/` (253M) holds the `_mtwdo_` headline cube; also `blue/ red/ provenance/` |
| `notebooks/` | 164 MB | incl. `notebooks/tutorials/` — gitignored, so **the tutorial suite + the 2026-07-27 README fix exist only on this disk and on Drive** |
| `logs/` | 47 MB | sweep logs |
| `figures/`, `results/figures/` | 12 + 64 MB | |
| repo total | **2.4 GB** | |

### 1.4 The sibling directory the repo depends on

Six symlinks at the repo root point **outside** the repo, into
`../velocity_dispersion_from_IFU/` (**4.7 GB, not a git repo at all**):

```
TEXT -> ../velocity_dispersion_from_IFU/TEXT                     (622 MB — veldis templates)
pipes -> ../velocity_dispersion_from_IFU/pipes                   (3.3 MB — Bagpipes posteriors)
spectra_{fsps,emiles,xsl}_9.0.npz -> ../velocity_dispersion_from_IFU/...
Nov17_2025_..._icubes_wcs.fits -> ../velocity_dispersion_from_IFU/...
```

All six resolve today. **Copy `velocity_dispersion_from_IFU/` alongside the repo or every link
breaks** and ppxf/Bagpipes stop working. Preserve the relative layout
(`~/Documents/AGEL/<both dirs>`), and copy with `cp -R` / rsync **without** `-L` if you want the
links to stay links, or resolve them deliberately.

`AGEL_0206_ApJL_Figures/` (250 MB) also reads
`../AGEL0206_kinematics_and_photometry/results/PAPER_VALUES.json` by relative path — same
constraint.

### 1.5 What is already safe

- **Google Drive** `My Drive/AGEL/ppxf_redshift_tutorials/` — full tutorial bundle incl. the
  265 MB cube; synced and cloud-backed. Verified md5-identical to local on 2026-07-27.
- `~/.claude/projects/.../memory/` — cross-session memory (separate from this repo; confirm it
  rides along or is re-synced).

### 1.6 Suggested order of operations

1. `git push` both repos (branch + main) — *needs your go-ahead; I have not pushed.*
2. Commit or explicitly stash the 11 modified + ~30 untracked files (the 12 scripts especially).
3. Rsync `results/`, `raw_KCWI/`, `notebooks/`, `logs/` to external disk **or** Drive.
4. Rsync all of `velocity_dispersion_from_IFU/` (4.7 GB).
5. Rsync `AGEL_0206_ApJL_Figures/` (250 MB) after pushing it.
6. Export conda envs (§4) — do **not** try to copy `/opt/anaconda3/envs` across an OS update.
7. Confirm the SSH key for GitHub exists on the new machine before you need it.

---

## 2. What changed since the 2026-07-06 handoff

### 2.1 Violin M_BH figures for the HST DDT proposal (2026-07-30, today)

New section at the bottom of `AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb`
(**markdown cell 30 + code cells 31–32**), titled *"Figs 3 & 4 — AGEL0206 M_BH as violins:
Archival HST vs Proposed DDT"*. Styled to match **Giovanni's HST DDT-proposal panel**:

- AGEL0206 M_BH drawn as **two violins** from `data_ddt_converged.npz` (keys `x_old/p_old`,
  `x_ddt/p_ddt`; 256-point PDFs): **Archival HST** = orange `#F0A44A` (current, broad,
  upper-limit-like posterior, violin spans central **4σ**) vs **Proposed DDT** = navy `#161670`
  (converged, spans central **3σ**). Tails deliberately clipped, not full-range.
- Only the **Greene+2020** relation retained — KH13 / van den Bosch / Reines / Sahu / Farrah
  dropped for this panel. Local ETGs gray open circles, lensing comparisons green `#8DC63F`.
- Outputs: `Mbh_sigma_relation_scatterband_violin.pdf`,
  `Mbh_Mstar_relation_clean_scatterband_violin.pdf`. Backup of the pre-edit notebook:
  `figures_paper4_2026-06-08.ipynb.bak.preViolin_20260730`.
- σ_e is read live from `../AGEL0206_kinematics_and_photometry/results/PAPER_VALUES.json`.
- Run in env **`agel0206_figs`**. All of this is **uncommitted**.

An earlier upper-limit variant exists from 2026-07-13
(`*_scatterband_upperlimit.pdf`, `.bak.preUL_20260713`).

### 2.2 ppxf tutorial bundle — 05/06 scoping fixed (2026-07-27)

Collaborator hand-off pass. Tutorials 05–06 were documented as merely "additionally needing"
repo files; they in fact **cannot run from the bundle at all**. Corrected in
`notebooks/tutorials/README.md` ("Two tiers" now enumerates every missing dependency in a table:
the three `scripts/` modules, `PAPER_VALUES.json`, the two N=500 caches, the `_mtwdo_` cube
flagged as a *different* cube from `data/`, HST imaging + arc masks), plus a new
"Which download do I want?" section, corrected troubleshooting, and matching edits to
`START_HERE.md` and Drive-only `START_HERE.txt` (whose "ignore the code zip" line was wrong).
Propagated to Drive + both zips; **all five copies verified md5-identical**.
⚠ `notebooks/tutorials/` is gitignored → this work is **not in git**.

### 2.3 Fact-critic pass on DRAFTING §3.4 (2026-07-02)

`FACT_CRITIC_DRAFTING_FACTS_paper_2026-05-29_2026-07-02.md`. Numbers sound (15/17 PASS, 2 typos
auto-fixed, conclusions unchanged). **Citation layer is the weak point** — see §3.

---

## 3. Outstanding TODOs

**New / highest priority (from the fact-critic, blocking §3.4 drafting):**

- **C4 — wrong paper in Zotero.** The "Lauer 2007" item (`3TUYW9FH`) is the *M•–σ selection-bias*
  paper, **not** the cusp/core one. The dichotomy paper is Lauer et al. 2007, **ApJ 664, 226,
  DOI 10.1086/519229**. Fix before citing.
- **C3 / C5 / C6 — verified real, missing from Zotero.** Faber+1997 (AJ 114, 1771,
  10.1086/118606); Kormendy, Fisher, Cornell & Bender 2009 (ApJS 182, 216,
  10.1088/0067-0049/182/1/216); Fisher & Drory 2008 (AJ 136, 773, 10.1088/0004-6256/136/2/773)
  and 2016 (ASSL 418, 41, 10.1007/978-3-319-19378-6_3) as the *primary* pseudobulge refs.
- **N17 / X1 — expected core radius unanchored.** §3.4 says r_b ~ 0.02–0.1 kpc (≈3–14 mas);
  `REPORT_cored_vs_coreless_2026-06-17.md` says 0.005–0.05″ (0.036–0.36 kpc). ~3× mismatch, both
  below PSF so the "unresolved" conclusion is unaffected, but pin it to ONE cited M•–r_b scaling
  (Lauer+2007 / Rusli+2013).

**Carried from 2026-07-06:**

- **M_BH stays PROVISIONAL** until Ferrami resolves free-vs-fixed BH. Registry has
  `lens_model.provisional = 1` — flip it when settled.
- **Registry reconcile (needs your OK):** `mbh_sigma_offset` −0.69 → mean −0.26; add
  `mbh_mstar_offset` −0.49 + `n_sigma` leaves.
- **`DLOGM_COSMO` cleanup:** regenerate `aperture_matched_photometry.npz` natively under Planck
  2015 instead of applying the +0.0282 dex rescale.
- **DRAFTING §2.4.2** pipeline walkthrough still pre-audit.
- **KCWI provenance:** all 3 K409 nights confirmed observed, but whether Aug 30 + Dec 29 were
  actually folded into the headline combined cube is **UNVERIFIED** — ask Kaustubh/Yuguang before
  claiming a "3-night combine". See `HANDOFF_kcwi_provenance_triple_check_2026-07-06.md`.
- **Hδ** targeted masking still flagged for revisit (`METHODS_AND_SYSTEMATICS.md` Part III.5 #0).

---

## 4. Environments to rebuild on the new machine

Conda root is `/opt/anaconda3` (25 envs). The three that matter here — **do not copy the
directories, re-create them**:

| env | size | used for |
|---|---|---|
| **`ISMgas`** | 1.4 GB | everything in this repo: ppxf, bootstraps, `scripts/*`, notebooks |
| **`agel0206_figs`** | 1.5 GB | `AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb` |
| **`ppxf_tutorials`** | 1.0 GB | the collaborator tutorial bundle (rebuildable from its `environment.yml`) |

`ISMgas` as of today: **Python 3.14.2, ppxf 9.4.5, astropy 7.2.0, numpy 2.2.0**.
Export before wiping:

```bash
for e in ISMgas agel0206_figs ppxf_tutorials; do
  conda env export -n $e > ~/env_backup_${e}_2026-07-30.yml
  conda list -n $e --explicit > ~/env_explicit_${e}_2026-07-30.txt
done
```

(`ppxf_tutorials` already has a curated `notebooks/tutorials/environment.yml` — prefer that one
for the collaborator bundle; the export above is for exact local reproduction.)

Note the ppxf version pin: **9.4.5** is what produced every number in the registry, and
`NOTES_ppxf_reporting_parameters_2026-06-22.md` records the exact call parameters
(moments=2, mdegree=0, regul=0, trig=False, degree 15–29, 2161 good pixels, FWHM_inst 0.692 Å).

---

## 5. Resume / verify on the new machine

```bash
conda activate ISMgas
cd ~/Documents/AGEL/AGEL0206_kinematics_and_photometry

# 0. symlinks resolve? (fails loudly if velocity_dispersion_from_IFU didn't come along)
ls -l TEXT pipes spectra_*.npz Nov17_2025_*.fits

# 1. registry + drift linter — must print GREEN
python scripts/paper_values.py --check DRAFTING_FACTS_paper_2026-05-29.md CLAUDE.md \
    TESTS_AND_DIAGNOSTICS.md METHODS_AND_SYSTEMATICS.md PROJECT_BRIEF.md

# 2. cheap analysis re-runs (seconds; no bootstrap)
python -m scripts.relation_offset_significance    # 0.57σ / 0.69σ
python -m scripts.greene_sample_redshift          # A1836-BCG z≈0.034
python -m scripts.core_sersic_test                # core unresolved

# 3. caches present? expect 1129 files / 601 MB
find results -type f | wc -l && du -sh results

# 4. figures: AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb (env agel0206_figs)
#    violin cells 30–32 at the bottom; needs ../AGEL0206.../results/PAPER_VALUES.json
```

State docs to re-read in order: `CLAUDE.md` → `TESTS_AND_DIAGNOSTICS.md` →
`HANDOFF_structure_relations_tutorials_figures_2026-07-06.md` → this file → memory `MEMORY.md`.
