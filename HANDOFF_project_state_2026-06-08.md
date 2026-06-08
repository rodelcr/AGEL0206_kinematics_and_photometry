<!-- pv-skip-file: dated snapshot — records historical (now-superseded) numbers; current headline in CLAUDE.md / results/PAPER_VALUES.json -->

# HANDOFF — AGEL0206 project state & remaining work — 2026-06-08

> **⚠ SUPERSEDED LATER THE SAME DAY (M12).** This is the *morning* remaining-work audit. The Tier-1/2
> tasks below were then executed → **current headline σ_e = 269.62 ± 13.27 km/s** (R_e-source ±6.13
> folded in). See **`HANDOFF_systematics_registry_2026-06-08.md`** for the closeout + the
> `scripts/paper_values.py` registry. The σ_e number in this file's body (±11.77) is pre-M12.

General audit of what's done and what's left, across the whole ApJL effort (kinematics +
photometry + drafting + figures). Supersedes the stale `project_roadmap` memory (dated 2026-05-11,
which still shows the pre-M8/M10/M11 σ_e and pre-principled M★). For the most recent work arc see
`HANDOFF_photometry_masking_Mstar_2026-05-29.md`.

## Current headline numbers (as of the morning audit — σ_e since revised to ±13.27 by M12)

- **σ_e(<R_e) = 269.62 ± 11.77 km/s** (asym −11.98 / +11.57) — wide arc-masked, NEW `_mtwdo_`
  reduction, He I + M10 sky masks, M11 budget. **[pre-M12; now ±13.27 — see banner above]**
- **log(M⋆/M☉) = 11.16 ± 0.08 (stat) +0.31 (sys)** at 10% flux errors [11.04 ± 0.14 +0.32 at 20%]
  — NEW principled arc masking (F200LP-located + IR-extended, raw aperture), one-sided +sys =
  under-arc light (fill-in reach 11.46). Explicit masking systematic ±0.16 dex. [2026-05-29]
- **R_e = 2.305″ = 16.23 kpc** (F140W+F200LP CoG mean, final_sigma_e.py).
- z_deflector = 0.67564; z_source = 1.302 (error bar still [TBD]).

## Git state

- On branch **`photometry-masking-mstar-2026-05-29`** (commit `122745c`) — NOT merged to main.
- Pre-existing unrelated working-tree mods left untouched (bootstrap_ppxf.py, run_isource_shape_sweep.py,
  nb09, HANDOVER.md, NOTES_methodology) — these predate this arc; decide separately.
- ApJL figures repo (`../AGEL_0206_ApJL_Figures`): Fig 2 + Fig 4 edited for the new M★ but NOT
  committed (had pre-existing staged state + my headless run left aplpy-error outputs in Fig-1 cells).

## REMAINING WORK (prioritized)

### Tier 1 — close out the photometry/M★ arc
1. **σ_e masking systematic** (the analogue of the ±0.16 dex M★ one): reproject each masking
   approach (F200LP-located, IR-extended per-band, global) to the IFU grid and re-run
   `scripts/run_wide_sigma_e.py --cube new_clean_hei` under each → quote the spread. *Not done.*
2. **Commit the ApJL figures** (Fig 2/Fig 4): clean the aplpy-error outputs (re-run interactively in
   the `lens` env), resolve the pre-staged state, commit in that repo. *User to do.*
3. **Merge branch → main** once the above are settled.
4. **CoG-algorithm reconciliation**: `measure_Re.hst_Re` (1.91″) vs `final_sigma_e.curve_of_growth`
   (2.52″) for F200LP disagree by ~0.6″ even at the same pixscale (different annulus binning).
   Headline uses final_sigma_e; reconcile or document which is canonical.

### Tier 2 — kinematics systematics still open
5. **R_e-source systematic at the wide window** (D7): currently 16.9 km/s measured at the *narrow*
   window and **NOT folded into the wide budget** — flagged for a wide-window re-measurement.
6. **Hδ (4101.74 Å) targeted masking** — open decision; TODO in `bootstrap_ppxf.py:262`
   (`_determine_goodpixels_no_balmer`) + METHODS Part III.5 #0 + DRAFTING_FACTS §line 281.
7. **Reduction-pass systematic** (±3.45, only 2 reductions) — refine if a 3rd reduction lands. Passive.

### Tier 3 — manuscript drafting (DRAFTING_FACTS gaps)
8. **Photometry / M⋆ section** — write with the new principled masking + 11.16 budget (replaces the
   old 11.33 aperture text).
9. **Kinematics section** — method (wide arc-masked + source-emission catalog) / results
   (269.62 ± 11.77) / systematics (M11 budget table). The SPS-spread-collapse 26→4 km/s is the
   strong methodological argument.
10. **Discussion** — rising outer σ profile (halo signature, not contamination); compare to
    Melo-Carneiro+2025 and Smith+2017.
11. **G1 — lens modeling (Ferrami et al.)**: BPL/point-mass/M(<r_break)/joint-posterior/PSF/shear —
    NOT in this repo; source from `202509_DESJ0206_modeling/` or the Ferrami draft.
12. **G2 — X-ray / "truly quiescent"**: no X-ray data in repo; source externally if claimed.
13. **G3 — "we do not account for peculiar velocities"**: true statement; just add the sentence.
14. **Source-z error bar** (z_source = 1.302 ± [TBD]) + its cache — placeholder to fill from the
    KCRM blue-arm cube.

### Tier 4 — housekeeping
15. **Memory:** `project_roadmap` updated this session (was stale); `project_apjl_paper` Fig 2 still
    says "IN PROGRESS" — now updated on both σ_e (wide) and M⋆ (principled) axes, pending figures-repo commit.
16. **Git org:** merge branch, push both repos, clean untracked (TEXT/, pipes/, logs/, raw_KCWI/ are
    data/outputs not meant for git).

## Key entry points for next session

- This file + `HANDOFF_photometry_masking_Mstar_2026-05-29.md` (recent arc).
- `DRAFTING_FACTS_paper_2026-05-29.md` — the drafting fact sheet (headline numbers + gaps).
- `TESTS_AND_DIAGNOSTICS.{md,pdf}` — test catalog (A–N series; §5 "what's NOT done").
- `CLAUDE.md` — current headline + systematic budgets + the arc-mask section.
- Memory `MEMORY.md` → `project_arc_mask_verification`, `project_2026-05-27_he_i_and_m10_mask_audit`.
- iCloud `~/Library/Mobile Documents/com~apple~CloudDocs/AGEL0206_Draft/` — phone-readable copies.
