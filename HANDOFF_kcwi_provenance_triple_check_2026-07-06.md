<!-- pv-skip-file: dated snapshot of a provenance investigation, not a numerical-results file -->

# HANDOFF — KCWI spectroscopic data provenance triple-check (2026-07-06)

## TL;DR

Set out to nail down clear provenance (dates, program numbers, PIs) for the KCWI data behind the
σ_e headline. Result: **PIs resolved for all three nights** (K409 = Alcorn, new; U002 = Jones;
U204/Dec-29 = Jones, inferred). **Dec-29 2024 as an observation is now airtight** (two independent
primary Keck documents, read directly, not just referenced by id). But a deeper check surfaced a
real gap: **there is no direct, file-level proof that the Aug-30 (K409) and Dec-29 (U204) nights
were actually folded into the headline combined cube** that produces the paper's σ_e number — only
Nov-17 has that direct proof. This is the most important open item from this session.

## What was resolved

### 1. K409 PI = Alcorn (new, 2026-07-06)

Source: AGEL Airtable base (`Observing Runs` table), record `recH2BYEfLJ6iwMO5`,
`Proposal_ID = ALCORN_2025B_K409_KCWI`, **directly linked** to `AGEL020613-011417A` in the
`Spectral Observations Tally` table. The record's own `Principle_Investigator` field is blank —
the PI name is embedded only in the `Proposal_ID` string, consistent with every other record in
the base (`JONES_2025B_KCWI`, `JONES_2020B_U044_NIRES`, etc., which do have the explicit field
filled in AND matching). High confidence, not a directly-populated DB field.

Also confirmed via Airtable: **U002 (Nov 17) = Jones**, record `recUNQalnuxphIAwa`,
`Proposal_ID = JONES_2025B_KCWI`, date "16 & 17 November", directly linked to this target.

### 2. Dec 29, 2024 / U204 — observation itself is now airtight

Read both primary Keck documents directly off the Google Drive mount (not just cited by Drive id
as before):

- **Machine-readable autolog** `kcwi_autolog_20241228.pdf` — header: `Project Name: U204`,
  `Date: 2024-12-28`, `OutDir: /s/sdata1400/kcwi7/2024dec29`. Object rows explicitly list
  `kr241229_00092–00095.fits` (RED, `Object=DESJ0206`) at **RA 02:06:13.47, Dec −01:14:17.4** —
  exact match to AGEL0206 — plus `kb241229_00085.fits` (BLUE, 1320 s).
- **Handwritten observer log** `kcwi-241228-written-log.pdf` — despite being flagged "image-only,
  no OCR" in the 2026-05-26 note, it IS directly readable: page 1 header "PROJECT: U204", same
  date/outdir, observer "Rhoades". Page 3: `B85/R92 — Object: "DES 206" — 1320 B / 300.4 R`,
  matching the autolog exactly.

**PI for U204 is still not directly stated in either document** (Keck nightly logs record
Project/Observer/SA/OA, not PI — PI is a semester-level TAC attribute, not logged per-night).
Best evidence: Airtable's only U204 record (`recwiGq3B0Pso8hvd`, PI=Jones) is for a *different*
night (Sept 6 2024, unrelated program use); Jones is also on the Dec-29 observer list. Reasonable
same-program/same-semester inference, not a directly-logged fact for this specific night.

## The important open item: was Dec-29 (and Aug-30) actually folded into the headline cube?

Pushed one level deeper than "was it observed" to "did it end up in the science-grade combined
cube," and found the evidence is much weaker than the standing notes implied.

**What's solid:**
- Found the actual Dec-29-only reduced cube on Drive:
  `Shared drives/KCWI arcs/Arc data/stacks/desj0206/desj0206_rl8000_icubes.fits` (+ a BLUE
  `bl4500` counterpart). Opened it directly — `OBJECT=DESJ0206`, `DATE-OBS=2024-12-29`,
  `HISTORY` shows `kr241229_00092_crmsk.fits cleaned cosmic rays` as the first combined input.
  This is Yuguang Chen's **standalone, single-night** reduction — solid proof the Dec-29 frames
  were successfully turned into a cube.

**What's NOT solid:**
- The actual headline cube (`raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits`)
  has a `HISTORY`/`COMMENT` block that is byte-identical to a **single-frame** kcwidrp reduction
  log of `kr251117_00129` (the first Nov-17 frame), with exactly one line appended: "Reduced and
  mosaicked by Kaustubh Rajesh Gupta... using KCWIKit v0.3Dev, KSkyWizard v0.2.2025-06-04." No
  NCOMBINE keyword, no input list, no per-night breakdown anywhere. **The KSkyWizard mosaic step
  does not preserve per-input provenance in its output header at all.**
- Checked Kaustubh's actual Drive working folder (`.../Kaustubh_Nov_17_2025_stack/`) — only output
  cubes + a PNG, no `.list`/`.par`/manifest files. His intermediate `desj0206-stack-{1,2,3}.fits`
  files — the literal 3 inputs `desj0206-red.py` expects — are **not discoverable anywhere on the
  Drive**, only (presumably) on his local machine.
- The only place "3 stacks combined" is textually confirmed is the header COMMENT of the
  **alternate, non-headline** cube `DESJ0206-red_medium_combined.fits` ("Files used:
  desj0206/desj0206-stack-1.fits, -2.fits, -3.fits") — but that combine step also strips what raw
  frames are inside each opaque stack.

**Conclusion:** "3-night, 3-program combine" is a structurally plausible inference (3 nights, 3
distinct RL central wavelengths — 7400/7150/8000 — each requiring its own per-config rectification
before mosaicking) but is **not independently verified at the file level**, for either the
alternate or the headline cube. Nov-17 has the strongest evidence (its first frame's HISTORY
literally appears in the headline cube). Aug-30 and Dec-29 inclusion rest on inference alone.

**Side finding (not blocking):** the local `dec29_2024_rl8000.par` alignment config references an
unrelated calibration field ("VS44," RA 114.28°/Dec 65.61° — nowhere near AGEL0206) as its WCS
cross-correlation anchor. Doesn't undermine the raw-frame identity (independently nailed down
above), but flags that this specific local file shouldn't be trusted as supporting evidence —
likely a stale/reused template in Yuguang's pipeline directory.

## Recommended next step (not done this session)

**Ask Kaustubh Gupta** (mosaicked the headline cube) **or Yuguang Chen** (reduced the Dec-29
single-night cube) directly which raw inputs went into the final `_mtwdo_` combine. This is the
only way to close the gap — no further local/Drive file inspection can resolve it, since both
combine steps (kcwiRedux `step2()` and KSkyWizard mosaic) strip per-input provenance from their
output headers.

## Files touched this session

- `DRAFTING_FACTS_paper_2026-05-29.md`, `METHODS_AND_SYSTEMATICS.md`, `TESTS_AND_DIAGNOSTICS.md` —
  replaced "K409 PI TBD" with "K409 PI = Alcorn (Airtable, 2026-07-06)" in all instances.
- `~/.claude/.../memory/reference_kcwi_data_properties.md` — added K409 PI resolution, the
  Dec-29 re-verification (item 4), and the new critical "folded into headline?" caveat (item 5).
- `~/.claude/.../memory/MEMORY.md` — updated the `reference_kcwi_data_properties.md` index line
  to flag the open verification gap.
- This file.

**Deliberately NOT touched:** the "3-night, 3-program combine" framing in `CLAUDE.md` and the main
prose of `DRAFTING_FACTS` — that claim is now known to be weaker than currently written, but
downgrading paper language is a decision for the user, not something to silently patch. Flag this
explicitly next session if Kaustubh/Yuguang haven't been asked yet.
