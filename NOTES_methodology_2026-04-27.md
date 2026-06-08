<!-- pv-skip-file: dated snapshot — records historical (now-superseded) numbers; current headline in CLAUDE.md / results/PAPER_VALUES.json -->

# Methodology session notes — 2026-04-27

Post-compaction working session focused on documenting the σ_e methodology choices. No new measurements made except a launched no-mask sensitivity test (running in background as of session end).

---

## 1. Cube provenance — corrected understanding

**Earlier error:** I had recorded that the cube used by the σ_e pipeline (`Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`, shape 3317×100×100) was the WRONG cube for the paper, and that we needed to re-run on the multi-night `DESJ0206-red_medium_combined.fits` (469 MB, 4269×120×120).

**Correction (per user):** The cube we used IS the final, most-reduced data product. Only its FITS header metadata is mislabeled (`DATE-OBS=2025-11-17`, `PROGNAME=U002`, `PROGPI=Jones` describes only the last input frame `kr251117_00129`, not the full input set). The headline σ_e numbers stand and do not need re-running.

**Action taken:** Softened the "WRONG cube" framing in `HANDOVER.md` and `reference_kcwi_data_properties.md` to "metadata mislabeled — cube content is the final reduction; σ_e numbers stand." The 469 MB alternate combine in `raw_KCWI/` was demoted from "must re-run" to "optional sanity-check."

### Multi-night provenance — UPDATED 2026-05-26 from the actual Keck observing logs

The earlier paper-cite recommendation here ("Aug 29 K409 + Sept 29 K409 + Dec 29 + Nov 17") was WRONG in two ways: (a) Sept 29 K409 has **zero DESJ0206 entries** in its observing log (verified against `provenance/K409_2025-09-29_UTC.html`) — that night did not contribute; (b) Dec 29 2024 is unconfirmed pending a raw FITS header dump. The verified contributing nights are:

| Night (UT) | Program | PI | DESJ0206 frames | Source / log file |
|---|---|---|---|---|
| **2025 Aug 30** | K409 | TBD (not in log header) | 12 RED `kr250830_00090–00101` + 4 BLUE `kb250830_00052–00055` (60 min RED + 66 min BLUE) | local `raw_KCWI/provenance/K409_2025-08-30_UTC.html` ↔ Drive id `1yoH_elZqEsFHPRc6f-XCB98l-7dGJewN` |
| **2025 Nov 17** | U002 | Jones | 20 RED `kr251117_00129–00148` + 5 BLUE `kb251117_00087–00091` (90 min RED + 110 min BLUE) | Drive `NightLog_KCWI_2025-11-17.pdf` id `1HVNoH2CSAc214b_PieBQ_UkeEasVl6q4` |
| 2024 Dec 29 | TBD | (Yuguang Chen reduction) | 4 RED `kr241229_00092–00095` per local `dec29_2024_rl8000.list` — **NOT independently confirmed as DESJ0206 pointings**; no machine-readable Drive log (only 44 MB image-only PDF id `16zWO8yaCIvRtoCGenNTndvB7dYdrCkMH`) | (raw frames on Yuguang Chen's machine at `~/obs/2024dec29/`; alignment QA PNGs in Drive id `1ELcTQiXV9m04Y8sA9st3OlpsNk7LE0D7`) |

Observers on both verified nights include **the user (Cordova)** along with Glazebrook, Jones, Kacprzak, Tran, Vasan G C, Alcorn; Nov 17 also includes Gupta (= Kaustubh, the data reducer), Rhoades, Chen, Barone. RL grating central wavelengths differ between nights (Aug 30: 7400 Å; Nov 17: 7150 Å); both rebinned onto a common 1.0 Å/pix grid spanning 5625–8941 Å by the kcwiRedux pipeline.

**For the paper:** cite the multi-night provenance as **K409 (Aug 30 UT) + U002 (Nov 17 UT)** with the Dec 29 2024 contribution **explicitly flagged as pending verification**. Do NOT cite Sept 29 K409 — it had no DESJ0206 entries. See `reference_kcwi_data_properties.md` (updated 2026-05-26) for the full inventory and Drive IDs.

---

## 2. Updated Figure 2 wiring (figures.ipynb)

Added a new cell under `# Updated Sigma measurement plot` in `/Users/rosador/Documents/AGEL/AGEL_0206_ApJL_Figures/figures.ipynb`. This rewires Figure 2 to use the §6cum cumulative I-weighted aperture spectrum at R<R_e (from `results/annular_bootstrap_07c/cumR_2p305_{sps}.npz`) and the joint combined-SPS posterior `cum_combined_cumR_2p305_samples` from `results/sigma_e_radial_07c.npz`.

The σ_e inset now reads **σ_e(<R_e) = 267 −24/+24 km/s** (N=22,500 = 3 SPS × 15 degrees × 500 boot). Spectrum panel shows the median ppxf model across all 45 (3 × 15) fits, residuals, and the same arc-mask exclusion shading.

Output PDF: `AGEL0206_sigma_e_SED_final.pdf`. Right-panel SED + filter curves preserved unchanged from the earlier cell.

---

## 3. Gültekin formula attribution — 2D vs 1D

Caught: the formula `σ_e² = ∫ I(r)[V²+σ²] 2πr dr / ∫ I(r) 2πr dr` was written in the memory file as "Gültekin 2009 eq. 1 / KH13 §6.6". User pointed out the 2πr factor isn't in the original Gültekin paper.

**The truth:** Gültekin+2009 eq. 1 is a **1D longslit** integral (no 2πr). Our IFU implementation generalizes it to 2D — for axisymmetric I(r), dA = 2πr dr in the 2D integral, and the 2π cancels in the ratio (the r factor is what's actually doing physical work). This 2D extension is the de facto IFU standard:

- Cappellari et al. 2006 (SAURON IV) eq. 1: σ_e² = Σ_n F_n (V_n² + σ_n²) / Σ_n F_n on Voronoi-binned spaxels
- Cappellari et al. 2013 (ATLAS3D paper III): same
- Westfall et al. 2019 (MaNGA DAP): same

**Action taken:** Updated `reference_gultekin_implementation.md` to clarify Gültekin's 1D-longslit form vs our 2D-IFU extension, and to cite SAURON/ATLAS3D/MaNGA. Removed misleading "Eq. 1" attribution from the paper-ready language.

---

## 4. §6cum vs §7 — judgment that §6cum is best

Wrote down seven specific reasons §6cum (cumulative I-weighted aperture ppxf) is the right primary estimator (now in `reference_cumulative_vs_annular_sigma_e.md`):

1. It IS what the M•–σ literature actually computes (KH13, Greene+20, SAURON, ATLAS3D, MaNGA).
2. No binning to defend (no equal-r vs equal-N debate at the measurement step).
3. Single LOSVD fit preserves line-shape information that §7's moment-pooling discards.
4. Bright-center I-weighting auto-suppresses arc contamination (verified by nb07e: < 0.1 km/s shift).
5. Fewer ppxf calls → smaller compounded noise budget; no discretization error.
6. No V_sys anchoring choice (§7 needs one; can shift σ_e by 5–15 km/s).
7. Robust to non-axisymmetry / mask gaps.

**§6cum is the measurement; §7 is the diagnostics.** §7 is for the σ(R) profile, per-annulus systematics inspection, and cross-check on the cumulative number. Never quote a §7 number as the headline.

---

## 5. Equal-r vs equal-N binning — judgment

For AGEL0206's geometry (Sersic-like I(r) + lensed arc at r ~ 1″), **equal-N (equal-spaxel-count) binning inside R_safe = 3R_e/4 + 1 outer flagged bin (nb07c) is the right choice for §7**, with equal-width binning (nb07) kept as a sensitivity test only.

- Equal-r is the *worst* combination for AGEL0206: noisy inner shell (few spaxels) + biased outer ring (most arc-contaminated, dominates the F-weighted sum).
- Equal-N gives balanced bootstrap variance per bin, narrow outer bins that cleanly diagnose arc contamination, and matches §6cum's underlying I-weighting logic in the N→large limit.
- The 257 (equal-width) vs 267 (equal-N) shift on σ_e(<R_e) is at the SPS-systematic level (±24 km/s).

**Don't quote equal-width §7 as headline.** Either use §6cum (no binning), or use equal-N §7 with the outer arc-bin flagged.

Action: documented across `reference_cumulative_vs_annular_sigma_e.md`, `project_sigma_e_gultekin.md`, `reference_gultekin_implementation.md`, and `CLAUDE.md`.

---

## 6. §6cum mechanics — clarification

User asked: does §6cum measure σ per spaxel, or build one spectrum and fit it?

**Answer:** Builds ONE I-weighted aperture spectrum (`flux_agg = Σ_spax (I_spax/ΣI) × cube[:, spax]`), then runs ppxf ONCE on it. Returns one (V, σ) per aperture per (SPS, degree). Per-spaxel kinematics are never measured.

Mathematically the result is the second moment of the I-weighted mixture of per-spaxel LOSVDs — which IS the Gültekin σ_e definition in continuous form. ppxf approximates this mixture with a single Gauss-Hermite LOSVD (fine for AGEL0206's low V/σ; would need attention for a strongly rotating system).

The per-spaxel weight is **I(spaxel)** (surface brightness), NOT F_j (per-annulus integrated flux — that's §7's weight). In the limit of N annuli → number-of-spaxels with proper intra-annulus I-weighting AND a joint LOSVD fit instead of post-fit moment pooling, §7 → §6cum exactly.

---

## 7. F200LP arc-mask sensitivity test — RESULT (run completed 2026-04-27 19:48)

**Question raised:** did we directly verify the F200LP arc mask is not a dominant systematic on σ_e?

**Existing evidence (orthogonal but not direct):**
- nb07d (mask dilation): more aggressive mask → σ inflates (over-masking inflates noise). Rules out "F200 mask too lenient".
- nb07e (arc-spectrum subtraction, N=500): subtracting α_j × arc shifts σ_e by < 0.1 km/s. Rules out "leftover arc light through F200 mask matters".
- nb07a (6 I-maps with different mask handling): σ_e spread ≤ 3–8% across mask-zero vs Sersic-fill variants.

**Direct test we hadn't done:** §6cum at N=500 with F200 mask DISABLED entirely, comparing to the masked headline.

**Action taken:** Added 3 cells (markdown header + extraction/bootstrap + comparison) to `notebooks/07c_sigma_e_equalN.ipynb` after the existing §6cum comparison plot (cells 34-36 in the new layout). Sets `arc_spax_mask_OFF = np.zeros_like(arc_spax_mask)`, re-runs the full cumulative pipeline at N=500 (parallelized with n_jobs=8), saves to `results/annular_bootstrap_07c_nomask/` (separate from production cache), prints/plots the masked-vs-nomask comparison per aperture.

**Result:** σ_e(<R_e) = 250.96 ± 23 km/s (no mask) vs 267.32 ± 24 km/s (F200LP masked) → **Δσ_e = −16.36 km/s**.

| R_max | masked | nomask | Δ |
|---|---|---|---|
| 0.762″ | 209.18 | 206.92 | −2.26 |
| 1.099″ | 225.77 | 221.18 | −4.59 |
| 1.321″ | 229.32 | 224.05 | −5.27 |
| 1.546″ | 235.81 | 227.04 | −8.77 |
| 1.724″ | 238.78 | 231.64 | −7.14 |
| **2.305″ (R_e)** | **267.32** | **250.96** | **−16.36** |

Per-SPS at R_e: FSPS −10, EMILES −20, XSL −17 km/s. Direction is consistent (no-mask < masked) and the shift grows monotonically with R_max — pattern of arc dilution (admitting flat rest-UV continuum from z=1.302 dilutes the deflector's Ca H+K + G-band absorption depth, ppxf reads slightly lower σ).

**Verdict:** prediction was wrong. The shift is ~16 km/s at R_e, ~68% of the 24 km/s SPS-systematic budget — NOT sub-dominant. The F200 mask is doing real work.

**Reconciliation with nb07e's "<0.1 km/s":** nb07e tested *residual arc light leaking through the F200 mask* — minimal (the mask is doing its job inside the masked region). This no-mask test removes the F200 mask entirely (admits ~30 arc-bright spaxels at R<R_e). Different tests; both correct.

**Implication for paper:** headline σ_e(<R_e) = 267 ± 24 (F200LP masked) stands. Add a footnote/sensitivity-table line: "Disabling the F200LP arc mask shifts σ_e(<R_e) to 251 ± 23 km/s (a −16 km/s arc-dilution shift), still consistent with the masked headline at < 1σ; the masked value is preferred because the lensed arc carries flat rest-UV continuum that dilutes the deflector's stellar absorption features."

Cache: `results/annular_bootstrap_07c_nomask/`. Plot: `results/figures/nb07c_s6cum_nomask_sensitivity.png`.

---

## Files touched this session

- `HANDOVER.md` — softened "wrong cube" framing
- `CLAUDE.md` — refreshed Method choice section
- `~/.claude/.../memory/reference_kcwi_data_properties.md` — corrected provenance framing
- `~/.claude/.../memory/reference_cumulative_vs_annular_sigma_e.md` — added §6cum-is-best section + binning judgment
- `~/.claude/.../memory/reference_gultekin_implementation.md` — fixed 2D/1D attribution + headline numbers update
- `~/.claude/.../memory/project_sigma_e_gultekin.md` — refreshed cumulative-vs-annular short version + binning paragraph
- `notebooks/07c_sigma_e_equalN.ipynb` — added §6cum-nomask sensitivity section (executing in background)
- `/Users/rosador/Documents/AGEL/AGEL_0206_ApJL_Figures/figures.ipynb` — added updated Figure 2 cell under "# Updated Sigma measurement plot"
