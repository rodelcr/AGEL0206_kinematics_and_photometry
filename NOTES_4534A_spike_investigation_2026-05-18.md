# Investigation: def-rest 4534 Å spike — telluric residual at O2 A-band edge

**Date:** 2026-05-18
**Scope:** read-only diagnostic of the empirically-masked def-rest 4525–4545 Å band (TESTS_AND_DIAGNOSTICS.md row J2, PROJECT_BRIEF.md "unidentified" line). No notebooks, figures, or production code were modified.

---

## TL;DR

**Verdict: telluric/sky-subtraction residual at the leading edge of the O2 A-band (7593–7700 Å).** The spike is NOT a source-emission line and is NOT localized to the deflector — it appears at the same observed wavelength, with comparable per-spaxel amplitude, in off-deflector spaxels that contain ZERO continuum from the deflector galaxy. The current mask (def-rest 4525–4545 Å) is the correct action; only the figure-caption/test-catalog label needs to change from "unidentified, ~7σ excess" to "telluric residual at O2 A-band edge (obs 7593–7700 Å)."

Diagnostic figure: `results/figures/spike_4534A_diagnostic.png`

---

## What I ran

Standalone script `/tmp/spike_4534_diagnostic2.py` (not committed). Operating on `Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits` (3317 × 100 × 100; CRVAL3=5625, Δλ=1.0 Å). Three per-spaxel-mean spectra over obs 7560–7700 Å, plus a spatial map of (flux integrated over obs 7593–7601 Å) − (local off-line continuum).

Regions (intersected with the valid-FOV mask `white > 0`):

| region | definition | N valid spaxels |
|---|---|---|
| (a) Deflector | R < 1.5″ around (cy,cx) = (50,50) | 69 |
| (b) Sky/noise patch | CLAUDE.md cube[:, 28:40, 45:70] | 300 |
| (c) Off-deflector | annulus 4″ < R < 8″ from (50,50) | 1680 |

Note: I replaced the originally-suggested corner cube[:, 5:15, 5:15] with an annulus 4–8″, because the four corners of this KCWI cube are zero-padded (only 3888 of 10000 spaxels have any data). The 4–8″ annulus is the right "in-FOV, off-deflector" comparison sample.

## What I found

Spike measurement (continuum from two ±15 Å sidebands at 7570–7585 and 7610–7625; FWHM from half-max pixels within ±10 Å):

| region | continuum | peak λ | excess (peak−cont) | S/N | FWHM |
|---|---|---|---|---|---|
| (a) Deflector R<1.5″ | 0.0110 | 7601.0 Å | +0.0170 | 12.3 | 3.0 Å |
| (b) Sky box ∩ FOV | 0.0002 | 7599.0 Å | +0.0088 | 8.0 | 4.0 Å |
| (c) Off-deflector 4–8″ | 0.0010 | 7599.0 Å | +0.0079 | 9.6 | 5.0 Å |

(units: per-spaxel mean, 10⁻⁸ erg/s/cm³/arcsec² — the cube's native BUNIT.)

Key facts:

1. **Same observed wavelength in all three regions.** The peak is at obs 7599–7601 Å in every aperture. Def-rest at z=0.67564: 4535–4536 Å, well inside the masked 4525–4545 Å band. The spike position does not depend on which spaxels you stack.

2. **The spike is NOT scaling with deflector continuum.**
   - White-light continuum (per spaxel): deflector/off-annulus = **17.5×**
   - 7600-Å spike amplitude (integrated, per spaxel): deflector/off-annulus = **1.11×**

   A real source-emission line at this rest wavelength (deflector or z=1.302 source) would scale with the source's spatial brightness; the ratio is essentially flat. This is the cleanest possible signature of a wavelength-dependent additive systematic (sky/telluric residual) that is approximately uniform across the FOV.

3. **Pearson r(spike map, white-light map) = +0.089.** Statistically no correlation — the spike map is approximately spatially uniform; the white-light map is strongly peaked on the deflector. See the right-column panels of the diagnostic figure.

4. **The wavelength match to a known telluric is exact.** The O2 atmospheric A-band (b¹Σg⁺ ← X³Σg⁻, 0–0 vibrational) absorbs from 7593 Å to ~7700 Å, with a sharp blue edge that goes from zero to deep absorption in a few Å around 7595–7600 Å. Any small mismatch between the telluric correction template and the actual airmass at the time of each KCWI exposure produces a positive residual at exactly this edge. The KCWI DRP's sky-subtraction step also leaves residuals at the strongest sky-emission features (and there is OH airglow nearby in 7720–7920 Å).

5. **FWHM ~3–5 Å is consistent with a sharp edge artifact, not a redshifted ISM emission line.** A z=1.302 narrow ISM line (e.g., [CII]λ2326 hits obs 5354 Å, not 7600 Å — and in general no plausible source-rest narrow emission line maps to obs 7596 Å for either z=0.67564 or z=1.302).

6. **Stacking does not average it out.** The deflector spectrum is built by I-weighted co-add of ~69 spaxels (R<1.5″ for this diagnostic; ~50 for the production R<R_e ppxf aperture). Because the telluric residual is wavelength-locked in the observed frame and approximately uniform across the FOV, weighted averaging preserves it — exactly as we see.

## Provenance of the original empirical detection

The mask was placed in `scripts/run_window_sweep.py:ARC_MASK_DEF_REST` as `(4525.0, 4545.0)` with the comment "user-flagged ~4550 (source rest 3300 Å)". TESTS_AND_DIAGNOSTICS.md row J2 currently records it as a "wild-bootstrap >7σ excess" with no line identification. The source-rest-3300-Å hypothesis (Balmer-edge continuum break in the z=1.302 source) was never confirmed and is now ruled out by the test above — the spike does not scale with source-arc surface brightness either (the arc lies off the deflector center and the spike map shows no arc-shaped excess).

## Recommendation

1. **Keep the mask as-is.** Def-rest 4525–4545 Å (obs 7593–7626 Å) covers the entire sharp blue edge of the O2 A-band. Narrowing it would risk leaking residual into the ppxf fit; widening it is unnecessary and would cost good pixels.

2. **Relabel in the catalog and figure caption.** Change "def-rest 4525–4545 Å = empirically identified excess; source rest" (TESTS_AND_DIAGNOSTICS.md §2/§3 around row J2) and PROJECT_BRIEF.md "unidentified" to:
   > def-rest 4525–4545 Å (obs 7593–7626 Å) — telluric residual at the leading edge of the atmospheric O₂ A-band; verified spatially uniform across the FOV (NOTES_4534A_spike_investigation_2026-05-18.md).

3. **Add the comment in `scripts/run_window_sweep.py:ARC_MASK_DEF_REST`.** Replace the inline "user-flagged ~4550 (source rest 3300 Å)" with "O2 A-band leading edge at obs 7593–7626 Å; see NOTES_4534A_spike_investigation_2026-05-18.md."

4. **Consider widening to def-rest 4525–4570 Å (obs 7593–7669 Å)** as a future test if any band-edge residual is later seen at the red side of the A-band. Currently the deflector ppxf residuals are clean on the red side (peak shape decays back to continuum by obs 7605 Å), so this is precautionary only and not required for the headline.

No paper numbers change. The headline σ_e and the wide-window cross-check are produced WITH this mask applied, so the masked-band reinterpretation is purely a label/provenance fix.

## Files

- Diagnostic figure: `/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry/results/figures/spike_4534A_diagnostic.png`
- Diagnostic script (uncommitted, in /tmp): `/tmp/spike_4534_diagnostic2.py`
- Cube: `/Users/rosador/Documents/AGEL/velocity_dispersion_from_IFU/Nov17_2025_DESJ0206_RL_combined_icubes_wcs.fits`
