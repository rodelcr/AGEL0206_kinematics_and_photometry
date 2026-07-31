# Is AGEL J0206's deflector a CORED or coreless elliptical? — report (2026-06-17)

**Short answer: it is fundamentally UNRESOLVABLE with HST/JWST imaging at this redshift.** The
depleted core a galaxy this massive is *expected* to have is a few–40 milliarcsec across — at or below
the sharpest PSF we can build (NIRCam SW HWHM ≈ 24 mas). We can neither confirm nor rule out a core; but
it does not matter for the science, because a depleted core removes <1% of the light (negligible for
M⋆). What the data *can* see in the centre is not a core but the imperfect single-Sérsic description of
the **outer cD envelope** (Test A) — which is already carried in the M⋆ systematic budget.

Data products: `notebooks/18_cored_vs_coreless_sersic.ipynb`; `scripts/{psf_star_census,build_psf_models,
core_sersic_test}.py`; `results/{psf_models/,core_sersic_test.npz}`; figures
`results/figures/{cored_test_A_boxsweep,cored_test_B_profile,core_sersic_test,cog_total_light}.png`.
TESTS rows A6 (PSF), A7 (core-Sérsic).

---

## 1. Why we expected a core

By the core/coreless dichotomy (Faber+1997, Lauer+2007, Kormendy+2009), massive slow-rotator
ellipticals above σ ≈ 240 km/s almost always host a central surface-brightness **deficit** ("depleted
core") scoured by a binary SMBH during dry mergers. AGEL J0206's deflector (σ_e = 267 km/s, log M⋆ =
11.5, passive, ~100 kpc from cluster ACT-CL J0206) is squarely in that regime → a core is *physically
expected*.

## 2. The decisive obstacle: angular size vs resolution

The break/core radius r_b correlates with M• and with galaxy mass. For depleted-core ellipticals,
r_b ranges from a few × 10 pc up to ~0.5 kpc, reaching ~1–4 kpc **only** in the most extreme
BCGs (M• ~ 10¹⁰; e.g. A2261-BCG). For our M• ≈ 5×10⁸ M⊙ the expected core is **r_b ~ 0.02–0.1 kpc**
(tens of pc to ~100 pc; Lauer+2007, Rusli+2013, Thomas+2016) — and even a generous 0.3 kpc is an upper
edge. At z = 0.67564 the scale is **7.276 kpc/arcsec**, so:

| quantity | physical | angular | vs sharpest PSF (F150W2, FWHM 0.049″ → HWHM 0.024″) |
|----------|----------|---------|----------------------------------------------------|
| expected core (typical) | 0.02–0.1 kpc | **0.003–0.014″ (3–14 mas)** | **5–8× too small to resolve** |
| expected core (generous upper) | 0.3 kpc | 0.041″ | ~marginal, ~1.7× HWHM |
| **resolution floor** to *detect* a core (r_b ≳ HWHM) | 0.18 kpc | **≈ 0.024″ (24 mas)** | — |
| PSF FWHM per band | — | F150W2 0.049″, F200LP 0.085″, F140W 0.103″, F322W2 0.112″ | — |

The expected core is **smaller than the finest PSF half-width by ~2–8×** (it would have to be an
extreme-BCG-scale ~0.2–0.3 kpc core to even approach detectability). To resolve the typical case you
would need a PSF FWHM of **≲ 5–10 mas** — roughly 5–10× sharper than NIRCam SW, i.e. beyond HST/JWST
imaging entirely (the realm of a 30 m-class ELT with AO, or interferometry). **So: not a tooling gap —
a hard diffraction/distance limit.** Building real PSFs (A6) confirmed this rather than removing it.

## 3. What the PSF-convolved fit actually found (A7) — and why it is NOT a core

We built real PSFs (A6: F140W from the STScI PSFSTD library; F200LP/F150W2/F322W2 empirical EPSFs from
full-frame field stars) and ran PSF-convolved single-Sérsic vs core-Sérsic (Trujillo+2004) fits on the
two deepest/sharpest bands (`scripts/core_sersic_test.py`):

| band | PSF FWHM | single-Sérsic n | core-Sérsic r_b | (kpc) | γ | ΔBIC (single−core) | inner residual deepens toward r=0? |
|------|----------|-----------------|------------------|-------|---|--------------------|-----------------------------------|
| F150W2 | 0.049″ | 0.92 | 0.50″ | 3.6 | 0.26 | +386 | **no** (flat ~−6%, oscillating) |
| F140W  | 0.103″ | 1.08 | 0.34″ | 2.5 | 0.17 | +43  | **no** |

A free core-Sérsic "prefers" the fit by large ΔBIC — **but this is a fitting artifact, not a core
detection:**
1. the fitted r_b = 2.5–3.6 kpc is **~10–50× larger** than any depleted core for this M• and ~10–20× the
   PSF HWHM — it is a galaxy-scale break, not a nucleus;
2. the single-Sérsic inner residual does **not deepen toward the centre** (the signature of a real
   core); it is flat (~−6%) and part of a global ±5–15% **oscillation**;
3. the single-Sérsic index rails low (n = 0.9–1.1 vs the validated 1.2–1.6).

Together these say the core-Sérsic's inner flattening is **absorbing the global single-Sérsic
mismatch** — i.e. the extended high-n / cD **outer envelope** that Test A already flagged as
box-limited — not a central deficit.

> **Methodological lesson (logged):** a free core-Sérsic ΔBIC "preference" is *not* a core detection
> when the baseline single-Sérsic mis-fits the galaxy globally. A real-core claim requires r_b that is
> both *resolved* (> PSF HWHM) **and** *physically small* (≲ 0.3 kpc), **and** an inner residual that
> deepens toward r=0. None hold here.

## 4. Verdict and impact

- **Cored or coreless?** *Indeterminate.* A genuine depleted core is unresolved at every band; the data
  cannot distinguish the two. **Upper limit: r_b < 0.024″ = 0.17 kpc** (F150W2 HWHM). A core, if present
  as expected, lives below this.
- **Does it bias the photometry / M⋆?** **No.** A depleted core removes <1% of the total light → ΔM⋆ ≪
  0.01 dex. The single-Sérsic is the appropriate *inner* light model.
- **Where the real structure is:** the *outer* envelope. Test A shows n and total light keep rising with
  fit box (non-convergent cD/ICL envelope); Test A2's curve-of-growth is +0.2 dex above the
  single-Sérsic total. That ambiguity is the one that matters for M⋆ and is **already captured** by the
  aperture-correction-model systematic (it spans the more-outer-light direction; adding a separate term
  would double-count — see §3.2 of DRAFTING_FACTS and notebook 18 synthesis).
- **M•–core cross-check:** the would-be "core scouring ↔ M•" consistency test (mass deficit vs M•) is
  **not feasible** here — it needs a resolved core.

**Bottom line:** this is a distance/diffraction limit, full stop. The investigation is complete and
correctly concluded; no further inner-profile work is warranted with HST/JWST. The only way forward
would be ~mas-resolution imaging (ELT-class), which is out of scope — and even then the science payoff
(a <1%-light core) does not justify it. The cored/coreless question is therefore **closed as
unresolvable**, with the scientifically relevant outer-envelope effect already folded into the M⋆ budget.
