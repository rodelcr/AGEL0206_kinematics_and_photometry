"""Deterministic paper-values registry — single source of truth (Task 8, 2026-06-08).

Computes EVERY headline number from the result caches on disk, with explicit
provenance (which cache + which formula), and emits a machine-readable
`results/PAPER_VALUES.json`. The figures notebooks and the summary
docs/PDFs should read FROM this registry (or its JSON) rather than hard-coding
literals — so no value is ever LLM-typed from local context.

Design
──────
- Every entry is a dict: {value, [err_lo, err_hi, err_sym], unit, provenance,
  formula}. `provenance` names the cache file; `formula` names the reduction.
- "live" components are recomputed from a cache here; "carried" components are
  M11 sweep results stored as documented constants with a provenance pointer to
  their own cache (wire them to live reads as a later refinement — flagged).
- `python scripts/paper_values.py`         → recompute + write JSON + print table
- `python scripts/paper_values.py --check <file>...` → scan files for headline
  numbers and flag any that disagree with the registry (drift linter).
- `python scripts/paper_values.py --render <file>...` → regenerate the marked
  blocks `<!-- PV:auto:NAME -->…<!-- /PV:auto:NAME -->` from the registry, so the
  numbers inside are GENERATED (never hand-typed). Blocks: `headline` (the
  headline-numbers list) and `budget` (the σ_e systematic-budget table).
  Idempotent: re-rendering an up-to-date file is a no-op.

Run from the repo root (the script chdir's there).
"""
from __future__ import annotations
import os, sys, json, argparse, glob, re
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); os.chdir(REPO)
RESULTS = REPO / 'results'


def _entry(value, unit, provenance, formula, **errs):
    d = {'value': float(value), 'unit': unit, 'provenance': provenance, 'formula': formula}
    d.update({k: float(v) for k, v in errs.items()})
    return d


# ─────────────────────────────────────────────────────────────────────────────
def sigma_e_pool():
    """Pool the 3-SPS wild-bootstrap σ samples from the headline R_e=2.097" best-mask
    aperture (resys_best_mean) → p16/p50/p84 (headline σ_e central + stat error).
    The aperture is r<R_e with R_e=2.097" (best-mask CoG, adopted 2026-06-11 'best
    mask throughout'); identical pipeline to new_clean_hei but at the new headline R_e."""
    caches = sorted(glob.glob(str(RESULTS / 'run_wide_sigma_e' / 'resys_best_mean'
                                  / 'wR3800_5400_arcmask_*_N500.npz')))
    caches = [c for c in caches if 'figure2_data' not in c]
    assert len(caches) == 3, f'expected 3 SPS caches, found {len(caches)}: {caches}'
    pool, per_sps = [], {}
    for c in caches:
        d = np.load(c, allow_pickle=True)
        s = d['sig_boot'].ravel(); s = s[np.isfinite(s)]
        pool.append(s)
        sps = str(d['sps']); per_sps[sps] = float(np.percentile(s, 50))
    pool = np.concatenate(pool)
    p16, p50, p84 = (float(x) for x in np.percentile(pool, [16, 50, 84]))
    prov = 'results/run_wide_sigma_e/resys_best_mean/wR3800_5400_arcmask_{fsps,emiles,xsl}_*_N500.npz'
    return p16, p50, p84, per_sps, prov


def _quad(*xs):
    return float(np.sqrt(sum(x * x for x in xs)))


def build():
    reg = {'_meta': {'generated_by': 'scripts/paper_values.py',
                     'note': 'Single source of truth. Do NOT hand-edit headline numbers elsewhere; '
                             'read from here. Regenerate with: python scripts/paper_values.py'}}

    # ── σ_e central + stat ───────────────────────────────────────────────────
    p16, p50, p84, per_sps, prov = sigma_e_pool()
    stat_lo, stat_hi = p50 - p16, p84 - p50
    stat_sym = (p84 - p16) / 2
    reg['sigma_e'] = {
        'central': _entry(p50, 'km/s', prov, 'p50 of pooled 3-SPS×15-deg×N500 sig_boot'),
        'stat_lo': _entry(stat_lo, 'km/s', prov, 'p50 - p16'),
        'stat_hi': _entry(stat_hi, 'km/s', prov, 'p84 - p50'),
        'stat_sym': _entry(stat_sym, 'km/s', prov, '(p84 - p16)/2'),
        'per_sps_p50': {k: round(v, 2) for k, v in per_sps.items()},
    }

    # ── σ_e systematic budget components ─────────────────────────────────────
    # cache-derived (live):
    d7 = np.load(RESULTS / 'sigma_e_Re_grid_N500.npz', allow_pickle=True)
    # ADOPTED (user 2026-06-11): best-mask light family {best F140W/mean/F200LP} only.
    # Full 7-pt grid (incl. CaHK+G 2.90") = ±9.98 is the conservative ceiling; CaHK+G
    # is kept as a noted cross-check (different I-map definition), NOT folded.
    re_source = float(d7['sys_re_bestlight'])
    msys = np.load(RESULTS / 'sigma_e_mask_systematic_N500.npz', allow_pickle=True)
    mask_approach = float(msys['sys_mask'])  # cross-check, NOT a separate budget line

    # named components (carried M11 constants pointed at their caches + live R_e-source)
    components = {
        'Ishape':    _entry(2.29, 'km/s', 'results/ishape_arcfree_pooled.npz (M11 + arc-free PSF-matched, 2026-06-12)',
                            '14-shape peak-to-peak/2 (266.79-271.37): the 10 raw maps + 4 arc-free '
                            'PSF-matched (Sérsic/filled→1.27" conv); was ±2.27 over 10', carried=1),
        # General arc-masking systematic (formerly "F200 mask"): the measured arc-dilution
        # sensitivity from the F200LP-located mask-weight sweep w∈{0..1}; subsumes the
        # 4-approach arc-mask-definition spread (±4.58, larger-of-two, no double-count).
        'F200mask':  _entry(6.65, 'km/s', 'results/maskweight_sweep_wR3800_5400_arcmask/ (M11)',
                            'arc masking: (w00-w100)/2 weight sweep; subsumes mask-approach 4.58 (larger-of-two)', carried=1),
        # NOTE (2026-06-14): the old "frame (vac/air) ±5" term was DROPPED. It is not a live
        # uncertainty: (1) applying each SPS in its native vac/air frame is a DETERMINISTIC
        # correction whose actual σ impact is <0.5 km/s (TESTS D4); (2) the ±5 was defined as
        # the "max σ shift across SPS" (TESTS C3) = the inter-SPS spread, which the pooled
        # 3-SPS bootstrap stat ALREADY marginalizes (between-SPS ±2.04 ⊂ pooled stat ±4.6).
        # Carrying it separately double-counted the pool. User decision 2026-06-14.
        'centering': _entry(4.0, 'km/s', 'HST WCS centroid', 'carried constant', carried=1),
        'fitwindow': _entry(3.82, 'km/s', 'results/run_wide_sigma_e/ 3-window (M11)',
                            'peak-to-peak/2 across 3 fit windows', carried=1),
        'reduction': _entry(3.45, 'km/s', 'results/run_wide_sigma_e/{new,headline}_clean_hei (M10)',
                            'half-Δ between NEW and OLD reductions', carried=1),
        'Re_source': _entry(re_source, 'km/s', 'results/sigma_e_Re_grid_N500.npz',
                            'peak-to-peak/2 across best-mask CoG family {F140W 1.91", mean 2.10", '
                            'F200LP 2.28"} at the 2.097" headline (user 2026-06-11; CaHK+G & full grid '
                            'are cross-checks, not folded)'),
    }
    sys_quad = _quad(*[c['value'] for c in components.values()])
    reg['sigma_e_budget'] = {
        'components': components,
        'sys_quad': _entry(sys_quad, 'km/s', 'quadrature of components above', 'sqrt(Σ component²)'),
        'total_lo': _entry(_quad(stat_lo, sys_quad), 'km/s', 'stat_lo ⊕ sys_quad', 'sqrt(stat_lo²+sys²)'),
        'total_hi': _entry(_quad(stat_hi, sys_quad), 'km/s', 'stat_hi ⊕ sys_quad', 'sqrt(stat_hi²+sys²)'),
        'total_sym': _entry(_quad(stat_sym, sys_quad), 'km/s', 'stat_sym ⊕ sys_quad', 'sqrt(stat_sym²+sys²)'),
        'cross_checks': {
            'mask_approach_sys': _entry(mask_approach, 'km/s',
                'results/sigma_e_mask_systematic_N500.npz',
                'peak-to-peak/2 across 4 arc-mask approaches; overlaps F200mask, NOT added'),
            'D7_light_family': _entry(float(d7['sys_re_light']),
                'km/s', 'results/sigma_e_Re_grid_N500.npz',
                'peak-to-peak/2 over light-CoG R_e estimators only (6-pt, CaHK+G excluded)'),
            'Re_source_fullgrid': _entry(float(d7['sys_re']),
                'km/s', 'results/sigma_e_Re_grid_N500.npz',
                'conservative ceiling: full 7-pt grid incl. CaHK+G 2.90" (NOT folded)'),
            'CaHK_G_deviation': _entry(float(d7['cahk_dev']),
                'km/s', 'results/sigma_e_Re_grid_N500.npz',
                'σ_e(CaHK+G R_e=2.90") − headline; different I-map definition, cross-check'),
        },
    }

    # ── M⋆ — matched 2 R_e aperture, 5 estimators (2026-06-10 aperture overhaul) ──
    # HEADLINE = aperture-corrected total (empirical aperture + single-Sersic model wings;
    # the GAMA/Taylor+2011 fluxscale approach). All five estimators + the Sersic-only
    # systematic budget are emitted. Provenance: results/aperture_matched_photometry.npz.
    #
    # COSMOLOGY (2026-06-12): the project adopts Planck 2015 (H0=67.7, Om0=0.302) to match
    # the companion lens-modeling paper (Ferrami et al.); the underlying npz was fit under
    # the old H0=70, Om=0.3. A flux-calibrated SED stellar mass scales EXACTLY as D_L^2
    # (the M/L is color-driven, i.e. cosmology-independent), so switching cosmology is a
    # uniform additive shift of +2·log10(D_L_P15/D_L_70) = +0.0282 dex at our fixed
    # z_l=0.67564 — applied here analytically to every estimator (NO Bagpipes re-fit).
    # TODO: drop DLOGM_COSMO once aperture_matched_photometry.npz is regenerated under
    #       Planck 2015 (else it would be double-counted).
    DLOGM_COSMO = 0.0282  # +2*log10(4214.06/4079.23), D_L at z=0.67564, P15 vs H0=70
    am = np.load(RESULTS / 'aperture_matched_photometry.npz', allow_pickle=True)
    provM = 'results/aperture_matched_photometry.npz'
    # HEADLINE M* PIPELINE (2026-06-15 audit). Two coupled fixes vs the prior 11.50:
    #  (1) VALIDATED-FIT aperture correction: the model-dependent estimators use a Sérsic
    #      whose SHAPE is fixed to the photutils-validated per-band table fit (amplitude+sky
    #      to data), NOT the biased auto fit_sersic (r_eff 2.2-2.7", ellip→0) — removes a
    #      ~0.15 mag/band beyond-aperture over-correction (scripts/aperture_correction_validated.py).
    #  (2) SPECTRUM-CONSISTENT (quiescence-constrained) SED prior: the 4-band SED is degenerate
    #      (age-dust-M/L; ΔlnZ=-0.18, identical χ², scripts/sed_quiescence_check.py) and the flat-age
    #      prior slides onto a young+dusty branch (SFR~57) the KCWI absorption-line spectrum rules
    #      out. The headline adopts a passive SFH prior (old age, short tau, low dust); the
    #      fiducial↔quiescent spread is carried as the SFH-prior systematic.
    #  Net: 11.50 → 11.47 (validated photom -0.11 dex, partly offset by old-population M/L +0.10).
    #  All five estimators + budget come from results/mstar_headline_quiescent.npz.
    amv = np.load(RESULTS / 'aperture_correction_validated.npz', allow_pickle=True)
    qui = np.load(RESULTS / 'mstar_headline_quiescent.npz', allow_pickle=True)
    provV = 'results/aperture_correction_validated.npz'
    provQ = 'results/mstar_headline_quiescent.npz'

    def _est(kind, nre=2.0):
        p = am[f'logM_{kind}_{nre:g}']
        return float(p[1]) + DLOGM_COSMO, float(p[1] - p[0]), float(p[2] - p[1])

    def _estq(kind):  # validated photometry + quiescent (spectrum-consistent) SED prior
        p = qui[f'logM_{kind}_qui']
        return float(p[1]) + DLOGM_COSMO, float(p[1] - p[0]), float(p[2] - p[1])

    tot, tlo, thi = _estq('total')        # HEADLINE: validated photom + quiescent SED prior
    raw, rlo, rhi = _estq('raw')          # empirical lower bound (validated photom, quiescent prior)
    rac, _, _ = _estq('raw_apcorr')
    fil, _, _ = _estq('filled')
    ser_tot = _estq('sersic')[0]          # pure-Sérsic total under the quiescent prior
    st = np.load(RESULTS / 'sersic_total_systematic.npz', allow_pickle=True)
    msk = np.load(RESULTS / 'Mstar_masking_systematic.npz', allow_pickle=True)
    masking_sys = float(msk['masking_sys_10pct'])              # ±0.086 (under-arc ⊕ mask-def)
    apcorr_sys = float(qui['sys_apcorr_model_dex'])            # auto vs validated apcorr (quiescent prior)
    sfh_sys = float(qui['sys_sfh_prior_dex'])                  # fiducial vs quiescent SED prior (age-dust-M/L)
    mstar_sys_quad = float(np.sqrt(masking_sys**2 + apcorr_sys**2 + sfh_sys**2))   # SYMMETRIC sys
    # Total-light / cD-envelope ambiguity is ALREADY captured by the existing budget: the apcorr-model
    # term (auto more-extended Sérsic = the more-outer-light direction) + the quadrature sys reach a +1σ
    # upper bound of ~11.65, and the empirical curve-of-growth total (CoG@8") lands at 11.69 — only +0.04
    # dex beyond +1σ. So the CoG is recorded as a CROSS-CHECK only (NOT a separate systematic; a separate
    # +0.22 term would double-count the apcorr-model contribution). User decision 2026-06-17.
    cogf = np.load(RESULTS / 'cog_total_light.npz', allow_pickle=True)
    logM_cog = float(cogf['logM_cog8_H070'][1]) + DLOGM_COSMO
    reg['logMstar'] = {
        'central_10pct': _entry(tot, 'dex', provQ, 'validated photom + quiescent SED prior, 2 R_e, 10% (HEADLINE)'),
        'stat_lo_10pct': _entry(tlo, 'dex', provQ, 'median - p16 (Bagpipes posterior, quiescent prior)'),
        'err_hi_10pct':  _entry(thi, 'dex', provQ, 'p84 - median (Bagpipes posterior, quiescent prior)'),
        'sys_quad_dex':  _entry(mstar_sys_quad, 'dex', provQ,
                                'folded systematic = masking ⊕ apcorr-model ⊕ SFH-prior (quadrature)'),
        'empirical_raw': _entry(raw, 'dex', provQ, 'empirical aperture (lower bound), 2 R_e, quiescent prior'),
        'raw_apcorr':    _entry(rac, 'dex', provQ, 'raw + validated model wings, 2 R_e, quiescent prior'),
        'filled':        _entry(fil, 'dex', provQ, 'masked pixels filled with validated model, quiescent prior'),
        'sersic_total':  _entry(ser_tot, 'dex', provQ, 'pure single-Sersic total to inf, quiescent prior (Planck15)'),
        'sersic_total_sys_dex': _entry(float(st['total']), 'dex',
                                'results/sersic_total_systematic.npz',
                                'Sersic-only stat+sys quadrature (mask-dominated 0.10)'),
        'aperture_arcsec': _entry(2 * 2.097, 'arcsec', provV, '2 R_e matched elliptical aperture (R_e=2.097")'),
    }
    reg['logMstar']['masking_sys_dex'] = _entry(masking_sys, 'dex',
        'results/Mstar_masking_systematic.npz', 'peak-to-peak/2 across masking approaches')
    reg['logMstar']['apcorr_model_sys_dex'] = _entry(apcorr_sys, 'dex', provQ,
        'auto-fit vs validated-fit total — aperture-correction Sérsic-model systematic (2026-06-15)')
    reg['logMstar']['sfh_prior_sys_dex'] = _entry(sfh_sys, 'dex', provQ,
        'fiducial (flat-age) vs quiescent SED prior — age-dust-M/L degeneracy systematic (2026-06-15)')
    reg['logMstar']['cog_total_10pct'] = _entry(logM_cog, 'dex',
        'results/cog_total_light.npz',
        'empirical CoG@8" total (quiescent prior) = 11.69; CROSS-CHECK only — lands +0.04 dex beyond the '
        '+1σ upper bound, i.e. the total-light/cD-envelope ambiguity is already captured by the apcorr-model + sys budget')

    # ── M⋆ error budget — two parallel framings (2026-06-14) ──────────────────
    # σ_e's budget is a single quadrature of named components; M⋆ has TWO framings:
    #  (1) HEADLINE total-light estimator: the reported error IS the asymmetric
    #      Bagpipes posterior (stat — the age–dust–M/L outshining low-tail, present
    #      in all 5 estimators incl. raw), with the masking-approach systematic
    #      (±0.086) reported alongside. The model-choice spread is NOT folded into
    #      the headline error — it is captured by framing (2), kept as a cross-check
    #      so model choices are not double-counted into the empirical-aperture mass.
    #  (2) Sérsic-only estimator: a named-component quadrature (mask-dominated) — the
    #      M⋆ analogue of the σ_e budget table. Cross-check, NOT the headline error.
    # Both emitted; --render writes the `mstar_budget` block from these.
    stc = {str(name): float(val) for name, val in st['components']}
    reg['logMstar_budget'] = {
        'headline': {
            'estimator': 'validated photom + quiescent (spectrum-consistent) SED prior (2 R_e, 10% floor)',
            'stat_lo': _entry(tlo, 'dex', provQ,
                'median - p16 (Bagpipes posterior, quiescent prior — tight; degeneracy tail removed)'),
            'stat_hi': _entry(thi, 'dex', provQ, 'p84 - median (Bagpipes posterior, quiescent prior)'),
            'masking_sys': _entry(masking_sys, 'dex',
                'results/Mstar_masking_systematic.npz',
                'peak-to-peak/2 across masking approaches = under-arc ⊕ mask-def'),
            'masking_under_arc': _entry(float(msk['under_arc_10pct']), 'dex',
                'results/Mstar_masking_systematic.npz', 'raw↔filled (under-arc fill)'),
            'masking_def': _entry(float(msk['maskdef_10pct']), 'dex',
                'results/Mstar_masking_systematic.npz', 'per-band↔global mask definition'),
            'apcorr_model_sys': _entry(apcorr_sys, 'dex', provQ,
                'auto-fit vs validated-fit total (aperture-correction Sérsic-model choice; 2026-06-15)'),
            'sfh_prior_sys': _entry(sfh_sys, 'dex', provQ,
                'fiducial (flat-age) vs quiescent SED prior (age-dust-M/L degeneracy; 2026-06-15)'),
            'sys_quad': _entry(mstar_sys_quad, 'dex', provQ,
                'folded systematic = masking ⊕ apcorr-model ⊕ SFH-prior (the apcorr-model term already spans the more-outer-light / cD-envelope direction)'),
            'total_lo': _entry(float(np.hypot(tlo, mstar_sys_quad)), 'dex', provQ, 'stat_lo ⊕ sys_quad'),
            'total_hi': _entry(float(np.hypot(thi, mstar_sys_quad)), 'dex', provQ, 'stat_hi ⊕ sys_quad'),
        },
        'sersic_only': {
            'central': _entry(float(st['central']) + DLOGM_COSMO, 'dex',
                'results/sersic_total_systematic.npz', 'Sérsic-only central (Planck15)'),
            'stat': _entry(float(st['stat']), 'dex',
                'results/sersic_total_systematic.npz', 'Bagpipes posterior width'),
            'sys_quad': _entry(float(st['sys_quad']), 'dex',
                'results/sersic_total_systematic.npz', 'quadrature of the components below'),
            'total': _entry(float(st['total']), 'dex',
                'results/sersic_total_systematic.npz', 'stat ⊕ sys'),
            'components': {
                'mask':         _entry(stc['mask'], 'dex', 'results/sersic_total_systematic.npz',
                                       'arc-mask choice (expert↔global Sérsic-total) — dominant'),
                'model_form':   _entry(stc['model_form'], 'dex', 'results/sersic_total_systematic.npz',
                                       'single-Sérsic vs alternative model form'),
                'flux_floor':   _entry(stc['flux_floor'], 'dex', 'results/sersic_total_systematic.npz',
                                       '10%↔20% flux floor'),
                'fit_param':    _entry(stc['fit_param'], 'dex', 'results/sersic_total_systematic.npz',
                                       'Sérsic-n / fit-parameter spread'),
                'apcorr_recon': _entry(stc['apcorr_recon'], 'dex', 'results/sersic_total_systematic.npz',
                                       'apcorr ↔ pure-model reconstruction'),
            },
        },
    }

    # ── R_e ──────────────────────────────────────────────────────────────────
    # Best-mask (single-Sérsic + color/morph gate + WCS reg) F140W+F200LP CoG.
    # Adopted 2026-06-11 ("best mask throughout"); photutils-validated to ±0.002"
    # (scripts/validate_Re_photutils.py). Was the expert-mask 2.305".
    cog = np.load(RESULTS / 'Re_cog_reconciliation_bestmask.npz', allow_pickle=True)
    # Planck 2015 (H0=67.7, Om0=0.302), z=0.67564 — matched to companion paper (Ferrami et al.),
    # 2026-06-12. Was 7.04 under H0=70,Om=0.3 (a pure cosmological rescale, ×1.0331). astropy:
    # FlatLambdaCDM(67.7,0.302).kpc_proper_per_arcmin(0.67564) → 7.2764 kpc/arcsec; D_A=1500.9 Mpc.
    KPC = 7.2764  # kpc/arcsec at z=0.67564 (Planck 2015)
    # Headline R_e KEPT at 2.097" (best-mask CoG; user decision 2026-06-11). The Sérsic
    # ELLIPTICITY fix (multi-start fitter) shifts the best-mask CoG to
    # cog['headline_mean_raw']≈2.093" — a −0.004" change, 25× below the ±0.100" method
    # systematic and below reporting precision (both round to 2.10") — documented as
    # robustness, NOT adopted (keeps σ_e=267.31 without re-running the N=500 grid).
    re_hl = 2.097
    reg['R_e'] = {
        'central_arcsec': _entry(re_hl, 'arcsec',
                                 'results/Re_cog_reconciliation_bestmask.npz / final_sigma_e.py',
                                 'best-mask mean(F140W,F200LP) raw CoG r_max=6"'),
        'central_kpc': _entry(re_hl * KPC, 'kpc', 'R_e × 7.04 kpc/arcsec', 'R_e[arcsec]*KPC'),
        'method_sys_arcsec': _entry(float(cog['re_systematic_arcsec']), 'arcsec',
            'results/Re_cog_reconciliation_bestmask.npz',
            'peak-to-peak/2 across {raw CoG, sky-sub CoG, Sérsic} on best mask'),
    }

    # ── Hδ decision ──────────────────────────────────────────────────────────
    hd = np.load(RESULTS / 'sigma_e_hdelta_test_N500.npz', allow_pickle=True)
    reg['hdelta'] = {
        'decision': 'keep unmasked (a)',
        'max_delta_kms': _entry(float(hd['max_delta']), 'km/s', 'results/sigma_e_hdelta_test_N500.npz',
            'max|Δσ_e| over mask variants = information content, NOT contamination'),
        'is_outlier': bool(float(hd['hdelta_z_max']) > float(hd['global_z_p99'])),
    }

    # ── fixed metadata ───────────────────────────────────────────────────────
    reg['constants'] = {
        'z_deflector': _entry(0.67564, '', 'notebook 04', 'line-fit systemic'),
        'z_source': _entry(1.30263, '', '[O II] λλ3726,3729 doublet (red cube); AGEL DR2 gives 1.302', 'source redshift'),
        'kpc_per_arcsec': _entry(KPC, 'kpc/arcsec', 'Planck 2015 (H0=67.7,Om0=0.302)', 'at z=0.67564'),
        'ra_deg': _entry(31.55611, 'deg', 'target coordinates (ICRS)', 'R.A.'),
        'dec_deg': _entry(-1.23817, 'deg', 'target coordinates (ICRS)', 'Decl.'),
    }

    # ── photometry — matched 2 R_e total AB mags (the SED / M⋆ fluxes) ─────────
    FLUX_FLOOR = 0.10                       # adopted 10% systematic flux floor
    mag_floor_err = 2.5 * np.log10(1 + FLUX_FLOOR)
    reg['photometry'] = {
        'aperture': _entry(2 * 2.097, 'arcsec', provM, 'matched 2 R_e elliptical aperture (b/a=0.75)'),
        'mag_floor_err': _entry(mag_floor_err, 'mag', provM,
                                '2.5·log10(1+0.10) — adopted 10% systematic flux floor (propagated errors 0.03–1.6%)'),
        'bands': {str(fn): {'pivot_AA': float(pv), 'm_AB_total': float(mg)}
                  for fn, pv, mg in zip(amv['filter_names'], amv['pivot'], amv['mag_total_2'])},
    }

    # ── lens model — PROVISIONAL (Ferrami et al. draft 2026-06-12; NOT final) ──
    # Combined free-BH posterior; Planck-2015 (= our adopted cosmology). Do NOT freeze
    # into the manuscript until Ferrami's final MultiNest runs land (DRAFTING §3.3 flag G1).
    # Central = free-BH combined posterior (highest evidence; user decision 2026-06-12). UNCERTAINTY
    # = the 1σ errors across the FINAL SELECTED (free-BH) EPL models = the 1σ model envelope
    # (min median−σlo to max median+σhi over the 3 free-BH solutions) = 3.81–6.91e8, via
    # bh_mass_combine.py (user 2026-06-16). Excludes the fixed-BH runs (not in the final selection).
    mbh, mbh_lo, mbh_hi = 5.2e8, 1.39e8, 1.71e8
    reg['lens_model'] = {
        'provisional': 1,
        'M_BH': _entry(mbh, 'Msun',
                       'Ferrami et al. draft Table 1; free-BH (selected EPL) 1σ model envelope, scripts/bh_mass_combine.py',
                       'free-BH central; uncertainty = 1σ errors across selected free-BH EPL models (3.8–6.9e8); PROVISIONAL', err_lo=mbh_lo, err_hi=mbh_hi),
        'logM_BH': _entry(np.log10(mbh), 'dex', 'derived from M_BH', 'log10(M_BH); PROVISIONAL',
                          err_lo=np.log10(mbh) - np.log10(mbh - mbh_lo),
                          err_hi=np.log10(mbh + mbh_hi) - np.log10(mbh)),
        'theta_E': _entry(1.36, 'arcsec', 'Ferrami et al. draft §3.3.2',
                          'main-halo Einstein radius (free-BH best model); PROVISIONAL'),
        'gamma': _entry(1.31, '', 'Ferrami et al. draft §3.3.2',
                        'total mass-density slope (sub-isothermal); PROVISIONAL', err_sym=0.08),
    }
    # M•–σ offset vs Greene+2020 — the relation drawn in the M•–σ figure
    # (figures_paper4 cell 11: log M_BH = 8.03 + 4.24·log10(σ/160) + 0.43, pivot σ=160).
    # Δ = log M_BH(obs) − relation(σ_e). PROVISIONAL (tracks the provisional M_BH);
    # offset error carried from the logM_BH posterior (relation-param + σ_e terms
    # are sub-dominant and the whole quantity is provisional).
    sig_e = reg['sigma_e']['central']['value']
    g_pred = 8.03 + 4.24 * np.log10(sig_e / 160.) + 0.43
    lm_off = np.log10(mbh) - g_pred
    reg['lens_model']['mbh_sigma_offset'] = _entry(
        lm_off, 'dex',
        'Greene+2020 M•–σ (α,β,ε=8.03,4.24,0.43 @ σ=160; figures_paper4 cell 11) − PROVISIONAL logM_BH',
        'log10(M_BH) − [8.03 + 4.24·log10(σ_e/160) + 0.43]; negative = undermassive BH; PROVISIONAL',
        err_lo=np.log10(mbh) - np.log10(mbh - mbh_lo),
        err_hi=np.log10(mbh + mbh_hi) - np.log10(mbh))

    # ── SED-derived stellar population (Bagpipes posterior) ──────────────────
    # From the HEADLINE-consistent validated-fit total run (AGEL0206_aperVAL_total_2Re),
    # NOT the superseded original SED fit — keeps the age tied to the adopted photometry.
    # full stellar-population posterior from the HEADLINE fit (validated photom + quiescent,
    # spectrum-consistent SED prior) — the same Bagpipes run that sets the headline log M*.
    # The fiducial flat-age fit's young+dusty SFR (~57) is a degeneracy artifact ruled out by
    # the KCWI absorption-line spectrum; these passive properties supersede it.
    spd = {str(k).replace('prop_', ''): (float(qui[k][0]), float(qui[k][1]), float(qui[k][2]))
           for k in qui.files if k.startswith('prop_')}
    provSP = 'results/mstar_headline_quiescent.npz (Bagpipes posterior, quiescent prior)'
    _U = {'mass_weighted_age': ('Gyr', 'mass-weighted stellar age'),
          'formed_mass': ('dex', 'log10 total mass ever formed (pre mass-loss)'),
          'sfr': ('Msun/yr', 'SFR at observation (delayed-exp SFH)'),
          'ssfr': ('log10(1/yr)', 'specific SFR; quiescent if <-11'),
          'metallicity': ('Zsun', 'stellar metallicity'),
          'dust_Av': ('mag', 'dust attenuation A_V'),
          'tau': ('Gyr', 'exponential SFH e-folding time'),
          'age_form': ('Gyr', 'time since SF onset')}
    reg['sed'] = {}
    for k, (unit, desc) in _U.items():
        if k in spd:
            p, lo, hi = spd[k]
            reg['sed'][k] = _entry(p, unit, provSP, desc + ' (posterior median; p16/p84)',
                                   err_lo=lo, err_hi=hi)
    return reg


def _flatten(reg, prefix=''):
    """Yield (dotted_key, value_float) for every numeric 'value' leaf."""
    for k, v in reg.items():
        if k.startswith('_'):
            continue
        kk = f'{prefix}.{k}' if prefix else k
        if isinstance(v, dict):
            if 'value' in v and isinstance(v['value'], (int, float)):
                yield kk, float(v['value'])
            else:
                yield from _flatten(v, kk)


def print_table(reg):
    print('═' * 78)
    print('PAPER_VALUES registry  (scripts/paper_values.py)')
    print('═' * 78)
    se = reg['sigma_e']; b = reg['sigma_e_budget']
    print(f"σ_e(<R_e)   = {se['central']['value']:.2f}  "
          f"−{b['total_lo']['value']:.2f}/+{b['total_hi']['value']:.2f} km/s   "
          f"[stat −{se['stat_lo']['value']:.2f}/+{se['stat_hi']['value']:.2f} ⊕ sys ±{b['sys_quad']['value']:.2f}]")
    print(f"            symmetric total ±{b['total_sym']['value']:.2f} km/s")
    print('  budget components (km/s):')
    for name, c in b['components'].items():
        tag = ' (carried M11)' if c.get('carried') else ' (live)'
        print(f'    {name:11s} ±{c["value"]:.2f}{tag}')
    cc = b['cross_checks']
    print(f"  cross-checks: mask-approach ±{cc['mask_approach_sys']['value']:.2f}  "
          f"D7 light-family ±{cc['D7_light_family']['value']:.2f}")
    M = reg['logMstar']
    print(f"\nlog M⋆      = {M['central_10pct']['value']:.2f}  "
          f"+{M['err_hi_10pct']['value']:.2f}/−{M['stat_lo_10pct']['value']:.2f} (aperture-corrected total, 2 R_e, 10%)")
    print(f"  estimators: raw {M['empirical_raw']['value']:.2f} (empirical) · "
          f"raw+apcorr {M['raw_apcorr']['value']:.2f} · filled {M['filled']['value']:.2f} · "
          f"total {M['central_10pct']['value']:.2f} (headline) · Sérsic-total {M['sersic_total']['value']:.2f}")
    print(f"  Sérsic-only budget ±{M['sersic_total_sys_dex']['value']:.2f} (mask-dom)  "
          f"mask-sys ±{M['masking_sys_dex']['value']:.2f}")
    MB = reg['logMstar_budget']; H = MB['headline']; S = MB['sersic_only']
    print(f"  budget (1) headline error  +{H['stat_hi']['value']:.2f}/−{H['stat_lo']['value']:.2f} (posterior) "
          f"⊕ masking ±{H['masking_sys']['value']:.3f} (under-arc ±{H['masking_under_arc']['value']:.3f} + mask-def ±{H['masking_def']['value']:.3f})")
    print(f"  budget (2) Sérsic-only     stat ±{S['stat']['value']:.3f} ⊕ "
          + ' ⊕ '.join(f"{n.replace('_','-')} ±{c['value']:.3f}" for n, c in S['components'].items())
          + f" → ±{S['total']['value']:.2f}")
    R = reg['R_e']
    print(f"R_e         = {R['central_arcsec']['value']:.3f}\" = {R['central_kpc']['value']:.2f} kpc  "
          f"(method sys ±{R['method_sys_arcsec']['value']:.2f}\")")
    H = reg['hdelta']
    print(f"Hδ          = {H['decision']}  (max Δσ_e {H['max_delta_kms']['value']:.1f} km/s, "
          f"outlier={H['is_outlier']})")
    print('═' * 78)


def check(files, reg, tol=0.05):
    """Drift linter (ANCHORED). Only validates the canonical *headline statements*
    — `σ_e(<R_e) = N ± M`, its `asym −lo/+hi` form, and `log(M⋆/M☉) = N` — by
    parsing the asserted value out of those specific sentence patterns and
    comparing to the registry. This avoids the false-positive flood you get from
    fuzzy-proximity matching every nearby sweep / historical / cross-check number.

    A flagged line is one that *states a headline value* but states it WRONG
    (stale). Historical mentions ("was ±11.77") are written as ranges/prose and
    do not match the equals-anchored patterns, so they are not flagged.
    """
    se = reg['sigma_e']['central']['value']
    se_sym = reg['sigma_e_budget']['total_sym']['value']
    se_lo = reg['sigma_e_budget']['total_lo']['value']
    se_hi = reg['sigma_e_budget']['total_hi']['value']
    lm = reg['logMstar']['central_10pct']['value']
    # (compiled pattern, [(registry_value, group_index), ...], label)
    checks = [
        (re.compile(r'σ_e\(<R_e\)\s*=\s*\*{0,2}([\d.]+)\s*±\s*([\d.]+)'),
         [(se, 1), (se_sym, 2)], 'σ_e = central ± sym'),
        (re.compile(r'σ_e\(<R_e\)\s*=\s*([\d.]+)\s*[−\-]([\d.]+)\s*/\s*\+([\d.]+)'),
         [(se, 1), (se_lo, 2), (se_hi, 3)], 'σ_e = central −lo/+hi'),
        (re.compile(r'log\(M[⋆\*]/M[☉o]\)\s*=\s*\*{0,2}([\d.]+)'),
         [(lm, 1)], 'log M⋆ central'),
    ]
    # Lines that legitimately state a NON-live value (cross-check, history,
    # rounded, superseded) reuse the same sentence shape; skip them so only a
    # genuinely-stale *live headline* is flagged. (Markers, case-insensitive.)
    skip = re.compile(r'(was |old |narrow|cross-?check|superseded|histor|pre-M|→|->|'
                      r'rounded|Δ|delta|nomask|no arc|legacy|earlier|prior|M8|M9|M10|M11|'
                      r'reduction|fill|reach|supersed|§6cum|§7|original cube|demoted|'
                      r'\bvs\b|upper|aperture|Sérsic-total|\(10%\)|\(20%\)|pv-skip)', re.I)

    def _mismatch(stated_str, regval):
        """Precision-aware: round the registry value to the # of decimals the
        doc states, so a legitimately rounded headline (270 ↔ 269.62) is NOT a
        mismatch, while a genuinely stale value (265.76, 11.77) still is."""
        nd = len(stated_str.split('.')[1]) if '.' in stated_str else 0
        return abs(float(stated_str) - round(regval, nd)) > max(tol, 0.5 * 10 ** (-nd))

    issues = 0
    n_checked = 0
    for f in files:
        try:
            lines = Path(f).read_text().splitlines()
        except Exception:
            continue
        # File-level exemption: dated/superseded snapshots carry a
        # `<!-- pv-skip-file -->` marker near the top — their bodies record
        # historically-correct (now-old) numbers and must NOT be "fixed".
        if any('pv-skip-file' in ln for ln in lines[:25]):
            continue
        n_checked += 1
        for ln_no, line in enumerate(lines, 1):
            if skip.search(line):
                continue
            for pat, specs, label in checks:
                for m in pat.finditer(line):
                    for regval, gi in specs:
                        if _mismatch(m.group(gi), regval):
                            print(f'  ✗ {f}:{ln_no}  [{label}] states {m.group(gi)} but registry = '
                                  f'{regval:.2f} (Δ={float(m.group(gi))-regval:+.2f})')
                            issues += 1
    status = 'OK — all headline statements match the registry' if issues == 0 else f'{issues} MISMATCH(es)'
    skipped = len(files) - n_checked
    print(f'check: {status}  ({n_checked} live file(s) checked, {skipped} pv-skip-file exempt)')
    return issues


# ─────────────────────────────────────────────────────────────────────────────
# Deterministic block rendering (--render): generate marked Markdown blocks from
# the registry so the doc's headline numbers are GENERATED, never hand-typed.
# A managed region looks like:
#     <!-- PV:auto:headline -->
#     ...generated content (do not edit by hand)...
#     <!-- /PV:auto:headline -->
# `--render <files>` replaces the inside of each region with the fresh block.
# ─────────────────────────────────────────────────────────────────────────────
from decimal import Decimal, ROUND_HALF_UP                                      # noqa: E402


def _rh(x, n):
    """Round-half-up to n decimals → string (avoids banker's-rounding surprises,
    so 269.625→'269.63', 2.3045→'2.305', matching how the doc quotes numbers)."""
    q = Decimal(1).scaleb(-n) if n > 0 else Decimal(1)
    s = str(Decimal(str(float(x))).quantize(q, rounding=ROUND_HALF_UP))
    return s


def render_headline(reg):
    se, b = reg['sigma_e'], reg['sigma_e_budget']
    M, R, C = reg['logMstar'], reg['R_e'], reg['constants']
    c = se['central']['value']
    sym, lo, hi = b['total_sym']['value'], b['total_lo']['value'], b['total_hi']['value']
    re_src = b['components']['Re_source']['value']
    tot, stlo, sthi = M['central_10pct']['value'], M['stat_lo_10pct']['value'], M['err_hi_10pct']['value']
    raw, rac = M['empirical_raw']['value'], M['raw_apcorr']['value']
    fil, ser = M['filled']['value'], M['sersic_total']['value']
    ser_sys, msk = M['sersic_total_sys_dex']['value'], M['masking_sys_dex']['value']
    ap = M['aperture_arcsec']['value']
    re_as = _rh(R['central_arcsec']['value'], 3)
    re_kpc = _rh(float(re_as) * C['kpc_per_arcsec']['value'], 2)  # round-then-mult, matches CLAUDE.md
    re_sys = R['method_sys_arcsec']['value']
    z_d, z_s = C['z_deflector']['value'], C['z_source']['value']
    return (
f"""**Headline numbers** — *generated from `results/PAPER_VALUES.json` by `scripts/paper_values.py --render`; do not hand-edit inside the markers.*

- σ_e(<R_e) = **{_rh(c,2)} ± {_rh(sym,2)} km/s** (asym −{_rh(lo,2)} / +{_rh(hi,2)}); often rounded **{_rh(c,0)} ± {_rh(sym,0)} km/s**. Measured at the best-mask R_e=2.097″ aperture; the R_e-source systematic (±{_rh(re_src,2)}, best-mask CoG family) is folded into the budget (§3.1).
- log(M⋆/M☉) = **{_rh(tot,2)} +{_rh(sthi,2)}/−{_rh(stlo,2)} (stat, 10%)** — aperture-corrected total at the matched 2 R_e = {_rh(ap,2)}″ elliptical aperture (empirical aperture + single-Sérsic model wings; GAMA/Taylor+2011 fluxscale). Five estimators: raw **{_rh(raw,2)}** (empirical) · raw+apcorr **{_rh(rac,2)}** · filled **{_rh(fil,2)}** · **total {_rh(tot,2)} (headline)** · Sérsic-total **{_rh(ser,2)}**. Sérsic-only systematic ±{_rh(ser_sys,2)} dex; masking-approach ±{_rh(msk,2)} dex. See §2.1.1b / §3.2.
- R_e = **{re_as}″ = {re_kpc} kpc** (best-mask F140W+F200LP CoG; method systematic ±{_rh(re_sys,2)}″). **Supersedes 2.305″** (expert mask).
- z_deflector = **{_rh(z_d,5)}**; z_source = **{_rh(z_s,3)}**.""")


def render_budget(reg):
    b = reg['sigma_e_budget']; se = reg['sigma_e']
    comp = b['components']
    rows = [
        ('stat (N=500)', se['stat_sym']['value'], f"asym −{_rh(se['stat_lo']['value'],2)}/+{_rh(se['stat_hi']['value'],2)}; pooled 3 SPS × 15 deg (marginalizes SPS + degree)"),
        ('I-shape', comp['Ishape']['value'], comp['Ishape']['formula']),
        ('arc masking', comp['F200mask']['value'], comp['F200mask']['formula']),
        ('centering', comp['centering']['value'], comp['centering']['formula']),
        ('fit-window', comp['fitwindow']['value'], comp['fitwindow']['formula']),
        ('reduction-pass', comp['reduction']['value'], comp['reduction']['formula']),
        ('R_e-source (best-mask grid)', comp['Re_source']['value'], comp['Re_source']['formula']),
    ]
    out = ["**σ_e systematic budget** — *generated by `scripts/paper_values.py --render`; do not hand-edit inside the markers.*",
           "", "| Component | ± km/s | Note |", "|---|---|---|"]
    for name, val, note in rows:
        out.append(f"| {name} | {_rh(val,2)} | {note} |")
    out.append(f"| **TOTAL (sym)** | **{_rh(b['total_sym']['value'],2)}** | quadrature (sys ±{_rh(b['sys_quad']['value'],2)} ⊕ stat) |")
    out.append(f"| **TOTAL (asym)** | **−{_rh(b['total_lo']['value'],2)} / +{_rh(b['total_hi']['value'],2)}** | preserves stat-side skew |")
    cc = b['cross_checks']
    out.append("")
    out.append(f"Cross-checks (not added): arc-mask-definition ±{_rh(cc['mask_approach_sys']['value'],2)} (overlaps F200 mask) · full-grid R_e-source ±{_rh(cc['Re_source_fullgrid']['value'],2)} (incl. CaHK+G, conservative ceiling) · CaHK+G(2.90″) deviation {cc['CaHK_G_deviation']['value']:+.2f} km/s.")
    return "\n".join(out)


def render_mstar_budget(reg):
    """M⋆ final error budget — a single named-component quadrature (now parallel to σ_e),
    plus the Sérsic-only quadrature kept as a model-systematic cross-check."""
    M, B = reg['logMstar'], reg['logMstar_budget']
    H, S = B['headline'], B['sersic_only']
    tot = M['central_10pct']['value']
    sysq = M['sys_quad_dex']['value']
    out = [
        "**M⋆ error budget** — *generated by `scripts/paper_values.py --render`; do not hand-edit inside the markers.* "
        "The headline error is a single named-component quadrature (like σ_e): the Bagpipes statistical "
        "posterior combined with the masking, aperture-correction-model, and SFH-prior systematics. The "
        "Sérsic-only quadrature is retained below as an independent model-path cross-check.",
        "",
        f"**(1) Headline — log(M⋆/M☉) = {_rh(tot,2)} = {_rh(tot,1)} (validated photometry + quiescent, "
        "spectrum-consistent SED prior; matched 2 R_e, 10% floor):**",
        "",
        "| Component | ± dex | Note |",
        "|---|---|---|",
        f"| stat (Bagpipes posterior) | +{_rh(H['stat_hi']['value'],3)} / −{_rh(H['stat_lo']['value'],3)} | "
        "quiescent prior — tight; the age–dust–M/L young-solution tail is removed by the spectrum-motivated prior |",
        f"| masking-approach | {_rh(H['masking_sys']['value'],3)} | "
        f"under-arc (raw↔filled) {_rh(H['masking_under_arc']['value'],3)} ⊕ mask-def (per-band↔global) {_rh(H['masking_def']['value'],3)} |",
        f"| aperture-correction model | {_rh(H['apcorr_model_sys']['value'],3)} | "
        "auto-fit vs validated-fit Sérsic shape (the over-correction this audit removed) |",
        f"| SFH prior | {_rh(H['sfh_prior_sys']['value'],3)} | "
        "fiducial flat-age vs quiescent SED prior (the age–dust–M/L degeneracy the 4-band SED cannot break) |",
        f"| **sys (quadrature)** | **{_rh(sysq,3)}** | masking ⊕ apcorr-model ⊕ SFH-prior |",
        f"| **REPORTED** | **+{_rh(H['stat_hi']['value'],2)} / −{_rh(H['stat_lo']['value'],2)} (stat) "
        f"± {_rh(sysq,2)} (sys)** | one-decimal: {_rh(tot,1)} ± {_rh(H['stat_hi']['value'],1)} (stat) ± {_rh(sysq,1)} (sys) |",
        f"| _CoG cross-check_ | _{_rh(M['cog_total_10pct']['value'],2)}_ | empirical curve-of-growth total-light "
        "(cD envelope); lands +0.04 dex beyond +1σ → total-light ambiguity already captured by the apcorr-model + sys budget, NOT added separately |",
        "",
        f"**(2) Sérsic-only estimator (log(M⋆/M☉) = {_rh(M['sersic_total']['value'],2)}) — independent "
        "model-path systematic envelope (cross-check):**",
        "",
        "| Component | ± dex | Note |",
        "|---|---|---|",
        f"| stat | {_rh(S['stat']['value'],3)} | Bagpipes posterior width |",
    ]
    for name, c in S['components'].items():
        out.append(f"| {name.replace('_','-')} | {_rh(c['value'],3)} | {c['formula']} |")
    out.append(f"| **TOTAL** | **{_rh(S['total']['value'],2)}** | "
               f"stat ±{_rh(S['stat']['value'],2)} ⊕ sys ±{_rh(S['sys_quad']['value'],2)} (mask-dominated) |")
    out.append("")
    out.append("Five M⋆ estimators (empirical→model, validated photometry + quiescent SED prior, matched 2 R_e, Planck 2015): "
               f"raw **{_rh(M['empirical_raw']['value'],2)}** (empirical lower bound) · "
               f"raw+apcorr **{_rh(M['raw_apcorr']['value'],2)}** · filled **{_rh(M['filled']['value'],2)}** · "
               f"**total {_rh(tot,2)} (headline)** · Sérsic-total **{_rh(M['sersic_total']['value'],2)}**.")
    return "\n".join(out)


BLOCKS = {'headline': render_headline, 'budget': render_budget, 'mstar_budget': render_mstar_budget}


def do_render(files, reg):
    pat = re.compile(r'(<!-- PV:auto:(\w+) -->)(.*?)(<!-- /PV:auto:\2 -->)', re.S)
    total = 0
    for f in files:
        p = Path(f)
        try:
            txt = p.read_text()
        except Exception:
            print(f'  ✗ {f}: cannot read'); continue
        def _sub(m):
            name = m.group(2)
            if name not in BLOCKS:
                print(f'  ⚠ {f}: unknown PV:auto block "{name}" — left as-is'); return m.group(0)
            return f'{m.group(1)}\n{BLOCKS[name](reg)}\n{m.group(4)}'
        new, n = pat.subn(_sub, txt)
        if n:
            p.write_text(new); total += n
            print(f'  ✓ {f}: rendered {n} block(s)')
        else:
            print(f'  – {f}: no PV:auto regions')
    print(f'render: {total} block(s) generated from the registry')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', nargs='*', default=None, help='files to lint for drift')
    ap.add_argument('--render', nargs='*', default=None,
                    help='files with <!-- PV:auto:NAME --> regions to regenerate from the registry')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    reg = build()
    out = RESULTS / 'PAPER_VALUES.json'
    out.write_text(json.dumps(reg, indent=2))
    rendering = args.render is not None
    checking = args.check is not None
    if not args.quiet and not checking and not rendering:
        print_table(reg)
        print(f'\n→ wrote {out}')
    if rendering:
        print('── render ──')
        do_render(args.render, reg)
    if checking:
        print('── drift check ──')
        check(args.check, reg)


if __name__ == '__main__':
    main()
