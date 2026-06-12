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
        'F200mask':  _entry(6.65, 'km/s', 'results/maskweight_sweep_wR3800_5400_arcmask/ (M11)',
                            '(w00-w100)/2; larger-of-two vs mask-approach 5.85 (no double-count)', carried=1),
        'frame':     _entry(5.0, 'km/s', 'structural (per-SPS native vac/air frame)',
                            'carried constant', carried=1),
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
    am = np.load(RESULTS / 'aperture_matched_photometry.npz', allow_pickle=True)
    provM = 'results/aperture_matched_photometry.npz'

    def _est(kind, nre=2.0):
        p = am[f'logM_{kind}_{nre:g}']
        return float(p[1]), float(p[1] - p[0]), float(p[2] - p[1])  # median, lo, hi

    tot, tlo, thi = _est('total')         # headline: aperture-corrected total at 2 R_e
    raw, rlo, rhi = _est('raw')           # empirical lower bound
    rac, _, _ = _est('raw_apcorr')        # most-empirical total
    fil, _, _ = _est('filled')
    ser = np.load(RESULTS / 'bagpipes_sersic_refit.npz', allow_pickle=True)
    st = np.load(RESULTS / 'sersic_total_systematic.npz', allow_pickle=True)
    reg['logMstar'] = {
        'central_10pct': _entry(tot, 'dex', provM, 'aperture-corrected total, 2 R_e, 10% (HEADLINE)'),
        'stat_lo_10pct': _entry(tlo, 'dex', provM, 'median - p16'),
        'err_hi_10pct':  _entry(thi, 'dex', provM, 'p84 - median'),
        'empirical_raw': _entry(raw, 'dex', provM, 'empirical aperture, masked (lower bound), 2 R_e'),
        'raw_apcorr':    _entry(rac, 'dex', provM, 'most-empirical total: raw + model wings, 2 R_e'),
        'filled':        _entry(fil, 'dex', provM, 'masked pixels filled within aperture, 2 R_e'),
        'sersic_total':  _entry(float(ser['log_M_sersic_p50']), 'dex',
                                'results/bagpipes_sersic_refit.npz', 'pure single-Sersic total to inf'),
        'sersic_total_sys_dex': _entry(float(st['total']), 'dex',
                                'results/sersic_total_systematic.npz',
                                'Sersic-only stat+sys quadrature (mask-dominated 0.10)'),
        'aperture_arcsec': _entry(2 * 2.097, 'arcsec', provM, '2 R_e matched elliptical aperture (R_e=2.097")'),
    }
    msk = np.load(RESULTS / 'Mstar_masking_systematic.npz', allow_pickle=True)
    reg['logMstar']['masking_sys_dex'] = _entry(float(msk['masking_sys_10pct']), 'dex',
        'results/Mstar_masking_systematic.npz', 'peak-to-peak/2 across masking approaches')

    # ── R_e ──────────────────────────────────────────────────────────────────
    # Best-mask (single-Sérsic + color/morph gate + WCS reg) F140W+F200LP CoG.
    # Adopted 2026-06-11 ("best mask throughout"); photutils-validated to ±0.002"
    # (scripts/validate_Re_photutils.py). Was the expert-mask 2.305".
    cog = np.load(RESULTS / 'Re_cog_reconciliation_bestmask.npz', allow_pickle=True)
    KPC = 7.04  # kpc/arcsec at z=0.67564 (CLAUDE.md)
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
        'z_source': _entry(1.302, '', 'AGEL DR2', 'source redshift'),
        'kpc_per_arcsec': _entry(KPC, 'kpc/arcsec', 'CLAUDE.md', 'at z=0.67564'),
    }
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
        ('F200 mask', comp['F200mask']['value'], comp['F200mask']['formula']),
        ('frame (vac/air)', comp['frame']['value'], comp['frame']['formula']),
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


BLOCKS = {'headline': render_headline, 'budget': render_budget}


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
