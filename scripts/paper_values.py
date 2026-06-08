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
    """Pool the 3-SPS wild-bootstrap σ samples from the headline new_clean_hei
    caches → p16/p50/p84 (the headline σ_e central value + stat error)."""
    caches = sorted(glob.glob(str(RESULTS / 'run_wide_sigma_e' / 'new_clean_hei'
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
    prov = 'results/run_wide_sigma_e/new_clean_hei/wR3800_5400_arcmask_{fsps,emiles,xsl}_*_N500.npz'
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
    d7 = np.load(RESULTS / 'sigma_e_Re_systematic_wide_N500.npz', allow_pickle=True)
    re_source = float(d7['sys_re'])  # full peak-to-peak/2 across 4 R_e estimators
    msys = np.load(RESULTS / 'sigma_e_mask_systematic_N500.npz', allow_pickle=True)
    mask_approach = float(msys['sys_mask'])  # cross-check, NOT a separate budget line

    # named components (carried M11 constants pointed at their caches + live R_e-source)
    components = {
        'Ishape':    _entry(2.27, 'km/s', 'results/ishape_sweep_wR3800_5400_arcmask/ (M11)',
                            '10-shape peak-to-peak/2 on NEW cube + M10', carried=1),
        'F200mask':  _entry(6.65, 'km/s', 'results/maskweight_sweep_wR3800_5400_arcmask/ (M11)',
                            '(w00-w100)/2; larger-of-two vs mask-approach 5.85 (no double-count)', carried=1),
        'frame':     _entry(5.0, 'km/s', 'structural (per-SPS native vac/air frame)',
                            'carried constant', carried=1),
        'centering': _entry(4.0, 'km/s', 'HST WCS centroid', 'carried constant', carried=1),
        'fitwindow': _entry(3.82, 'km/s', 'results/run_wide_sigma_e/ 3-window (M11)',
                            'peak-to-peak/2 across 3 fit windows', carried=1),
        'reduction': _entry(3.45, 'km/s', 'results/run_wide_sigma_e/{new,headline}_clean_hei (M10)',
                            'half-Δ between NEW and OLD reductions', carried=1),
        'Re_source': _entry(re_source, 'km/s', 'results/sigma_e_Re_systematic_wide_N500.npz',
                            'peak-to-peak/2 across 4 R_e estimators (full spread, user 2026-06-08)'),
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
            'D7_light_family': _entry((lambda a: (a.max()-a.min())/2)(
                np.array([float(d7['p50'][list(d7['labels']).index(l)])
                          for l in ('F140W', 'mean', 'F200LP')])),
                'km/s', 'results/sigma_e_Re_systematic_wide_N500.npz',
                'peak-to-peak/2 over light-R_e estimators only (CaHK+G excluded)'),
        },
    }

    # ── M⋆ ───────────────────────────────────────────────────────────────────
    m = np.load(RESULTS / 'photometry_systematics_Mstar.npz', allow_pickle=True)
    s10 = m['samples_perband_raw_10pct']; f10 = float(np.median(m['samples_perband_filled_10pct']))
    med10 = float(np.median(s10)); p16_10 = float(np.percentile(s10, 16)); p84_10 = float(np.percentile(s10, 84))
    err_hi10 = _quad(p84_10 - med10, f10 - med10)  # stat ⊕ one-sided fill-in
    s20 = m['samples_perband_raw_20pct']
    med20 = float(np.median(s20)); st20 = float((np.percentile(s20, 84) - np.percentile(s20, 16)) / 2)
    f20 = float(np.median(m['samples_perband_filled_20pct']))
    provM = 'results/photometry_systematics_Mstar.npz'
    reg['logMstar'] = {
        'central_10pct': _entry(med10, 'dex', provM, 'median(samples_perband_raw_10pct)'),
        'stat_lo_10pct': _entry(med10 - p16_10, 'dex', provM, 'median - p16'),
        'err_hi_10pct':  _entry(err_hi10, 'dex', provM, 'quad(p84-med, fillmed-med) [one-sided +sys]'),
        'fill_reach_10pct': _entry(f10, 'dex', provM, 'median(samples_perband_filled_10pct)'),
        'central_20pct': _entry(med20, 'dex', provM, 'median(samples_perband_raw_20pct)'),
        'stat_sym_20pct': _entry(st20, 'dex', provM, '(p84-p16)/2'),
        'fill_reach_20pct': _entry(f20, 'dex', provM, 'median(samples_perband_filled_20pct)'),
    }
    msk = np.load(RESULTS / 'Mstar_masking_systematic.npz', allow_pickle=True)
    reg['logMstar']['masking_sys_dex'] = _entry(float(msk['masking_sys_10pct']), 'dex',
        'results/Mstar_masking_systematic.npz', 'peak-to-peak/2 across masking approaches')

    # ── R_e ──────────────────────────────────────────────────────────────────
    cog = np.load(RESULTS / 'Re_cog_reconciliation.npz', allow_pickle=True)
    KPC = 7.04  # kpc/arcsec at z=0.67564 (CLAUDE.md)
    re_hl = float(cog['headline_mean_raw'])
    reg['R_e'] = {
        'central_arcsec': _entry(re_hl, 'arcsec', 'results/Re_cog_reconciliation.npz / final_sigma_e.py',
                                 'mean(F140W,F200LP) raw CoG r_max=6"'),
        'central_kpc': _entry(re_hl * KPC, 'kpc', 'R_e × 7.04 kpc/arcsec', 'R_e[arcsec]*KPC'),
        'method_sys_arcsec': _entry(float(cog['re_systematic_arcsec']), 'arcsec',
            'results/Re_cog_reconciliation.npz', 'peak-to-peak/2 across {raw CoG, sky-sub CoG, Sérsic}'),
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
          f"+{M['err_hi_10pct']['value']:.2f}/−{M['stat_lo_10pct']['value']:.2f} (10%)   "
          f"[{M['central_20pct']['value']:.2f}±{M['stat_sym_20pct']['value']:.2f} (20%)]  "
          f"fill→{M['fill_reach_10pct']['value']:.2f}  mask-sys ±{M['masking_sys_dex']['value']:.2f}")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', nargs='*', default=None, help='files to lint for drift')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    reg = build()
    out = RESULTS / 'PAPER_VALUES.json'
    out.write_text(json.dumps(reg, indent=2))
    # In --check mode, default to quiet (suppress the table) so the file list
    # after --check is never ambiguous with the table output.
    checking = args.check is not None
    if not args.quiet and not checking:
        print_table(reg)
        print(f'\n→ wrote {out}')
    if checking:
        print('── drift check ──')
        check(args.check, reg)


if __name__ == '__main__':
    main()
