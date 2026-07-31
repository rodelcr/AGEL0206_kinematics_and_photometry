#!/usr/bin/env python
"""Combine the three free-BH lens-model posteriors (Ferrami et al., DESJ0206 SMBH paper)
into a single 'final posterior' for the central black-hole mass.

User decision 2026-06-12: the canonical M_BH is the FREE-BH solution (BH position free
within 0.5" of the light centre), which is the highest-evidence family in Ferrami Table 1
(cEPL+free Dln E=0, bEPL+free -0.08, sEPL+free -1.43). The three medians agree to ~10%,
so we report the model-averaged posterior across all three.

CAVEAT: this combines the Table-1 SUMMARY statistics (median + asymmetric 1-sigma) as a
split-normal mixture, NOT Giovanni's actual MultiNest chains. The rigorous version stacks
the real posterior samples; because the three free-BH posteriors are near-identical, the
summary-level result is robust. Reproduces the numbers quoted in DRAFTING_FACTS §3.3.1.

Output (units 1e8 Msun):
  equal-weight    median ~5.22  -1.31/+1.47   (canonical headline: 5.2 -1.3/+1.5 e8)
  evidence-wtd    median ~5.29  -1.35/+1.45   (identical within rounding)
  conservative 2sigma envelope across the 3 models: ~2.6-8.3 e8
"""
import numpy as np

# free-BH posteriors from Ferrami Table 1, units 1e8 Msun: (median, sigma_lo, sigma_hi, dlnE)
MODELS = {
    'cEPL+free': (5.51, 1.40, 1.40,  0.00),
    'sEPL+free': (5.01, 1.20, 1.50, -1.43),
    'bEPL+free': (5.11, 1.30, 1.50, -0.08),
}
# fixed-BH (point-mass at the light centre) solutions from the SAME Table 1, units 1e8 Msun.
# Lower Bayesian evidence than free-BH (Δlnℰ −3.0 to −3.65), but Ferrami's draft PROSE/best-model
# figure favour them — the fixed-vs-free choice is UNRESOLVED in the draft. They sit 3–4× higher and
# must be included in the provisional uncertainty (user decision 2026-06-16: "full uncertainty across
# the table"). cEPL 1.61±0.21e9, sEPL 2.10 +0.19/−0.30 e9, bEPL 1.86±0.20e9.
MODELS_FIXED = {
    'cEPL+fixed': (16.1, 2.1, 2.1),
    'sEPL+fixed': (21.0, 3.0, 1.9),
    'bEPL+fixed': (18.6, 2.0, 2.0),
}
MBH_CENTRAL = 5.2   # adopted central = free-BH combined (highest-evidence; user decision 2026-06-12)
N = 2_000_000
SEED = 20260612


def _split_normal(med, slo, shi, n, rng):
    """Two-piece normal matching median + asymmetric 1-sigma errors."""
    z = rng.standard_normal(n)
    return np.where(z < 0, med + z * slo, med + z * shi)


def mixture(weights, rng):
    w = np.array([weights[k] for k in MODELS], float)
    w /= w.sum()
    counts = rng.multinomial(N, w)
    return np.concatenate([_split_normal(*MODELS[k][:3], c, rng)
                           for k, c in zip(MODELS, counts)])


def _report(name, s):
    p16, p50, p84 = np.percentile(s, [16, 50, 84])
    p2_5, p97_5 = np.percentile(s, [2.5, 97.5])
    print(f"{name:14s} M_BH = {p50:.2f} -{p50-p16:.2f}/+{p84-p50:.2f} e8 Msun"
          f"   [2sigma {p2_5:.2f}-{p97_5:.2f}]   sym ~+-{(p84-p16)/2:.2f}")


def main():
    rng = np.random.default_rng(SEED)
    _report('equal-weight', mixture({k: 1 for k in MODELS}, rng))
    _report('evidence-wtd', mixture({k: np.exp(MODELS[k][3]) for k in MODELS}, rng))
    los, his = [], []
    for k, (m, slo, shi, _) in MODELS.items():
        s = _split_normal(m, slo, shi, 500_000, rng)
        lo, hi = np.percentile(s, [2.5, 97.5])
        los.append(lo); his.append(hi)
    print(f"\nconservative 2sigma envelope (min/max of 2sigma pctiles, 3 free-BH models): "
          f"{min(los):.2f}-{max(his):.2f} e8 Msun")

    # ADOPTED provisional uncertainty (user 2026-06-16, final): the 1σ errors across the FINAL
    # SELECTED (free-BH) EPL models — the 1σ model envelope = min(median−σlo) to max(median+σhi) over
    # the three free-BH solutions (3.81–6.91 e8 ≈ 5.2 −1.4/+1.7). Excludes the FIXED-BH runs (not in
    # the final selection). (2σ envelope and free+fixed envelope kept below as references only.)
    p_lo = min(m - slo for m, slo, shi, _ in MODELS.values())   # 1σ model envelope low
    p_hi = max(m + shi for m, slo, shi, _ in MODELS.values())   # 1σ model envelope high
    elo, ehi = MBH_CENTRAL - p_lo, p_hi - MBH_CENTRAL
    lc = np.log10(MBH_CENTRAL * 1e8)
    print(f"\n*** ADOPTED: free-BH (selected EPL) 1σ model envelope ***")
    print(f"    span {p_lo:.2f}-{p_hi:.2f} e8  (={p_lo*1e8:.2e} to {p_hi*1e8:.2e} Msun)")
    print(f"    M_BH = {MBH_CENTRAL:.1f} +{ehi:.2f}/-{elo:.2f} e8   (central = free-BH combined)")
    print(f"    log M_BH = {lc:.2f} +{np.log10(p_hi*1e8)-lc:.2f}/-{lc-np.log10(p_lo*1e8):.2f}")
    print(f"    -> paper_values.py: mbh=5.2e8, mbh_lo={elo*1e8:.3e}, mbh_hi={ehi*1e8:.3e}")
    # over-conservative reference (free+fixed) — NOT adopted:
    allmod = {**{k: v[:3] for k, v in MODELS.items()}, **MODELS_FIXED}
    fp16 = min(m - slo for m, slo, shi in allmod.values())
    fp84 = max(m + shi for m, slo, shi in allmod.values())
    print(f"    [ref only, free+fixed envelope (NOT adopted): {fp16:.2f}-{fp84:.2f} e8]")


if __name__ == '__main__':
    main()
