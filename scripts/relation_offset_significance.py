"""Significance of AGEL0206's offset from the Greene+2020 M•–σ and M•–M⋆ relations (2026-06-17).

"How many sigma is our point from each relation?" — combined in quadrature ("pythagorean"):

    N_sigma = |Δ| / sqrt( σ_M•(meas)²  +  (β·σ_x)²  +  ε_int²  [+ σ_α² + (σ_β·log(x/pivot))²] )

where, for each relation log M• = α + β·log(x/pivot):
  Δ           = log M•(obs) − [α + β·log(x_obs/pivot)]      (offset from the MEAN relation; NEGATIVE = undermassive)
  σ_M•(meas)  = our log M_BH uncertainty (registry; side toward the relation)
  β·σ_x       = our σ_e / log M⋆ uncertainty propagated through the slope (side that closes the gap)
  ε_int       = the relation's INTRINSIC SCATTER (dex) — the dominant term; this is the population
                dispersion about the mean, which is what "distance from the relation" must be measured in
  σ_α, σ_β    = relation-parameter (zero-point, slope) uncertainties — added in the "full" variant only

Greene+2020 parameters are taken verbatim from the paper's figure code
(AGEL_0206_ApJL_Figures/figures_paper4_2026-06-08.ipynb, M•–σ/M•–M⋆ cells):
  M•–σ  : α=8.03±0.06, β=4.24±0.25, ε_int=0.43±0.04, pivot σ =160 km/s
  M•–M⋆ : α=7.89±0.09, β=1.33±0.12, ε_int=0.65±0.05, pivot M⋆=3×10¹⁰ M☉

NOTE (2026-06-17): the registry's `mbh_sigma_offset` (−0.69) was computed against
[α + β·log(σ/160) + ε_int] — i.e. the intrinsic scatter was ADDED to the mean relation. That is not
the offset from the relation (and not where the relation is drawn in the figure, which is the mean
locus). The correct mean-relation offset is computed here (≈ −0.26 dex).

Run from repo root:  python -m scripts.relation_offset_significance
"""
import json
import numpy as np

REG = json.load(open("results/PAPER_VALUES.json"))
LN10 = np.log(10.0)

# Greene+2020 relation parameters — VERIFIED against Greene, Strader & Ho 2020 (ARA&A 58, 257;
# arXiv:1911.09678) Supplemental Table 5, "Early" (early-type) subsample — the correct sample for the
# massive passive E/S0 deflector. Confirmed value-for-value against the in-repo machine table
# ../AGEL_Mbh_sigma/Greene20_Supple_5.csv AND the supplemental PDF ../AGEL_Mbh_sigma/aa58_greene_supmat.pdf
# (central values + uncertainties). Pivots: σ/160 km/s (main text), M⋆/3e10 M⊙ (Suppl. Tab. 5 caption).
# Full provenance + field-standard methodology: NOTES_relation_offset_provenance_2026-06-17.md
REL = {
    "M•–σ": dict(alpha=8.03, beta=4.24, eps=0.43, pivot=160.0,
                 a_err=0.06, b_err=0.25),
    "M•–M⋆": dict(alpha=7.89, beta=1.33, eps=0.65, pivot=3e10,
                  a_err=0.09, b_err=0.12),
}


def main():
    L = REG["lens_model"]
    logMbh = L["logM_BH"]["value"]
    sMbh_hi = L["logM_BH"]["err_hi"]      # upward (toward the relation, since point is below)
    sMbh_lo = L["logM_BH"]["err_lo"]

    se = REG["sigma_e"]; seb = REG["sigma_e_budget"]
    sigma = se["central"]["value"]
    sigma_err = seb["total_sym"]["value"]                       # km/s, stat⊕sys
    sig_logsigma = sigma_err / (sigma * LN10)                   # dex

    M = REG["logMstar"]; Mb = REG["logMstar_budget"]["headline"]
    logMstar = M["central_10pct"]["value"]
    sig_logMstar = max(Mb["total_hi"]["value"], Mb["total_lo"]["value"])   # dex (stat⊕sys)

    xobs = {"M•–σ": sigma, "M•–M⋆": 10 ** logMstar}
    sig_logx = {"M•–σ": sig_logsigma, "M•–M⋆": sig_logMstar}

    out = {}
    print(f"log M_BH (obs) = {logMbh:.3f}  (+{sMbh_hi:.3f}/−{sMbh_lo:.3f})   [PROVISIONAL]")
    print(f"σ_e = {sigma:.1f} ± {sigma_err:.1f} km/s  (σ_logσ={sig_logsigma:.4f} dex)")
    print(f"log M⋆ = {logMstar:.3f} ± {sig_logMstar:.3f} dex\n")
    print(f"{'relation':8s} {'pred':>7s} {'Δ(dex)':>8s} {'σ_M•':>6s} {'β·σ_x':>7s} {'ε_int':>6s} "
          f"{'σ_tot':>7s} {'Nσ':>6s} {'Nσ(+α,β)':>9s} {'Nσ(meas)':>9s}")
    print("-" * 90)
    for name, p in REL.items():
        x = xobs[name]; slx = sig_logx[name]
        pred = p["alpha"] + p["beta"] * np.log10(x / p["pivot"])      # MEAN relation
        delta = logMbh - pred
        sMbh = sMbh_hi if delta < 0 else sMbh_lo                       # error toward the relation
        term_x = p["beta"] * slx
        term_eps = p["eps"]
        sig_tot = np.sqrt(sMbh ** 2 + term_x ** 2 + term_eps ** 2)
        nsig = abs(delta) / sig_tot
        # variant including relation-parameter (zero-point, slope) uncertainties
        term_a = p["a_err"]; term_b = p["b_err"] * abs(np.log10(x / p["pivot"]))
        sig_full = np.sqrt(sig_tot ** 2 + term_a ** 2 + term_b ** 2)
        nsig_full = abs(delta) / sig_full
        # scatter-FREE distance: vs the singular mean line, MEASUREMENT errors only (no ε_int)
        sig_meas = np.sqrt(sMbh ** 2 + term_x ** 2)
        nsig_meas = abs(delta) / sig_meas
        print(f"{name:8s} {pred:7.3f} {delta:+8.3f} {sMbh:6.3f} {term_x:7.3f} {term_eps:6.3f} "
              f"{sig_tot:7.3f} {nsig:6.2f} {nsig_full:9.2f} {nsig_meas:9.2f}")
        key = "sigma" if name == "M•–σ" else "mstar"
        out[f"{key}_pred"] = pred; out[f"{key}_offset"] = delta
        out[f"{key}_sigtot"] = sig_tot; out[f"{key}_nsigma"] = nsig
        out[f"{key}_nsigma_full"] = nsig_full
        out[f"{key}_sig_meas"] = sig_meas; out[f"{key}_nsigma_measonly"] = nsig_meas
    print("\nColumns: Nσ = vs the relation INCL. intrinsic scatter (the population dispersion — the "
          "correct test of single-object consistency); Nσ(+α,β) also folds the relation-parameter "
          "errors; Nσ(meas-only) = distance to the SCATTER-FREE mean line using only our measurement "
          "errors.\nInterpretation: both offsets are negative (mildly undermassive BH); INCL. scatter "
          "the point is <1σ → fully consistent with both local relations (scatter dominates). Against "
          "the scatter-free mean line it is ~1.8σ (measurement-limited). Registry mbh_sigma_offset=−0.69 "
          "added ε to the mean; the mean-relation offset is −0.26 (M•–σ) / −0.49 (M•–M⋆).")
    np.savez("results/relation_offset_significance.npz", **out)
    print("\nSaved -> results/relation_offset_significance.npz")


if __name__ == "__main__":
    main()
