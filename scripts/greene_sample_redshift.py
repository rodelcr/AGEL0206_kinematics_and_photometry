"""Redshift reach of the Greene+2020 calibrating sample (2026-06-17).

Backs §5 of NOTES_relation_offset_provenance_2026-06-17.md: the most distant galaxy in the
Greene+2020 M•–σ and M•–M⋆ "Early" samples, to show how local the calibration is relative to
AGEL J0206 (z=0.67564). The supplemental tables list distances (Mpc), not redshifts; we convert
with z ≈ H0·D/c (H0=67.7, Planck 2015) — exact only in the Hubble-flow limit (fine for the
farthest members, a recession-velocity proxy for the nearest).

Samples assembled exactly as in the paper figure code (figures_paper4):
  M•–σ  = Suppl-2 (HT in {E,S0}) + Suppl-3 (all) + Suppl-4 (HT=S0)
  M•–M⋆ = Suppl-3 (HT=E) + Suppl-4 (HT=S0), finite logMstar

Run from repo root:  python -m scripts.greene_sample_redshift
"""
import numpy as np
from astropy.table import Table, vstack

C_KMS = 299792.458
H0 = 67.7
MBH_DIR = "../AGEL_Mbh_sigma"


def z_of(D_mpc):
    return H0 * np.asarray(D_mpc, float) / C_KMS


def main():
    s2 = Table.read(f"{MBH_DIR}/Greene20_Supple_2.csv")
    s3 = Table.read(f"{MBH_DIR}/Greene20_Supple_3.csv")
    s4 = Table.read(f"{MBH_DIR}/Greene20_Supple_4.csv")

    msig = vstack([s2[(s2["HT"] == "E") | (s2["HT"] == "S0")], s3, s4[s4["HT"] == "S0"]])
    mst = vstack([s3[s3["HT"] == "E"], s4[s4["HT"] == "S0"]])
    lm = np.array([x if x not in ("", None) else np.nan for x in mst["logMstar"]], float)
    mst = mst[np.isfinite(lm)]

    for label, tab in [("M•–σ", msig), ("M•–M⋆", mst)]:
        D = np.array(tab["Distance_Mpc"], float)
        order = np.argsort(D)[::-1]
        print(f"\n=== {label}  (N = {len(tab)}) ===  most distant 5:")
        for j in order[:5]:
            print(f"  {tab['Galaxy'][j]:14s} {D[j]:7.1f} Mpc   z ≈ {z_of(D[j]):.4f}")
        i = order[0]
        print(f"  -> MAX: {tab['Galaxy'][i]} = {D[i]:.1f} Mpc, z ≈ {z_of(D[i]):.4f}")

    allg = vstack([s2, s3, s4])
    Da = np.array(allg["Distance_Mpc"], float); k = np.nanargmax(Da)
    print(f"\nfull compilation (N={len(allg)}): max {allg['Galaxy'][k]} "
          f"{Da[k]:.1f} Mpc, z ≈ {z_of(Da[k]):.4f}")
    print(f"\nAGEL J0206 is at z = 0.67564 — ~{0.67564/z_of(152.4):.0f}× the redshift of the "
          f"most distant calibrator.")


if __name__ == "__main__":
    main()
