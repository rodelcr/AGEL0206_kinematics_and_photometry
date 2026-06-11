"""Validate the (multi-start) Sérsic fitter by synthetic-source recovery (2026-06-11).

No external 2D-Sérsic fitter (GALFIT/imfit/petrofit/statmorph) is installed, so we
validate against astropy's `Sersic2D` — the published, peer-reviewed reference
implementation of the Sérsic (1968) profile (Astropy Collaboration 2013/2018/2022).
This is the same injection-recovery test statmorph (Rodriguez-Gomez+2019) and petrofit
(Geda+2022) use to validate themselves: render a KNOWN Sérsic galaxy, add realistic
noise, run `sersic_total_photometry.fit_sersic2d`, and check the recovered parameters
are unbiased — in particular that the fixed multi-start fitter recovers a known
ELLIPTICITY (the b/a-railing bug it was built to cure) rather than collapsing to b/a=1.

Grid: n∈{1,2,4} × b/a∈{0.70,0.85,1.00} × PA∈{20,70,120}° at fixed r_eff=1.8" and the
deflector's brightness, on the F150W2 pixel grid, at the band's empirical noise level.

Output: results/validate_sersic_fitter_synthetic.npz (+ printed bias/scatter table)
Usage: conda activate ISMgas; python scripts/validate_sersic_fitter_synthetic.py
"""
import os, sys, itertools
import numpy as np
from astropy.modeling.models import Sersic2D

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band
from scripts.sersic_total_photometry import (fit_sersic2d, sersic_total_flux_analytic,
                                             b_n, R_E_INIT)

BAND = "F150W2"      # bright deflector → ellipticity well-posed
REFF_AS = 1.8
N_GRID = (1.0, 2.0, 4.0)
BA_GRID = (0.70, 0.85, 1.00)
PA_GRID = (20.0, 70.0, 120.0)   # deg


def main():
    b = load_band(BAND)
    ps = b["pix_scale"]; zp = b["ab_zp"]
    ny, nx = 220, 220
    cx, cy = nx / 2.0, ny / 2.0
    reff_pix = REFF_AS / ps
    # noise level = per-pixel sky sigma from a BLANK corner of the real band (the
    # whole-image MAD is inflated by the bright deflector → use source-free corners)
    img = np.nan_to_num(b["img"].astype(float))
    h, w = img.shape; c = 40
    corners = np.concatenate([img[:c, :c].ravel(), img[:c, -c:].ravel(),
                              img[-c:, :c].ravel(), img[-c:, -c:].ravel()])
    sig = float(np.median(np.abs(corners - np.median(corners))) * 1.4826)
    # amplitude so the injected total mag ~ the real deflector (m≈18.2 in F150W2)
    m_target = 18.2
    F_target = 10 ** (-0.4 * (m_target - zp))

    rng = np.random.default_rng(20260611)
    yy, xx = np.mgrid[:ny, :nx]
    rows = []
    nomask = np.zeros((ny, nx), bool)
    for n, ba, pa in itertools.product(N_GRID, BA_GRID, PA_GRID):
        ell = 1 - ba; th = np.radians(pa)
        # set amplitude (=I_e) to hit F_target for this (n, ell)
        bn = b_n(n)
        unit = (2 * np.pi * n * np.exp(bn) / bn ** (2 * n)
                * 1.0 * reff_pix ** 2 * (1 - ell))   # F_tot per unit amplitude
        from scipy.special import gamma as G
        unit *= G(2 * n)
        amp = F_target / unit
        truth = Sersic2D(amplitude=amp, r_eff=reff_pix, n=n, x_0=cx, y_0=cy,
                         ellip=ell, theta=th)
        clean = np.asarray(truth(xx, yy), dtype=float)
        noisy = clean + rng.normal(0, sig, size=clean.shape)
        fit, _, info = fit_sersic2d(noisy, nomask, (cx, cy), R_E_INIT, 6.0, ps)
        rn, rreff, rell, rth = (float(fit.n.value), float(fit.r_eff.value) * ps,
                                float(fit.ellip.value), float(fit.theta.value))
        rba = 1 - rell
        rmag = -2.5 * np.log10(sersic_total_flux_analytic(
            float(fit.amplitude.value), float(fit.r_eff.value), rn, rell)) + zp
        # PA recovered mod 180, compare only if elliptical enough to be meaningful
        dpa = (np.degrees(rth) - pa + 90) % 180 - 90 if ba < 0.97 else np.nan
        rows.append((n, ba, pa, rn, rreff, rba, rmag - m_target, dpa))

    rows = np.array(rows, dtype=float)
    print(f"  Injection-recovery on the {BAND} grid (astropy Sersic2D truth → our fitter):")
    print(f"  {'n_in':>4s} {'ba_in':>5s} {'PA_in':>5s} | {'n_rec':>5s} {'re_rec':>6s} "
          f"{'ba_rec':>6s} {'Δmag':>6s} {'ΔPA':>6s}")
    for r in rows:
        print(f"  {r[0]:>4.0f} {r[1]:>5.2f} {r[2]:>5.0f} | {r[3]:>5.2f} {r[4]:>6.3f} "
              f"{r[5]:>6.2f} {r[6]:>+6.3f} {('%+5.0f'%r[7]) if np.isfinite(r[7]) else '  -- '}")

    # summary bias/scatter (exclude the round b/a=1 cases from the b/a-bias stat,
    # where the recovered ellipticity is expected to be ~0 by construction)
    ell_mask = rows[:, 1] < 0.97
    dn = rows[:, 3] - rows[:, 0]
    dba = rows[:, 5] - rows[:, 1]
    dmag = rows[:, 6]
    dpa = rows[ell_mask, 7]
    print("\n  recovery bias ± scatter:")
    print(f"    n     : {np.mean(dn):+.3f} ± {np.std(dn):.3f}")
    print(f"    b/a   : {np.mean(dba):+.3f} ± {np.std(dba):.3f}   (ellip-only: "
          f"{np.mean(dba[ell_mask]):+.3f} ± {np.std(dba[ell_mask]):.3f})")
    print(f"    m_tot : {np.mean(dmag):+.3f} ± {np.std(dmag):.3f} mag")
    print(f"    PA    : {np.nanmean(dpa):+.1f} ± {np.nanstd(dpa):.1f} deg  (elliptical inputs)")
    # the key test: are the b/a=0.70/0.85 inputs recovered as elliptical (not railed to 1)?
    rec_ell = rows[ell_mask, 5]
    n_railed = int(np.sum(rec_ell > 0.97))
    print(f"\n  ellipticity-railing check: {n_railed}/{ell_mask.sum()} elliptical inputs "
          f"mis-recovered as round (b/a>0.97)  → {'PASS' if n_railed == 0 else 'FAIL'}")

    np.savez("results/validate_sersic_fitter_synthetic.npz", rows=rows,
             cols=np.array(["n_in", "ba_in", "PA_in", "n_rec", "re_rec_as",
                            "ba_rec", "dmag", "dPA"]),
             n_bias=float(np.mean(dn)), n_scat=float(np.std(dn)),
             ba_bias=float(np.mean(dba[ell_mask])), ba_scat=float(np.std(dba[ell_mask])),
             mag_bias=float(np.mean(dmag)), mag_scat=float(np.std(dmag)),
             n_railed=n_railed, n_ell_inputs=int(ell_mask.sum()))
    print("\n  → results/validate_sersic_fitter_synthetic.npz")


if __name__ == "__main__":
    main()
