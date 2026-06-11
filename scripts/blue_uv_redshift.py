"""Deflector redshift from rest-UV resonance lines in the KCRM BLUE-arm cube (2026-06-11).

The red cube (z from Ca H/K, Balmer, G-band) is the headline; this adds the blue-arm
rest-UV absorption features the red cube can't reach: the **Mg II λλ2796,2803 doublet**
and the **Fe II UV multiplet** (2344/2374/2383/2587/2600). At z=0.676 these land at
obs 3.9-4.7 kµm — inside the blue cube (obs 3278-5872 Å).

AIR/VAC handling (critical): the BLUE cube is **air** (FITS CTYPE3='AWAV'), unlike the
RED cube which is **vacuum** (CTYPE3='WAVE'). UV resonance lines are standardly quoted
in VACUUM (Morton 2003), so we convert the blue cube air→vacuum (Ciddor 1996 index) and
fit vacuum rest wavelengths: z = λ_obs,vac/λ_rest,vac − 1. (nb04's blue cells mislabeled
the air cube as '_vac' and applied vac_to_air → a latent double-conversion; avoided here.)

Output: results/blue_uv_redshift.npz (+ printed per-line table)
Usage: conda activate ISMgas; python scripts/blue_uv_redshift.py
"""
import os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit

RA_DEFL, DEC_DEFL = 31.55611, -1.23817

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)

BLUE = "DESJ0206_medium_combined_blue.fits"
Z_GUESS = 0.67564

# Morton (2003) VACUUM rest wavelengths (Å) of the deflector rest-UV features.
MGII = (2796.352, 2803.531)            # Mg II resonance doublet
FEII = {"Fe II 2344": 2344.214, "Fe II 2374": 2374.461, "Fe II 2383": 2382.765,
        "Fe II 2587": 2586.650, "Fe II 2600": 2600.173}


def ciddor_n(lam_air):
    s = 1e4 / lam_air
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s**2) + 1.5998e-4 / (38.9 - s**2)


def air_to_vac(lam_air):
    # one Newton step is plenty at these wavelengths (n varies slowly)
    return lam_air * ciddor_n(lam_air)


def gauss_abs(x, amp, cen, sig, cont):
    return cont - amp * np.exp(-0.5 * ((x - cen) / sig) ** 2)


def mgii_doublet(x, amp1, amp2, z, sig, cont):
    c1, c2 = MGII[0] * (1 + z), MGII[1] * (1 + z)
    return cont - amp1 * np.exp(-0.5 * ((x - c1) / sig) ** 2) \
                - amp2 * np.exp(-0.5 * ((x - c2) / sig) ** 2)


def main():
    with fits.open(BLUE) as h:
        hdr = h[0].header; cube = h[0].data.astype(float)
    crval = hdr["CRVAL3"]; cd = hdr.get("CD3_3", hdr.get("CDELT3"))
    crpix = hdr.get("CRPIX3", 1.0); npix = hdr["NAXIS3"]
    lam_air = crval + cd * (np.arange(npix) + 1 - crpix)
    lam = air_to_vac(lam_air)              # → vacuum
    print(f"  blue cube CTYPE3={hdr.get('CTYPE3')} (air) → converted to vacuum; "
          f"obs {lam[0]:.0f}-{lam[-1]:.0f} Å vac")

    # deflector located by its KNOWN sky position via the cube WCS (the blue cube is
    # dominated by the blue lensed arc, so a white-light peak would grab the arc, not the
    # red deflector). Refine to the local continuum peak within ±2 spaxels.
    ny, nx = cube.shape[1], cube.shape[2]
    wx, wy = WCS(hdr, naxis=2).world_to_pixel_values(RA_DEFL, DEC_DEFL)
    cx0, cy0 = int(round(float(wx))), int(round(float(wy)))
    white = np.nansum(cube, axis=0)
    sub = white[cy0-2:cy0+3, cx0-2:cx0+3]
    dyx = np.unravel_index(np.nanargmax(sub), sub.shape)
    cy, cx = cy0 - 2 + int(dyx[0]), cx0 - 2 + int(dyx[1])
    yy, xx = np.mgrid[:ny, :nx]
    aper = np.hypot(xx - cx, yy - cy) <= 3        # 0.9" radius
    print(f"  deflector WCS pixel=({cx0},{cy0}) → continuum-refined ({cx},{cy})")
    sky = np.hypot(xx - cx, yy - cy)
    skymask = (sky > 12) & (sky < 20)
    flux = np.nansum(cube[:, aper], axis=1)
    noise = np.nanstd(cube[:, skymask], axis=1) * np.sqrt(aper.sum())
    cont = np.nanmedian(flux[(lam > MGII[0]*(1+Z_GUESS)-200) & (lam < MGII[0]*(1+Z_GUESS)+200)])
    sn = cont / np.nanmedian(noise[(lam > 4500) & (lam < 4900)])
    print(f"  deflector centre (y,x)=({cy},{cx})  aper=29 spax  continuum S/N≈{sn:.1f} at Mg II")

    # --- JOINT multi-line fit: ALL UV lines share one (z, σ); per-line amplitude ≥0.
    # The robust low-S/N approach — independent fits scatter on noise. Continuum-normalise
    # first (median filter), then fit normalised = 1 − Σ_i A_i exp(−½((λ−rest_i(1+z))/σ_λ)²). ---
    from scipy.ndimage import median_filter
    allrest = list(MGII) + [FEII[k] for k in FEII]
    lo, hi = min(allrest) * (1 + Z_GUESS) - 60, max(allrest) * (1 + Z_GUESS) + 60
    reg = (lam > lo) & (lam < hi)
    lr, fr, nr = lam[reg], flux[reg], noise[reg]
    cont_sm = median_filter(fr, size=151)
    norm = fr / cont_sm; nnorm = nr / cont_sm

    def joint(x, z, sigA, *amps):
        m = np.ones_like(x)
        for a, r0 in zip(amps, allrest):
            m = m - a * np.exp(-0.5 * ((x - r0 * (1 + z)) / sigA) ** 2)
        return m
    try:
        p0 = [Z_GUESS, 4.0] + [0.05] * len(allrest)
        lb = [Z_GUESS - 0.005, 1.0] + [0.0] * len(allrest)
        ub = [Z_GUESS + 0.005, 12.0] + [0.8] * len(allrest)
        jpopt, jpcov = curve_fit(joint, lr, norm, p0=p0, sigma=nnorm,
                                 absolute_sigma=True, bounds=(lb, ub), maxfev=40000)
        z_joint, z_joint_err = jpopt[0], np.sqrt(np.diag(jpcov))[0]
        det = [(allrest[i], jpopt[2 + i], np.sqrt(np.diag(jpcov))[2 + i]) for i in range(len(allrest))]
        print(f"\n  JOINT multi-line fit (shared z, {len(allrest)} lines):")
        print(f"    z = {z_joint:.5f} ± {z_joint_err:.5f}   σ_λ = {jpopt[1]:.2f} Å")
        for r0, a, ae in det:
            print(f"    rest {r0:8.2f}:  depth {a:.3f} ± {ae:.3f}  ({a/ae:.1f}σ)" if ae > 0 else
                  f"    rest {r0:8.2f}:  depth {a:.3f}")
        ndet = sum(1 for _, a, ae in det if ae > 0 and a / ae > 2)
        print(f"    → {ndet}/{len(allrest)} lines >2σ; cz offset from red z = "
              f"{(z_joint-Z_GUESS)*3e5:.0f} km/s")
    except Exception as e:
        z_joint, z_joint_err, ndet = np.nan, np.nan, 0
        print("  joint fit failed:", e)

    rows = []
    # --- Mg II doublet (shared z) ---
    win = (lam > MGII[0]*(1+Z_GUESS) - 40) & (lam < MGII[1]*(1+Z_GUESS) + 40)
    try:
        c0 = np.nanmedian(flux[win])
        p0 = [0.1*c0, 0.07*c0, Z_GUESS, 3.0, c0]
        popt, pcov = curve_fit(mgii_doublet, lam[win], flux[win], p0=p0,
                               sigma=noise[win], absolute_sigma=True, maxfev=20000)
        zerr = np.sqrt(np.diag(pcov))[2]
        rows.append(("Mg II 2796/2803", "doublet", popt[2], zerr,
                     popt[0]/np.nanmedian(noise[win])))
    except Exception as e:
        rows.append(("Mg II 2796/2803", "doublet", np.nan, np.nan, 0)); print("  MgII fit:", e)

    # --- Fe II single lines ---
    for nm, lr in FEII.items():
        obs = lr * (1 + Z_GUESS)
        w = (lam > obs - 25) & (lam < obs + 25)
        if w.sum() < 6:
            continue
        try:
            c0 = np.nanmedian(flux[w])
            popt, pcov = curve_fit(gauss_abs, lam[w], flux[w],
                                   p0=[0.1*c0, obs, 2.5, c0], sigma=noise[w],
                                   absolute_sigma=True, maxfev=20000)
            z = popt[1]/lr - 1; zerr = np.sqrt(np.diag(pcov))[1]/lr
            amp_snr = popt[0]/np.nanmedian(noise[w])
            rows.append((nm, f"{lr:.2f}", z, zerr, amp_snr))
        except Exception as e:
            print(f"  {nm} fit:", e)

    print(f"\n  {'feature':18s} {'rest_vac':>9s} {'z':>8s} {'z_err':>8s} {'amp/σ':>6s}")
    good = []
    for nm, rl, z, ze, snr in rows:
        flag = "" if (snr > 2 and np.isfinite(z) and abs(z-Z_GUESS) < 0.01) else "  (weak/unreliable)"
        print(f"  {nm:18s} {rl:>9s} {z:>8.5f} {ze:>8.5f} {snr:>6.1f}{flag}")
        if not flag:
            good.append((z, ze))
    if good:
        zs = np.array([g[0] for g in good]); ws = 1/np.array([g[1] for g in good])**2
        zc = np.sum(zs*ws)/np.sum(ws); zce = 1/np.sqrt(np.sum(ws))
        print(f"\n  Blue-arm UV weighted z = {zc:.5f} ± {zce:.5f}  ({len(good)} reliable lines)")
        print(f"  (red-cube headline z = 0.67564 ± 0.00033)")
    else:
        zc, zce = np.nan, np.nan
        print("\n  No reliable UV detection — deflector too faint in rest-UV (report as attempted).")

    np.savez("results/blue_uv_redshift.npz",
             rows=np.array(rows, dtype=object), z_blue_uv=zc, z_blue_uv_err=zce,
             sn_continuum=sn, center_yx=np.array([cy, cx]))
    print("\n  → results/blue_uv_redshift.npz")


if __name__ == "__main__":
    main()
