"""Lensed-source (arc) redshift from rest-UV ISM resonance lines (2026-06-11).

The z=1.302 lensed arc is a bright star-forming galaxy → strong **Mg II λλ2796,2803**
and **Fe II UV** resonance ABSORPTION (the classic SF-galaxy ISM lines). At z=1.302:
  - Fe II 2344/2374/2382  → obs 5396/5466/5485 Å  → in the BLUE cube (3279-5874 Å, air)
  - Mg II 2796/2803       → obs 6442/6459 Å        → in the RED cube  (5142-9416 Å, vacuum)
So the "Fe triplet" comes from the blue arm, the "Mg doublet" from the red arm.

AIR/VAC: rest UV lines are vacuum (Morton 2003). The blue cube is air (AWAV) → convert
to vacuum; the red cube is already vacuum (WAVE). z = λ_obs,vac/λ_rest,vac − 1.

Arc extraction: the arc is the brightest BLUE continuum source (≠ the red deflector); in
the red cube it is located from its [O II] λ3727 (z=1.302 → 8580 Å) emission map.

Output: results/source_uv_redshift.npz (+ printed table)
Usage: conda activate ISMgas; python scripts/source_uv_redshift.py
"""
import os, sys
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)

BLUE = "DESJ0206_medium_combined_blue.fits"
RED = "raw_KCWI/New_red/Nov17_2025_DESJ0206_RL_combined_mtwdo_icubes_wcs.fits"
Z_SRC = 1.302
FE_TRIPLET = {"Fe II 2344": 2344.214, "Fe II 2374": 2374.461, "Fe II 2382": 2382.765}
MG_DOUBLET = (2796.352, 2803.531)


def ciddor_n(la):
    s = 1e4 / la
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s**2) + 1.5998e-4 / (38.9 - s**2)


def lam_axis(hdr, to_vac):
    cd = hdr.get("CD3_3", hdr.get("CDELT3"))
    lam = hdr["CRVAL3"] + cd * (np.arange(hdr["NAXIS3"]) + 1 - hdr.get("CRPIX3", 1.0))
    return lam * ciddor_n(lam) if to_vac else lam       # air→vac if blue


def gauss_abs(x, amp, cen, sig, cont):
    return cont - amp * np.exp(-0.5 * ((x - cen) / sig) ** 2)


def make_joint(rests):
    """Joint absorption model with FIXED rest wavelengths (only z, σ, cont, per-line
    amplitudes free → a well-posed z error, unlike freeing the rests)."""
    def joint_abs(x, z, sig, cont, *amps):
        m = np.full_like(x, cont)
        for a, r0 in zip(amps, rests):
            m = m - a * np.exp(-0.5 * ((x - r0 * (1 + z)) / sig) ** 2)
        return m
    return joint_abs


def extract_blue_arc():
    with fits.open(BLUE) as h:
        hdr, cube = h[0].header, h[0].data.astype(float)
    lam = lam_axis(hdr, to_vac=True)
    white = np.nansum(cube, axis=0)
    cy, cx = np.unravel_index(np.nanargmax(white), white.shape)
    yy, xx = np.mgrid[:white.shape[0], :white.shape[1]]
    aper = np.hypot(xx - cx, yy - cy) <= 3
    sky = (np.hypot(xx - cx, yy - cy) > 12) & (np.hypot(xx - cx, yy - cy) < 20)
    flux = np.nansum(cube[:, aper], axis=1)
    noise = np.nanstd(cube[:, sky], axis=1) * np.sqrt(aper.sum())
    print(f"  BLUE arc at white-light peak (y,x)=({cy},{cx})")
    return lam, flux, noise


def extract_red_arc():
    with fits.open(RED) as h:
        hdr, cube = h[0].header, h[0].data.astype(float)
    lam = lam_axis(hdr, to_vac=False)
    # [O II] 3727 (z=1.302 → 8580 vac) narrow-band → locate the arc
    nb = (lam > 8575) & (lam < 8595)
    cont = (lam > 8500) & (lam < 8560)
    oii = np.nansum(cube[nb], axis=0) - np.nanmedian(cube[cont], axis=0) * nb.sum()
    cy, cx = np.unravel_index(np.nanargmax(np.nan_to_num(oii)), oii.shape)
    yy, xx = np.mgrid[:oii.shape[0], :oii.shape[1]]
    aper = np.hypot(xx - cx, yy - cy) <= 3
    sky = (np.hypot(xx - cx, yy - cy) > 12) & (np.hypot(xx - cx, yy - cy) < 20)
    flux = np.nansum(cube[:, aper], axis=1)
    noise = np.nanstd(cube[:, sky], axis=1) * np.sqrt(aper.sum())
    print(f"  RED arc at [O II] peak (y,x)=({cy},{cx})")
    return lam, flux, noise


def fit_lines(lam, flux, noise, rests, label):
    lo, hi = min(rests.values()) * (1 + Z_SRC) - 60, max(rests.values()) * (1 + Z_SRC) + 60
    reg = (lam > lo) & (lam < hi) & np.isfinite(flux)
    lr, fr, nr = lam[reg], flux[reg], noise[reg]
    cont = median_filter(fr, size=121); norm = fr / cont; nn = nr / cont
    names, r0s = list(rests.keys()), list(rests.values())
    joint_abs = make_joint(r0s)
    p0 = [Z_SRC, 4.0, 1.0] + [0.1] * len(r0s)
    lb = [Z_SRC - 0.006, 1.0, 0.7] + [0.0] * len(r0s)
    ub = [Z_SRC + 0.006, 14.0, 1.3] + [0.9] * len(r0s)
    try:
        popt, pcov = curve_fit(joint_abs, lr, norm, p0=p0, sigma=nn, absolute_sigma=True,
                               bounds=(lb, ub), maxfev=60000)
        z, ze = popt[0], np.sqrt(np.diag(pcov))[0]
        amps = popt[3:3 + len(r0s)]; aerr = np.sqrt(np.diag(pcov))[3:3 + len(r0s)]
        print(f"\n  {label}: joint z = {z:.5f} ± {ze:.5f}  (σ_λ={popt[1]:.1f} Å)")
        det = []
        for nm, a, ae in zip(names, amps, aerr):
            snr = a / ae if ae > 0 else 0
            print(f"    {nm:12s} depth {a:.3f} ± {ae:.3f}  ({snr:.1f}σ)")
            det.append(snr)
        return z, ze, max(det), popt[1]
    except Exception as e:
        print(f"  {label}: fit failed — {e}")
        return np.nan, np.nan, 0, np.nan


def main():
    lamb, fb, nb_ = extract_blue_arc()
    z_fe, ze_fe, snr_fe, sig_fe = fit_lines(lamb, fb, nb_, FE_TRIPLET, "Fe II triplet (blue arm)")

    lamr, fr, nr = extract_red_arc()
    mg = {"Mg II 2796": MG_DOUBLET[0], "Mg II 2803": MG_DOUBLET[1]}
    z_mg, ze_mg, snr_mg, sig_mg = fit_lines(lamr, fr, nr, mg, "Mg II doublet (red arm)")

    print(f"\n  ── Source (arc) UV redshift summary ──")
    print(f"    Fe II triplet (blue): z = {z_fe:.5f} ± {ze_fe:.5f}  (max {snr_fe:.1f}σ)")
    print(f"    Mg II doublet (red):  z = {z_mg:.5f} ± {ze_mg:.5f}  (max {snr_mg:.1f}σ)")
    print(f"    [O II] doublet (red, nb04): z = 1.30263 ± 0.00003  (systemic)")
    print(f"    NB: a Mg II–[O II] velocity offset is NOT necessarily an outflow — emission")
    print(f"        infill, differential lensing/aperture, and profile asymmetry all mimic it.")
    np.savez("results/source_uv_redshift.npz",
             z_fe=z_fe, z_fe_err=ze_fe, snr_fe=snr_fe, sig_fe=sig_fe,
             z_mg=z_mg, z_mg_err=ze_mg, snr_mg=snr_mg, sig_mg=sig_mg, z_oii=1.30263)
    print("\n  → results/source_uv_redshift.npz")


if __name__ == "__main__":
    main()
