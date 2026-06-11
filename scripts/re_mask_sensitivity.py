import os, sys, numpy as np
from astropy.io import fits
sys.path.insert(0, os.getcwd())
from scripts.final_sigma_e import curve_of_growth, find_center, RA_DEFL, DEC_DEFL
from scripts.mask_method_comparison import load_band
from astropy.wcs import WCS
RE_HEAD = 2.305
sysz = np.load("results/photometry_systematics.npz", allow_pickle=True)
apm = np.load("results/aperture_2re_masks.npz", allow_pickle=True)
VDI = "../velocity_dispersion_from_IFU"
cfg = {"F140W": f"{VDI}/AGEL020613-011417A_F140W_WFC3_cutout_L3.fits",
       "F200LP": f"{VDI}/AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits"}
mcfg = {k: v.replace(".fits", "_mask.fits") for k, v in cfg.items()}

def re_for(masks):
    res = {}
    for n in ["F140W", "F200LP"]:
        with fits.open(cfg[n]) as h:
            img = h[0].data.astype(float); ps = abs(WCS(h[0].header).proj_plane_pixel_scales()[0].value) * 3600
        w = WCS(fits.getheader(cfg[n]))
        xc, yc = find_center(img, masks[n], w, RA_DEFL, DEC_DEFL, 3.0, ps)
        res[n] = curve_of_growth(img, (xc, yc), ps, mask=masks[n])
    return res, 0.5 * (res["F140W"] + res["F200LP"])

expert = {n: fits.getdata(mcfg[n]).astype(bool) for n in cfg}
glob = {n: sysz[f"{n}_global_mask"].astype(bool) for n in cfg}
globcomp = {n: (sysz[f"{n}_global_mask"].astype(bool) | apm[f"{n}_2Re_companion"].astype(bool)) for n in cfg}

for lab, masks in [("expert (headline)", expert), ("new global (color/morph+reg)", glob),
                   ("new global + companion", globcomp)]:
    r, mean = re_for(masks)
    print(f"  {lab:<30} F140W={r['F140W']:.3f}\"  F200LP={r['F200LP']:.3f}\"  mean={mean:.3f}\"  (Δ vs 2.305 = {mean-RE_HEAD:+.3f}\")")
