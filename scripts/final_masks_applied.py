"""The FINAL arc masks applied in the deflector aperture photometry, all 4 bands.

Single-Sersic deflector model + band-specific IR-extension gate (HST color gate /
deep-JWST morphology guard) + WCS cross-correlation registration. Shows, per band:
the image (asinh), the per-band mask (lime fill+contour), the global mask (orange
contour), and the photometry aperture (cyan dashed) + sky annulus (cyan dotted).
Annotates mask pixel count and the raw aperture AB mag actually used.

Usage: conda activate ISMgas; python scripts/final_masks_applied.py
Output: results/figures/final_masks_applied.png
"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.visualization import AsinhStretch, ManualInterval

REPO = "/Users/rosador/Documents/AGEL/AGEL0206_kinematics_and_photometry"
sys.path.insert(0, REPO); os.chdir(REPO)
from scripts.mask_method_comparison import load_band

ORDER = ["F200LP", "F140W", "F150W2", "F322W2"]
C = np.load("results/photometry_systematics.npz", allow_pickle=True)


def main():
    bands = {n: load_band(n) for n in ORDER}
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.6))
    for ax, n in zip(axes, ORDER):
        b = bands[n]; ps = b["pix_scale"]
        p = json.load(open(b["cfg"]["params"])); pg = p["pixel_geometry"]
        xa, ya = b["wcs"].world_to_pixel_values(p["target_ra"], p["target_dec"])
        half = int(np.ceil(4.5 / ps))
        cy, cx = int(ya), int(xa)
        sl = (slice(max(0, cy - half), cy + half), slice(max(0, cx - half), cx + half))
        y0, x0 = sl[0].start, sl[1].start
        img = np.nan_to_num(b["img"])[sl]
        lo, hi = np.nanpercentile(img, 35), np.nanpercentile(img, 99.2)
        ax.imshow(AsinhStretch(0.05)(ManualInterval(lo, hi)(img)), origin="lower", cmap="gray")
        pb = C[f"{n}_perband_mask"][sl].astype(float)
        gl = C[f"{n}_global_mask"][sl].astype(float)
        ax.contourf(pb, levels=[.5, 1.5], colors="lime", alpha=0.30)
        ax.contour(pb, levels=[.5], colors="lime", linewidths=1.2)
        ax.contour(gl, levels=[.5], colors="orange", linewidths=0.8)
        # aperture + sky annulus
        ax.add_patch(Ellipse((xa - x0, ya - y0), 2 * pg["a"], 2 * pg["b"],
                             angle=np.degrees(pg["theta_rad"]), fill=False,
                             ec="cyan", lw=1.5, ls="--"))
        for aa in ("annulus_a_in", "annulus_a_out"):
            ratio = pg["b"] / pg["a"]
            ax.add_patch(Ellipse((xa - x0, ya - y0), 2 * pg[aa], 2 * pg[aa] * ratio,
                                 angle=np.degrees(pg["theta_rad"]), fill=False,
                                 ec="cyan", lw=0.8, ls=":"))
        i = ORDER.index(n)
        ax.set_title(f"{n}: per-band mask={int(C[f'{n}_perband_mask'].sum())}px  "
                     f"raw={C['perband_raw'][i]:.2f} AB\n"
                     f"lime=mask  orange=global  cyan--=aperture  cyan:=sky annulus",
                     fontsize=9)
        ax.set_xlim(0, img.shape[1]); ax.set_ylim(0, img.shape[0])  # keep image filling panel
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("FINAL arc masks applied in the deflector aperture photometry "
                 "(single-Sersic + color/morph gate + WCS-registered)", fontsize=13, y=1.02)
    plt.tight_layout()
    out = "results/figures/final_masks_applied.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
