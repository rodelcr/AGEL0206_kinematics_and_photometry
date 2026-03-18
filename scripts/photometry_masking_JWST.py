"""
Interactive JWST Aperture Photometry Tool

Performs elliptical aperture photometry on JWST imaging data with:
- Interactive GUI for aperture/annulus geometry adjustment (matplotlib sliders)
- Click-to-mask contaminating sources (left-click add, right-click undo)
- Noise map support (ERR extension or external RMS maps) with error propagation
- AB magnitude and flux density (Jy, erg/s/cm2/Hz) output
- Save/load parameter files (JSON) and mask files (FITS)
- Import and reproject masks/params from other filters

Three execution modes:
  1. Reload: Load existing params.json + mask.fits and recalculate
  2. Import: Import parameters/masks from another filter and reproject
  3. Interactive: Launch GUI to define aperture and mask from scratch

Usage:
  python photometry_masking_JWST.py
  (Edit the 'fname' variable near the top to select the target FITS file)

JWST-specific notes:
  - Uses PIXAR_SR header keyword (pixel area in steradians) for AB zero point
  - AB_ZP = -6.10 - 2.5*log10(PIXAR_SR)
  - Supports JWST multi-extension files (SCI in ext 1, ERR in ext 2)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import json
import sys
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils.aperture import EllipticalAperture, EllipticalAnnulus, CircularAperture, ApertureStats

# We need reproject to align the large noise map to the small cutout
try:
    from reproject import reproject_interp
    HAS_REPROJECT = True
except ImportError:
    print("WARNING: 'reproject' package not found.")
    print("You will not be able to import masks or align noise maps.")
    HAS_REPROJECT = False

# --- 1. SETTINGS ---
#fname = 'AGEL020613-011417A_F200LP_WFC3_cutout_L3.fits' 
#fname = 'AGEL020613-011417A_F140W_WFC3_cutout_L3.fits'
fname = 'jw05594-o101_t103_nircam_clear-f150w2_i2d.fits'
#fname = 'jw05594-o101_t103_nircam_clear-f322w2_i2d.fits'

filename = fname
print(f"Target file: {filename}")

target_ra = 31.55611000 + 0.00005
target_dec = -1.23817000 + 0.00005

# Initial Defaults
init_semimajor_axis = 0.6  
init_ellipticity = 0.2     
init_theta_deg = 45.0      
init_annulus_gap = 0.2     
init_annulus_width = 0.4   
init_mask_radius = 0.2     

# Derived Filenames
json_filename = filename.replace('.fits', '_params.json')
mask_filename = filename.replace('.fits', '_mask.fits')

# -----------------------------------------------------------------------------
# HELPER: PHOTOMETRY MATH
# -----------------------------------------------------------------------------
def run_photometry_math(data, header, aperture, annulus, total_mask, error_map=None):
    """ 
    Performs statistics, noise calculation, and magnitude conversion. 
    error_map: A sigma (RMS) map matching 'data' dimensions.
    """
    
    # 1. Calculate Stats (Masked pixels ignored)
    # We pass the error map here so photutils calculates sum_err (propagation of pixel errors)
    ap_stats = ApertureStats(data, aperture, error=error_map, mask=total_mask)
    ann_stats = ApertureStats(data, annulus, error=error_map, mask=total_mask)

    raw_flux = ap_stats.sum
    # This is the propagated error of the pixels inside the aperture (sqrt(sum(sigma^2)))
    raw_flux_err = ap_stats.sum_err if error_map is not None else 0.0
    
    bkg_median = ann_stats.median
    eff_area = ap_stats.sum_aper_area.value
    
    # 2. Net Flux
    net_flux = raw_flux - (bkg_median * eff_area)

    # 3. Net Error Calculation
    # Total Error^2 = (Source Error)^2 + (Background Error)^2
    # Background Error involves the uncertainty in the *determination* of the median background.
    # Formula: sigma_bkg_total = Area * (sigma_bkg_pixel / sqrt(N_annulus_pixels))
    
    if error_map is not None:
        n_ann_pix = ann_stats.sum_aper_area.value
        
        # --- FIX: Manually extract noise values within annulus ---
        # ApertureStats doesn't calculate 'median error value' automatically.
        ann_mask = annulus.to_mask(method='center')
        
        if ann_mask is not None:
            # Get error values and bad-pixel mask values within the annulus geometry
            ann_error_values = ann_mask.get_values(error_map)
            ann_bad_pixel_values = ann_mask.get_values(total_mask)
            
            # Filter: Keep only error values where the mask is False (Good Data)
            # Flatten in case of dimension weirdness, though get_values returns 1D usually
            valid_error_values = ann_error_values[~ann_bad_pixel_values.astype(bool)]
            
            # Take the median of the noise values in the annulus
            if len(valid_error_values) > 0:
                sigma_bkg_pix = np.nanmedian(valid_error_values)
            else:
                sigma_bkg_pix = 0.0
        else:
            sigma_bkg_pix = 0.0
        # ---------------------------------------------------------

        bkg_determination_err = eff_area * (sigma_bkg_pix / np.sqrt(n_ann_pix))
        total_flux_err = np.sqrt(raw_flux_err**2 + bkg_determination_err**2)
    else:
        total_flux_err = 0.0

    print("\n--- RESULTS ---")
    print(f"Raw Flux (e-/s):    {raw_flux:.4f}")
    print(f"Background (med):   {bkg_median:.5f}")
    print(f"Effective Area:     {eff_area:.2f} pix")
    
    if error_map is not None:
        print(f"Net Flux (e-/s):    {net_flux:.4f} +/- {total_flux_err:.4f}")
        if total_flux_err > 0:
            snr = net_flux / total_flux_err
            print(f"S/N Ratio:          {snr:.2f}")
    else:
        print(f"Net Flux (e-/s):    {net_flux:.4f}")

    # 4. Retrieve Photometry Keywords
    pixar_sr = None    
    try:
        pixar_sr = header['PIXAR_SR']

    except KeyError:
        print("\nWARNING: PIXAR_SR keyword missing in current file header.")
        ans = input("Do you want to load these keywords from another FITS file? [y/n]: ").strip().lower()

        if ans == 'y':
            alt_path = input("Enter path to FITS file with valid header: ").strip().replace("'", "").replace('"', "")
            try:
                with fits.open(alt_path) as hdul_alt:
                    found = False
                    for ext_id in [1, 0]:
                        if len(hdul_alt) > ext_id and 'PIXAR_SR' in hdul_alt[ext_id].header:
                            pixar_sr = hdul_alt[ext_id].header['PIXAR_SR']
                            found = True
                            break
            except Exception:
                pass

    # 5. Calculate Magnitude and Flux Units
    if pixar_sr is not None:
        # AB Zero Point Formula
        ab_zp = -6.10 - (2.5 * np.log10(pixar_sr))
        
        # References for conversion from STSci documentation:
        #https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#gsc.tab=0

        if net_flux > 0:
            ab_mag = ab_zp - 2.5 * np.log10(net_flux)

            # Setup for Flux conversion (Jy and erg/s/cm^2/Hz)
            fnu_units = u.erg/u.s/u.Hz/(u.cm**2)
            Fnu = (ab_mag * u.ABmag).to(fnu_units)
            Janskys = Fnu.to(u.Jy)
            
            

            # Mag Error ~ (2.5 / np.log(10)) * (FluxErr / Flux)
            if total_flux_err > 0:
                mag_err = (2.5 / np.log(10)) * (total_flux_err / net_flux)
                mag_str = f"{ab_mag:.4f} +/- {mag_err:.4f}"

                # Propagate error to Flux density
                # Upper/Lower bounds in flux space (since mag is log)
                fnu_err_upper = ((ab_mag - mag_err) * u.ABmag).to(fnu_units) - Fnu
                fnu_err_lower = Fnu - ((ab_mag + mag_err) * u.ABmag).to(fnu_units)
                
                # Note: Magnitude error adds/subtracts inversely to flux
                # Lower Mag = Higher Flux
                
                fnu_str = f"{Fnu.value:.3e} (+{fnu_err_upper.value:.3e}/-{fnu_err_lower.value:.3e}) {fnu_units}"
            else:
                mag_str = f"{ab_mag:.4f}"
                fnu_str = f"{Fnu.value:.3e} {fnu_units}"
                
            print("--------------------------------")
            print(f"AB Zero Point: {ab_zp:.3f}")
            print(f"AB Magnitude:  {mag_str}")
            print(f"Flux Density:  {fnu_str}")
            print(f"Janskys:       {Janskys.value:.3e} Jy")
            print("--------------------------------")
        else:
            print("Object too faint (negative flux).")

# -----------------------------------------------------------------------------
# INTERACTIVE CLASS
# -----------------------------------------------------------------------------
class PhotometryTool:
    def __init__(self, data, wcs, x_t, y_t, pixel_scale_deg, existing_mask=None):
        self.data = data
        self.wcs = wcs
        self.x_t = x_t
        self.y_t = y_t
        self.pix_scale_deg = pixel_scale_deg
        self.base_mask = existing_mask if existing_mask is not None else np.zeros_like(data, dtype=bool)
        self.mask_data = [] 
        self.finished = False

        self.fig = plt.figure(figsize=(10, 12))
        self.ax = self.fig.add_axes([0.1, 0.35, 0.8, 0.6])
        
        mean_val = np.nanmean(data)
        std_val = np.nanstd(data)
        self.ax.imshow(data, cmap='Greys', origin='lower', 
                       vmin=mean_val - 1*std_val, vmax=mean_val + 5*std_val)
        
        if np.any(self.base_mask):
            self.ax.imshow(self.base_mask, cmap='Blues', origin='lower', alpha=0.3)

        self.ax.set_title("Left Click: Mask (Red) | Right Click: Undo | Blue: Imported Mask")

        # Sliders
        ax_a = self.fig.add_axes([0.15, 0.25, 0.65, 0.03])
        ax_ell = self.fig.add_axes([0.15, 0.21, 0.65, 0.03])
        ax_theta = self.fig.add_axes([0.15, 0.17, 0.65, 0.03])
        ax_gap = self.fig.add_axes([0.15, 0.13, 0.65, 0.03])
        ax_wid = self.fig.add_axes([0.15, 0.09, 0.65, 0.03])
        ax_mask = self.fig.add_axes([0.15, 0.05, 0.65, 0.03])
        
        self.s_a = Slider(ax_a, 'Semi-Major (")', 0.1, 5.0, valinit=init_semimajor_axis)
        self.s_ell = Slider(ax_ell, 'Ellipticity', 0.0, 0.9, valinit=init_ellipticity)
        self.s_theta = Slider(ax_theta, 'Angle (deg)', 0.0, 360.0, valinit=init_theta_deg)
        self.s_gap = Slider(ax_gap, 'Sky Gap (")', 0.0, 12.0, valinit=init_annulus_gap) # Need a wider range to encompass the area beyond the deflector
        self.s_wid = Slider(ax_wid, 'Sky Width (")', 0.1, 3.0, valinit=init_annulus_width)
        self.s_mask = Slider(ax_mask, 'Click Mask Radius (")', 0.05, 1.0, valinit=init_mask_radius)

        self.s_a.on_changed(self.update_plot)
        self.s_ell.on_changed(self.update_plot)
        self.s_theta.on_changed(self.update_plot)
        self.s_gap.on_changed(self.update_plot)
        self.s_wid.on_changed(self.update_plot)

        ax_done = self.fig.add_axes([0.8, 0.01, 0.1, 0.04])
        self.b_done = Button(ax_done, 'Calculate')
        self.b_done.on_clicked(self.finish)

        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.update_plot(None)
        plt.show()

    def get_geometry_pixels(self):
        scale = 1.0 / (self.pix_scale_deg * 3600)
        a = self.s_a.val * scale
        ell = self.s_ell.val
        theta_rad = np.deg2rad(self.s_theta.val)
        b = a * (1.0 - ell)
        a_in = a + (self.s_gap.val * scale)
        a_out = a_in + (self.s_wid.val * scale)
        b_out = a_out * (1.0 - ell)
        return a, b, theta_rad, a_in, a_out, b_out

    def update_plot(self, val):
        [p.remove() for p in reversed(self.ax.patches)]
        a, b, theta, a_in, a_out, b_out = self.get_geometry_pixels()
        
        ap = EllipticalAperture((self.x_t, self.y_t), a=a, b=b, theta=theta)
        ap.plot(ax=self.ax, color='#00FF00', lw=2)
        
        ann = EllipticalAnnulus((self.x_t, self.y_t), a_in=a_in, a_out=a_out, 
                                b_out=b_out, theta=theta)
        ann.plot(ax=self.ax, color='#FFFF00', lw=1.5, ls='--')
        
        if self.mask_data:
            for mx, my, mr in self.mask_data:
                c = plt.Circle((mx, my), mr, color='red', alpha=0.5)
                self.ax.add_patch(c)
        self.fig.canvas.draw_idle()

    def on_click(self, event):
        if event.inaxes != self.ax: return
        if event.button == 1:
            scale = 1.0 / (self.pix_scale_deg * 3600)
            current_r_pix = self.s_mask.val * scale
            self.mask_data.append((event.xdata, event.ydata, current_r_pix))
            self.update_plot(None)
        elif event.button == 3:
            if self.mask_data:
                self.mask_data.pop()
                self.update_plot(None)

    def finish(self, event):
        self.finished = True
        plt.close(self.fig)

# -----------------------------------------------------------------------------
# MAIN SCRIPT
# -----------------------------------------------------------------------------

print(f"Loading {filename}...")
try:
    hdul = fits.open(filename)
except FileNotFoundError:
    print("Error: File not found.")
    sys.exit()

# Handle Cutouts (Ext 0 vs 1)
ext = 1
if len(hdul) > 1 and 'SCI' in hdul[1].name:
    ext = 1
elif len(hdul) == 1:
    ext = 0
header = hdul[ext].header
data = hdul[ext].data
wcs = WCS(header)
data = np.nan_to_num(data, nan=0.0)

pixel_scale_deg = np.sqrt(wcs.proj_plane_pixel_area().value)
target_coords = SkyCoord(ra=target_ra, dec=target_dec, unit='deg', frame='icrs')
x_t, y_t = wcs.world_to_pixel(target_coords)

# =============================================================================
# OPTIONAL: NOISE MAP LOADING
# =============================================================================
aligned_error_map = None

print("\n--- NOISE / ERROR MAP SETUP ---")
ans_noise = input("Do you want to use a Noise/Weight map for error calculation? [y/n]: ").strip().lower()

if ans_noise == 'y':
    
    raw_noise = None
    is_weight_map = False # Default assumption is RMS/Sigma (common for JWST ERR)

    # 1. Check for Internal Extension (e.g., JWST Ext 2)
    use_internal = False
    if len(hdul) > 2:
        print(f"Current file has {len(hdul)} extensions.")
        ans_int = input("Do you want to use Extension [2] of the CURRENT file? (e.g. JWST ERR/WHT) [y/n]: ").strip().lower()
        if ans_int == 'y':
            use_internal = True

    # --- BRANCH A: INTERNAL EXTENSION ---
    if use_internal:
        print("Loading data from Extension [2]...")
        raw_noise = hdul[2].data
        
        # Check header for clues about units (ERR vs WHT)
        try:
            extname = hdul[2].header.get('EXTNAME', '').upper()
            print(f"Extension Name: {extname}")
            if 'WHT' in extname or 'WGT' in extname:
                is_weight_map = True
        except Exception:
            pass
        
        # Ask to confirm if logic didn't auto-detect WHT
        if not is_weight_map:
            ans_type = input(f"Is Ext [2] a WEIGHT map (1/var) [y] or an RMS/SIGMA map [n]? ").strip().lower()
            if ans_type == 'y':
                is_weight_map = True

    # --- BRANCH B: EXTERNAL FILE ---
    else:
        noise_path = input("Enter path to external noise map file: ").strip().replace("'", "").replace('"', "")
        if os.path.exists(noise_path):
            print("Loading external noise map...")
            try:
                with fits.open(noise_path) as hdul_noise:
                    # Heuristic: Find data ext (0 or 1)
                    next_id = 0
                    if len(hdul_noise) > 1 and 'SCI' not in hdul_noise[0].header:
                        next_id = 1
                    
                    ext_data = hdul_noise[next_id].data
                    
                    # Check alignment
                    if ext_data.shape != data.shape:
                        print("External map shape differs from Science. Reprojecting...")
                        if HAS_REPROJECT:
                            # Reproject to match science header
                            raw_noise, _ = reproject_interp(hdul_noise[next_id], header)
                        else:
                            print("Error: 'reproject' needed for alignment. Skipping noise.")
                            raw_noise = None
                    else:
                        raw_noise = ext_data

                    # Guess type based on filename
                    if 'wht' in noise_path.lower() or 'weight' in noise_path.lower():
                        is_weight_map = True

            except Exception as e:
                print(f"Error reading external file: {e}")
                raw_noise = None
        else:
            print("File not found.")

    # --- FINAL PROCESSING (Units & NaNs) ---
    if raw_noise is not None:
        # Sanity check shape (Internal ext could theoretically differ, though unlikely)
        if raw_noise.shape != data.shape:
            print("Warning: Noise map shape mismatch. Reprojecting...")
            if HAS_REPROJECT:
                # If internal, we need a WCS to reproject. Assuming hdul[2] has WCS.
                raw_noise, _ = reproject_interp(hdul[2], header)
            else:
                print("Error: Cannot align map.")
                raw_noise = None

    if raw_noise is not None:
        # Fix NaNs
        raw_noise = np.nan_to_num(raw_noise, nan=0.0)

        if is_weight_map:
            print("Converting Weight map to RMS (Sigma)...")
            # Weight = 1/Sigma^2  --> Sigma = 1/sqrt(Weight)
            mask_good = raw_noise > 0
            aligned_error_map = np.zeros_like(raw_noise)
            aligned_error_map[mask_good] = 1.0 / np.sqrt(raw_noise[mask_good])
        else:
            print("Using as RMS/Sigma map.")
            aligned_error_map = raw_noise
        
        print("Noise map prepared.")

# =============================================================================
# DECISION LOGIC: PARAMS & MASK
# =============================================================================

loaded_params = None
loaded_mask = None

# --- PART 1: APERTURE PARAMETERS ---
if os.path.exists(json_filename):
    print(f"\nFound existing parameter file: {json_filename}")
    ans = input("Use these existing parameters? [y/n]: ").strip().lower()
    if ans == 'y':
        print("Loading local parameters...")
        with open(json_filename, 'r') as f:
            loaded_params = json.load(f)

if loaded_params is None:
    ans = input("\nDo you want to IMPORT parameters from another file? [y/n]: ").strip().lower()
    if ans == 'y':
        import_path = input("Enter path to JSON parameter file: ").strip().replace("'", "").replace('"', "")
        if os.path.exists(import_path):
            print("Loading imported parameters...")
            with open(import_path, 'r') as f:
                loaded_params = json.load(f)
            loaded_params['source'] = 'import'

# --- PART 2: MASK FILE ---
if os.path.exists(mask_filename):
    print(f"\nFound existing mask file: {mask_filename}")
    ans = input("Use this existing mask? [y/n]: ").strip().lower()
    if ans == 'y':
        print("Loading local mask...")
        loaded_mask = fits.getdata(mask_filename).astype(bool)

if loaded_mask is None:
    ans = input("\nDo you want to IMPORT a mask from another file? [y/n]: ").strip().lower()
    if ans == 'y':
        if not HAS_REPROJECT:
            print("Error: 'reproject' package not installed.")
        else:
            mask_import_path = input("Enter path to MASK FITS file: ").strip().replace("'", "").replace('"', "")
            if os.path.exists(mask_import_path):
                print("Importing and reprojecting mask...")
                ext_mask_hdu = fits.open(mask_import_path)[0] 
                reprojected_mask, _ = reproject_interp(ext_mask_hdu, header)
                loaded_mask = (np.nan_to_num(reprojected_mask) > 0.1)
            else:
                print("File not found.")

# =============================================================================
# EXECUTION
# =============================================================================

# CASE A: WE HAVE PARAMETERS
if loaded_params is not None:
    
    geo = loaded_params['geometry']
    scale = 1.0 / (pixel_scale_deg * 3600)
    a_pix = geo['semimajor_axis_arcsec'] * scale
    ell = geo['ellipticity']
    b_pix = a_pix * (1.0 - ell)
    theta_rad = np.deg2rad(geo['theta_deg'])
    
    a_in_pix = (geo['semimajor_axis_arcsec'] + geo['annulus_gap_arcsec']) * scale
    a_out_pix = a_in_pix + (geo['annulus_width_arcsec'] * scale)
    b_out_pix = a_out_pix * (1.0 - ell)
    
    aperture = EllipticalAperture((x_t, y_t), a=a_pix, b=b_pix, theta=theta_rad)
    annulus = EllipticalAnnulus((x_t, y_t), a_in=a_in_pix, a_out=a_out_pix, 
                                b_out=b_out_pix, theta=theta_rad)
    
    if loaded_mask is None:
        final_mask = np.zeros_like(data, dtype=bool)
    else:
        final_mask = loaded_mask

    # 3. VERIFICATION PLOT
    print("\nParameters and Mask prepared.")
    ans_plot = input("Do you want to visualize the aperture and mask on this image before calculating? [y/n]: ").strip().lower()
    
    if ans_plot == 'y':
        plt.figure(figsize=(10, 10))
        mean_val = np.nanmean(data)
        std_val = np.nanstd(data)
        plt.imshow(data, cmap='Greys', origin='lower', 
                   vmin=mean_val - 1*std_val, vmax=mean_val + 5*std_val)
        if np.any(final_mask):
            plt.imshow(final_mask, cmap='Reds', origin='lower', alpha=0.3)
        aperture.plot(color='#00FF00', lw=2, label='Aperture')
        annulus.plot(color='#FFFF00', lw=1.5, ls='--', label='Annulus')
        plt.legend()
        plt.title(f"Visual Check: {filename}")
        plt.show()

    # Save logic
    mask_hdu = fits.PrimaryHDU(data=final_mask.astype('int16'), header=header)
    mask_hdu.writeto(mask_filename, overwrite=True)
    
    loaded_params['pixel_geometry'] = {
        "a": a_pix, "b": b_pix, "theta_rad": theta_rad,
        "annulus_a_in": a_in_pix, "annulus_a_out": a_out_pix
    }
    with open(json_filename, 'w') as f:
        json.dump(loaded_params, f, indent=4)

    run_photometry_math(data, header, aperture, annulus, final_mask, error_map=aligned_error_map)

# CASE B: INTERACTIVE MODE
else:
    print("\nStarting Interactive Tool...")
    tool = PhotometryTool(data, wcs, x_t, y_t, pixel_scale_deg, existing_mask=loaded_mask)

    if not tool.finished:
        print("Window closed without calculation.")
        sys.exit()

    final_a_arcsec = tool.s_a.val
    final_ell = tool.s_ell.val
    final_theta_deg = tool.s_theta.val
    final_gap_arcsec = tool.s_gap.val
    final_width_arcsec = tool.s_wid.val
    
    scale = 1.0 / (pixel_scale_deg * 3600)
    a_pix = final_a_arcsec * scale
    b_pix = a_pix * (1.0 - final_ell)
    theta_rad = np.deg2rad(final_theta_deg)
    
    a_in_pix = (final_a_arcsec + final_gap_arcsec) * scale
    a_out_pix = a_in_pix + (final_width_arcsec * scale)
    b_out_pix = a_out_pix * (1.0 - final_ell)

    final_mask = tool.base_mask.copy()
    if tool.mask_data:
        for (mx, my, mr) in tool.mask_data:
            temp_ap = CircularAperture((mx, my), r=mr)
            ap_mask = temp_ap.to_mask(method='center')
            if ap_mask is not None:
                 img_mask = ap_mask.to_image(data.shape)
                 final_mask = final_mask | (img_mask.astype(bool))

    # Save
    mask_hdu = fits.PrimaryHDU(data=final_mask.astype('int16'), header=header)
    mask_hdu.writeto(mask_filename, overwrite=True)

    params = {
        "target_ra": target_ra,
        "target_dec": target_dec,
        "geometry": {
            "semimajor_axis_arcsec": final_a_arcsec,
            "ellipticity": final_ell,
            "theta_deg": final_theta_deg,
            "annulus_gap_arcsec": final_gap_arcsec,
            "annulus_width_arcsec": final_width_arcsec
        },
        "pixel_geometry": {
            "a": a_pix, "b": b_pix, "theta_rad": theta_rad,
            "annulus_a_in": a_in_pix, "annulus_a_out": a_out_pix
        },
        "masks_pixel_coords": tool.mask_data
    }
    with open(json_filename, 'w') as f:
        json.dump(params, f, indent=4)
        
    print(f"Files saved: {json_filename}, {mask_filename}")

    aperture = EllipticalAperture((x_t, y_t), a=a_pix, b=b_pix, theta=theta_rad)
    annulus = EllipticalAnnulus((x_t, y_t), a_in=a_in_pix, a_out=a_out_pix, 
                                b_out=b_out_pix, theta=theta_rad)
    
    run_photometry_math(data, header, aperture, annulus, final_mask, error_map=aligned_error_map)

hdul.close()