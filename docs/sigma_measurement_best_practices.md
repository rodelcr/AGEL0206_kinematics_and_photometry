# Velocity Dispersion Measurement: Best Practices for M-sigma

## Literature Review Summary (compiled for AGEL0206 ApJL paper)

### The standard quantity: sigma_e (effective velocity dispersion)

**Definition:** The luminosity-weighted second velocity moment within the effective
radius R_e. For IFU data, this is typically measured by co-adding all spectra within
R_e and fitting the integrated spectrum with ppxf.

**Key references:**
- Cappellari et al. (2006, 2013): Defined sigma_e as the luminosity-weighted
  line-of-sight velocity dispersion within 1 R_e, using eq. 29 of Cappellari+2013.
- Emsellem et al. (2007): Flux-weighted method for sigma_e from IFU data.

### Two common conventions for M-sigma

**1. sigma_e (within R_e):**
Used by: Cappellari+2013, van den Bosch 2016, Sahu+2024
- sigma_e^2 = <V_rms^2> = <V^2 + sigma^2>_lum-weighted within R_e
- Accounts for both rotation and random motions
- Standard for dynamical modeling (JAM, Schwarzschild)

**2. sigma_c (central, within R_e/8):**
Used by: Ferrarese & Merritt 2000, Kormendy & Ho 2013, Greene+2020
- Corrected to a standard aperture R_e/8 using Jorgensen+1995 formula
- Power-law correction: sigma_c/sigma_aper = (R_aper/R_e/8)^alpha
- alpha ≈ -0.04 to -0.06 for early-type galaxies (Jorgensen+1995)
- Recent work (Sohn+2024, MaNGA data) shows alpha varies with M_r, color, Sersic n

### Aperture correction formula (Jorgensen, Franx & Kjaergaard 1995)

    sigma_corrected / sigma_measured = (R_aperture / R_standard)^alpha

    alpha = -0.04 (JFK95 original, for early-type galaxies)
    R_standard = R_e/8 (for sigma_c) or R_e (for sigma_e)

For our case: if we measure sigma within R_aperture=1" and R_e=2.6",
then R_aperture/R_e = 0.38, and:
- sigma_e = sigma(<R_e) directly (no correction needed if we integrate within R_e)
- sigma_c = sigma(<1") * (1"/0.325")^(-0.04) = sigma * 1.048 ≈ sigma * 1.05

### What Melo-Carneiro+2025 did (Cosmic Horseshoe, our closest analog)

1. **Measured R_e = 2.10"** from MGE surface brightness model
2. **Voronoi binned** IFU data (target SNR=15, min pixel SNR=2.5)
3. **Masked source emission** regions during binning
4. **Fitted ppxf** on each Voronoi bin (5600-7600 Å rest-frame)
5. **Computed V_rms = sqrt(V^2 + sigma^2)** per bin for radial profile
6. **For M-sigma:** co-added ALL spectra within R_e, fitted integrated spectrum
   → sigma_e directly, no aperture correction
7. **Found V_rms profile:** nearly flat in center, rises at edges
8. **Validated:** flux-weighted method (Emsellem+2007) gives same sigma_e
9. **Result:** sigma_e = 366 ± 6 km/s → M_BH = 3.6 × 10^10 M_sun

### Aperture correction from MaNGA IFU data (Sohn et al. 2024, arXiv:2307.12251)

- 10,000 MaNGA galaxies with spatially-resolved kinematics
- sigma_aper/sigma_e = (R_aper/R_e)^alpha
- alpha varies with: M_r (absolute magnitude), g-i color, Sersic n
- For massive, red, high-n galaxies: alpha ≈ -0.03 to -0.06
- For spirals/late-types: alpha can be very different (even positive)
- **Key warning:** "aperture correction derived from early-type galaxies
  cannot be applied to galaxies with intermediate Sersic indices"

### Recommendations for AGEL0206

**Our situation:**
- Spiral lens galaxy at z=0.675
- KCWI IFU with 0.3"/spaxel, seeing FWHM=1.27"
- Lensed arc contamination at R > 1"
- R_e ≈ 2.6" (from white-light curve of growth)
- sigma(R<1") ≈ 210 km/s (robust, method-independent)
- sigma(R<R_e) ≈ 330-430 km/s (contamination-dependent)

**Option A: sigma within R_e (Melo-Carneiro approach)**
- Co-add spectra within R_e with contamination weighting
- Report sigma_e directly
- Problem: heavy arc contamination at R>1" makes this unreliable

**Option B: Central sigma with aperture correction**
- Use sigma(R<1") ≈ 210 km/s as sigma_aper
- Apply JFK95 correction to R_e/8: sigma_c = 210 * (1.0/0.33)^(-0.04) ≈ 206 km/s
- Problem: correction assumes early-type profile; this is a spiral

**Option C: sigma within cleanest aperture (R<R_e/2)**
- sigma(R<R_e/2) ≈ 223-238 km/s depending on masking
- This is well inside the contamination zone
- Report as sigma(R<R_e/2) with explicit aperture specification
- Apply Sohn+2024 correction for the galaxy's Sersic n

**Option D: V_rms approach**
- Our Test 17 shows V=75 km/s was purely a z offset (true z ≈ 0.67555)
- At correct z: V≈0, so V_rms ≈ sigma ≈ 210 km/s
- Melo-Carneiro+2025 found V_rms ≈ sigma for their system too

**Recommended approach:**
1. Report sigma(R<1") = 210 ± 15 km/s as the primary measurement
   (contamination-free, z-independent, method-robust)
2. Note the aperture: R=1" = 7 kpc = 0.38 R_e
3. Apply minimal aperture correction to R_e/8 if needed: ≈ 206 km/s
4. Flag the spiral morphology as a caveat for JFK95 correction
5. Report the radial profile showing the rising sigma for context
6. Distinguish sigma from V_rms (they're equal for this galaxy at correct z)

### Key references

- Cappellari & Emsellem (2004): pPXF method
- Cappellari (2017): Improved pPXF with Gauss-Hermite functions
- Cappellari (2023): pPXF with photometry, MNRAS 526, 3273 (arXiv:2208.14974)
- Emsellem et al. (2007): Flux-weighted sigma_e from SAURON IFU
- Ferrarese & Merritt (2000): M-sigma, sigma within R_e/8
- Gebhardt et al. (2000): M-sigma, sigma within R_e/8
- Greene et al. (2020): Updated M-sigma relations for different galaxy types
- Jorgensen, Franx & Kjaergaard (1995): Aperture correction formula
- Kormendy & Ho (2013, ARAA 51): Comprehensive M-sigma review, sigma_c convention
- Melo-Carneiro et al. (2025, MNRAS 541, 2853): Cosmic Horseshoe, sigma_e from IFU
- Sohn et al. (2024, arXiv:2307.12251): Aperture correction from 10,000 MaNGA galaxies
- van den Bosch (2016): M-sigma with fundamental plane, sigma_e
