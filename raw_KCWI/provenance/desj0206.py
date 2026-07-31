from ISMgas.kcwi.kcwiReduction import kcwiRedux
import astropy.units as u

filenames = [
    "redux/kb250830_00052_icubes.fits",
    "redux/kb250830_00053_icubes.fits",
    "redux/kb250830_00054_icubes.fits",
    "../../KCWI_arctomo/KCWI_Nov16/redux/kb251117_00087_icubes.fits",
    "../../KCWI_arctomo/KCWI_Nov16/redux/kb251117_00089_icubes.fits",
    "../../KCWI_arctomo/KCWI_Nov16/redux/kb251117_00089_icubes.fits",
    "../../KCWI_arctomo/KCWI_Nov16/redux/kb251117_00090_icubes.fits",
    "../../KCWI_arctomo/KCWI_Nov16/redux/kb251117_00091_icubes.fits",    
]

objid = "DESJ0206"
ra  = 31.5561
dec = -1.2382
resolution = 0.3*u.arcsecond
size= 100
slicer='medium'

obj  = kcwiRedux(
    objid, 
    ra, dec, 
    resolution, size, 
    filenames, 
    slicer, 
    autocorrelate=True, 
    grab=False,
)

obj.step1(refCombine=0)
obj.step2()