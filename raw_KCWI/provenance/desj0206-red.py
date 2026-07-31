from ISMgas.kcwi.kcwiReduction import kcwiRedux
import astropy.units as u

filenames = [
    "desj0206/desj0206-stack-1.fits",
    "desj0206/desj0206-stack-2.fits",
    "desj0206/desj0206-stack-3.fits",

]

objid = "DESJ0206-red"
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
    grab=True 
)

obj.step1(refCombine=1)
obj.step2()