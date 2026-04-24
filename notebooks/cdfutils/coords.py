"""
coords.py

A collection of functions modified from coords.c and fitswcs.c, or new
functions that act as wrappers for functions in MWA's wcs.py
"""

import numpy as np
from math import pi, sin, cos, atan, atan2, sqrt
from astropy import wcs
from astropy import units as u
try:
    from astropy.io import fits as pf
except ImportError:
    import pyfits as pf
try:
    from astropy.coordinates import SkyCoord
except ImportError:
    from astropy.coordinates import ICRS as SkyCoord

# ---------------------------------------------------------------------------


class fileWCS(object):
    """

    Container for WCS information and some associated variables that are
    useful for doing calculations on and selection from astromical imaging
    data that have been stored in a fits file
    """

    def __init__(self, hdr, verbose=True):
        """
        Takes the input fits header (hdr) and extracts useful WCS information
         from it, if that the information is present.  This is essentially
         a simple instatiation of a astropy.wcs.WCS class, but it adds a few
         additional variables that can be derived from the WCS class.

        Inputs:
         hdr     - fits header that may contain WCS information
         verbose - if True (the default) print out possibly useful messages

        """

        """ Initialize some variables """
        self.wcsinfo = None
        self.radec = None
        self.pixscale = None
        self.impa = None

        """ Get the WCS information out of the header"""
        try:
            self.wcsinfo = wcs.WCS(hdr)
        except:
            if verbose:
                print('No WCS information in image header')
            self.wcsinfo = None
            raise KeyError

        """
        Make sure that the WCS information is actually WCS-like and not,
        for example, pixel-based
        """

        imwcs = self.wcsinfo.wcs
        rafound = False
        decfound = False
        count = 0
        for ct in imwcs.ctype:
            if ct[0:2] == 'RA':
                rafound = True
                raax = count
                self.raaxis = raax + 1
                rakey = 'naxis%d' % self.raaxis
            if ct[0:3] == 'DEC':
                decfound = True
                decax = count
                self.decaxis = decax + 1
                deckey = 'naxis%d' % self.decaxis
            count += 1
        if rafound is False or decfound is False:
            if verbose:
                print('No valid WCS information in image header')
                print(' CTYPE keys are not RA/DEC')
            self.wcsinfo = None
            raise KeyError

        """ Get the RA and Dec of the center of the image """
        xcent = hdr[rakey] / 2.
        ycent = hdr[deckey] / 2.
        imcent = np.ones((1, hdr['naxis']))
        imcent[0, raax] = xcent
        imcent[0, decax] = ycent
        imcentradec = self.wcsinfo.wcs_pix2world(imcent, 1)
        self.radec = radec_to_skycoord(imcentradec[0, raax],
                                       imcentradec[0, decax])

        """ Calculate the pixel scale and image PA (E of N) """
        w = self.wcsinfo.wcs
        rad2deg = 180. / pi
        if imwcs.has_cd():
            if verbose:
                print('Using CD matrix to determine pixel scale')
            cdelt1, cdelt2, crot = cdmatrix_to_rscale(w.cd, raax, decax)
            self.pixscale = [3600. * abs(cdelt1), 3600. * abs(cdelt2)]
            self.impa = crot * rad2deg
        elif imwcs.has_pc():
            if verbose:
                print('Using CDELT to determine pixel scale')
            self.pixscale = [3600. * abs(w.cdelt[raax]),
                             3600. * abs(w.cdelt[decax])]
            self.impa = atan(-1. * w.pc[raax, decax] /
                             w.pc[decax, decax]) * rad2deg
        elif isinstance(imwcs.cdelt, np.ndarray):
            self.pixscale = [abs(w.cdelt[raax]) * 3600.,
                             abs(w.cdelt[decax]) * 3600.]
        else:
            print('Warning: no WCS info in header')
            self.wcsinfo = None

        """ Report results if desired """
        if self.wcsinfo is not None and verbose:
            print('')
            print('Rough WCS info')
            print('--------------------------------------------------')
            print('Pixel scale: %7.3f arcsec/pix' % self.pixscale[0])
            print('Instrument FOV (arcsec): %7.1f %7.1f' %
                  (self.pixscale[0] * hdr[rakey],
                   self.pixscale[1] * hdr[deckey]))
            print('Image position angle (E of N): %+7.2f' % self.impa)


# -----------------------------------------------------------------------


def radec_to_skycoord(ra, dec):
    """
    Converts a (RA, dec) pair into the astropy.coordinates SkyCoord
    format

    Required inputs:
     ra  - RA in one of three formats:
            Decimal degrees: ddd.ddddddd  (as many significant figures
              as desired)
            Sexigesimal:     hh mm ss.sss (as many significant figures
              as desired)
            Sexigesimal:     hh:mm:ss.sss (as many significant figures
              as desired)
     dec - Dec in one of three formats:
            Decimal degrees: sddd.ddddddd, where "s" is + or -
            Sexigesimal:     sdd mm ss.sss (as many significant figures
               as desired)
            Sexigesimal:     sdd:mm:ss.sss (as many significant figures
               as desired)
    """

    """ Get RA format """
    if type(ra) == float or type(ra) == np.float32 or \
            type(ra) == np.float64:
        rafmt = u.deg
    else:
        rafmt = u.hourangle

    """ Dec format is always degrees, even if in Sexigesimal format """
    decfmt = u.deg

    """ Do the conversion """
    radec = SkyCoord(ra, dec, unit=(rafmt, decfmt), frame='fk5')
    return radec

# -----------------------------------------------------------------------


def make_header(radec, pixscale, nx, ny=None, rot=None):
    """

    Makes a header with wcs information.

    Required Inputs:
      radec    - The desired (RA, Dec) pair to be put into the CRVAL
                  keywords.
                 This can be in one of two forms:
                  1. A (RA, Dec) pair as either:
                     a. a two-element numpy array
                     b. a two-element tuple
                     c. a two-element list
                  2. a SkyCoord object, as defined in astropy.coordinates
      pixscale - Desired pixel scale in arcsec/pix
      nx       - image size along the x-axis --or-- if the image is square
                   (indicated by ny=None) then this is also the y-axis size

    Optional Inputs:
      ny       - y-axis size, if different from the x-axis size
                   ny=None means that the two axes have the same size
      rot      - desired rotation angle, in DEGREES E of N.

    """

    """
    Convert the input (ra, dec) pair to a SkyCoord object if it is not in
    this format already
    """
    if isinstance(radec, np.ndarray) or isinstance(radec, tuple) or \
            isinstance(radec, list):
        inradec = radec_to_skycoord(radec[0], radec[1])
    elif isinstance(radec, SkyCoord):
        inradec = radec
    else:
        msg = '\nERROR: radec needs to either be a two-element array, list '
        msg += ' or tuple, or\n a SkyCoord object'
        raise ValueError(msg)

    """ Create a blank 2d WCS container """
    w = wcs.WCS(naxis=2)

    """ Get the image size and central pixel """
    cp1 = nx / 2.
    if ny is None:
        cp2 = cp1
    else:
        cp2 = ny / 2.

    """ Make the pixel scale into a two-element quantity """
    tmppix = np.atleast_1d(pixscale)
    if tmppix.size < 2:
        px = np.array([pixscale/3600., pixscale/3600.])
    else:
        px = tmppix / 3600.

    """ Fill the WCS header with appropriate values and save it """
    w.wcs.crpix = [cp1, cp2]
    w.wcs.crval = [inradec.ra.degree, inradec.dec.degree]
    w.wcs.cdelt = [(-1.*px[0]), px[1]]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.equinox = 2000.

    """
    Implement the rotation via the PC matrix, which will have elements:
          cos(theta)    sin(theta)
         -sin(theta)    cos(theta)
     where theta is the angle north through east.
    If you are more used to the CD matrix format, there is a simple expression
     relating the two:

         CD = pixscale * PC

     assuming that the pixel scales are the same along the two axes.
    """
    if rot is not None:
        rotrad = rot * pi
        pc11 = cos(rotrad)
        pc22 = pc11
        pc12 = sin(rotrad)
        pc21 = -pc12
        w.wcs.pc = np.array([[pc11, pc12], [pc21, pc22]])

    """ Convert to a fits header format """
    hdr = w.to_header()
    return hdr, w

# -----------------------------------------------------------------------


def rot_to_pcmatrix(rot, verbose=True):
    """
    Converts an image PA into a PC (rotation) matrix.  Note that the
    PC matrix encodes a rotation differently than a standard rotation matrix.
    For a rotation, theta, the matrix is

       PC1_1     PC1_2          cos(theta)    sin(theta)
                            =
       PC2_1     PC2_2          -sin(theta)   cos(theta)

    Inputs:
      rot       - rotation angle, IN DEGREES (North through East)
    """

    """ Initialize, and convert rotation to radians """
    pcmatx = np.zeros((2, 2))
    if verbose:
        print('')
        print('Creating a new PC matrix with:')
        print('  rotation (degrees N->E):    %+7.2f' % rot)
    rot *= pi / 180.

    """ Set the values """
    pcmatx[0, 0] = cos(rot)
    pcmatx[1, 0] = -1.0 * sin(rot)
    pcmatx[0, 1] = sin(rot)
    pcmatx[1, 1] = cos(rot)

    return pcmatx

# -----------------------------------------------------------------------


def rscale_to_cdmatrix(pixscale, rot, pixscale2=0.0, verbose=True):
    """
    Converts a pixel scale and rotation (in degrees) into a cd matrix,
    which is returned.

    Inputs:
      pixscale  - pixel scale in arcsec/pix (corresponds to CDELT1)
      rot       - rotation angle, IN DEGREES (North through East)
      pixscale2 - [optional] second pixel scale, in arcsec/pix, if
                    CDELT2 != CDELT1.  Default=0.0 means just use pixscale
      rot2      - [optional] second rotation if axes have different rotations.
                    Default=0.0 means that only one rotation is used.
    """

    """ Initialize, and convert rotations to arcsec """
    cdmatx = np.zeros((2, 2))
    if pixscale2 == 0.0:
        pixscale2 = pixscale
    if verbose:
        print('')
        print('Creating a new CD matrix with:')
        print('  pixel scales (arcsec/pix): %7.4f  %7.4f'
              % (pixscale, pixscale2))
        print('  rotation (degrees N->E):    %+7.2f' % rot)

    rot *= pi / 180.

    """ Set the values """
    cdmatx[0, 0] = -1.0 * pixscale * cos(rot)
    cdmatx[1, 0] = -1.0 * pixscale * sin(rot)
    cdmatx[0, 1] = -1.0 * pixscale2 * sin(rot)
    cdmatx[1, 1] = pixscale2 * cos(rot)

    """ Convert cd matrix to degrees and return"""
    cdmatx /= 3600.
    return cdmatx

# -----------------------------------------------------------------------


def matrix_to_rot(matrix, raax=0, decax=1, verbose=True):
    """

    Converts a PC or CD matrix into a rotation

    Inputs:
       matrix  - a CD or PC matrix
       raax    - coordinate axis corresponding to RA
       decax   - coordinate axis corresponding to Dec
       verbose - report progress?  Default is True

    """

    """
    Find the determinant, which sets the sign used in the calculation
    Note that the matrix typically encodes a left-handed coordinate system,
    since RA increases going to the left.
    """

    if np.linalg.det(matrix) < 0:
        cdsgn = -1
    else:
        cdsgn = 1
        if verbose:
            print('  WARNING: matrix_to rot')
            print('    Astrometry is for a right-handed coordinate system.')

    """
    Now convert the PCn_m or CDn_m headers into a rotation.
    Use the following (under the default assumption of axis assignments),
    which is taken from a fits WCS document (fitswcs_draft.pdf) by
    Hanisch and Wells [Look for a published version!]:
      sign(CDELT1*CDELT2) = sign(CD1_1 * CD2_2 - CD1_2 * CD2_1)
         This is determined as cdsgn above
      CROTA2 = arctan(sign * CD1_2 / CD2_2)

    NOTE: CROTA2 is deprecated as a fits keyword, but it is the value of
     the image PA
    """

    crota2 = atan2(cdsgn * matrix[raax, decax], matrix[decax, decax])
    crota2 *= 180. / pi

    return crota2


# -----------------------------------------------------------------------


def cdmatrix_to_rscale(cdmatrix, raax=0, decax=1, verbose=True):
    """
    Converts a CD matrix from a fits header to pixel scales and rotations.

    Inputs:
      cdmatrix - the CD matrix, as a numpy array.  For a standard
                 two-dimensional image, the CD matrix that we care about
                 will be constructed from the fits header keywords as
                  [[CD1_1, CD1_2], [CD2_1, CD2_2]]
                 This would correspond to the default values for the
                  (zero-indexed) raax (default=0) and decax (default=1)
                 However, if the input file has more than two dimensions, then
                  the RA axis is not always the first and the Dec axis is
                  not always the second.  Some 2d images also for some reason
                  have the two axes reversed.  Therefore, use the raax and
                  decax (zero-indexed) are used to indicate which axes
                  correspond to RA and Dec.
      raax     - index for the RA axis.  The default, raax=1, means that CD1_1
                  is associated with RA
      decax    - index for the Dec axis.  The default, decax=1, means that
                  CD2_2 is associated with Dec
      verbose  - set to True (the default) for verbose output
    """

    """
    Find the determinant, which sets the sign used in the calculation
    Note that the CD matrix typically encodes a left-handed coordinate system,
    since RA increases going to the left.
    """

    if np.linalg.det(cdmatrix) < 0:
        cdsgn = -1
    else:
        cdsgn = 1
        if verbose:
            print('  WARNING: cdmatrix_to rscale')
            print('    Astrometry is for a right-handed coordinate system.')

    """
    Now convert CDn_m headers into pixel scales, expresess as CDELTn values.
    Use the following (under the default assumption of axis assignments),
    which is taken from a fits WCS document (fitswcs_draft.pdf) by
    Hanisch and Wells [Look for a published version!]:
      sign(CDELT1*CDELT2) = sign(CD1_1 * CD2_2 - CD1_2 * CD2_1)
         This is determined as cdsgn above
      CDELT1 = sqrt(CD1_1^2 + CD2_1^2)
      CDELT2 = sqrt(CD1_2^2 + CD2_2^2)
      CROTA2 = arctan(sign * CD1_2 / CD2_2)

    NOTE: CROTA2 is deprecated as a fits keyword, but it is the value of
     the image PA
    """

    cdelt1 = cdsgn * sqrt(cdmatrix[raax, raax]**2 + cdmatrix[raax, decax]**2)
    cdelt2 = sqrt(cdmatrix[decax, raax]**2 + cdmatrix[decax, decax]**2)
    crota2 = matrix_to_rot(cdmatrix, raax, decax, verbose)

    """ Return the calculated values """
    return cdelt1, cdelt2, crota2

# -----------------------------------------------------------------------


def update_cdmatrix(fitsfile, pixscale, rot, pixscale2=0.0,
                    hext=0, verbose=True):
    """
    A wrapper for rscale_to_cdmatrix.  In this case, however, the inputs
    include a fits filename, and the cdmatrix generated by the call to
    rscale_to_cdmatrix is written into the header cards of the input fits file.

    Inputs:
      fitsfile  - name of fits file for which the CD matrix will be updated
      pixscale  - pixel scale in arcsec/pix (corresponds to CDELT1)
      rot       - rotation angle, in DEGREES (North through East)
      pixscale2 - [optional] second pixel scale, in arcsec/pix, if
                    CDELT2 != CDELT1.  Default=0.0 means just use pixscale
      rot2      - [optional] second rotation if axes have different rotations.
                    Default=0.0 means that only one rotation is used.
      hext      - [optional] HDU extension if not the default of HDU 0
      verbose   - [optional] flag set to True (the default) for verbose output
    """

    """ Open the input file """
    hdulist = pf.open(fitsfile, mode='update')
    hdr = hdulist[hext].header

    """ Generate a CD matrix based on the inputs """
    cdmatx = rscale_to_cdmatrix(pixscale, rot, pixscale2, verbose)

    """ Transfer new CD matrix to fits header cards and save """
    hdr.update('cd1_1', cdmatx[0, 0])
    hdr.update('cd1_2', cdmatx[0, 1])
    hdr.update('cd2_1', cdmatx[1, 0])
    hdr.update('cd2_2', cdmatx[1, 1])
    hdulist.flush()
    del hdr, hdulist

# -----------------------------------------------------------------------


def sky_to_darcsec(alpha0, delta0, alpha, delta):
    """
    Given a central position on the sky (alpha0,delta0), where these
    coordinates are (for now) expected to be in decimal degrees, and
    one or more other positions (alpha,delta), returns the offsets betwee
    the central position and each of the other positions in arcsec
    on a projected plane.
    Formulas taken from AIPS Memo 27
    For now, this uses only the SIN projection.

    Inputs:
      alpha0 - RA of central position in decimal degrees
      delta0 - Dec of central position in decimal degrees
      alpha  - array containing other RAs in decimal degrees
      delta  - array containing other Decs in decimal degrees
    """

    """ Make sure that inputs are in numpy array format """
    a0 = np.atleast_1d(alpha0).copy()
    d0 = np.atleast_1d(delta0).copy()
    a = np.atleast_1d(alpha).copy()
    d = np.atleast_1d(delta).copy()

    """ Convert inputs into radians """
    degs2rad = np.pi / 180.
    for i in [a0, d0, a, d]:
        i *= degs2rad

    """
    Calculate the direction cosines depending on projection method.
    For now, sin is the only projection.
    """
    proj = 'sin'

    if proj.lower() == 'sin':
        L = np.cos(d) * np.sin(a - a0)
        M = np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0)

    """ Convert direction cosines to arcseconds """
    x_offset = L * 3600. / degs2rad
    y_offset = M * 3600. / degs2rad

    return x_offset, y_offset

# -----------------------------------------------------------------------


def darcsec_to_dpix(dalpha, ddelta, cdmatrix):
    """
    Given a series of (RA, Dec) offsets in arcsec, represented by dalpha and
    ddelta, converts them to pixel offsets for a given fits image.  The
    conversion from arcsec to pixels is done via the CD matrix that has
    been obtained from the fits header.  In particular, if invcd is the
    multiplicative inverse of the cdmatrix, then

      dx = invcd[0,0]*dalpha + invcd[0,1]*ddelta
      dy = invcd[1,0]*dalpha * invcd[1,1]*ddelta

    Inputs:
      dalpha   - vector containing RA offsets in arcsec
      ddelta   - vector containing Dec offsets in arcsec
      cdmatrix - CD matrix obtained from the fits header, where:
                    cdmatx[0,0] = CD1_1
                    cdmatx[0,1] = CD1_2
                    cdmatx[1,0] = CD2_1
                    cdmatx[1,1] = CD2_2
    """

    """ Initialize """
    da = np.atleast_1d(dalpha)
    dd = np.atleast_1d(ddelta)
    if da.size != dd.size:
        print('')
        print('ERROR: darcsec_to_dpix. Input offset vectors must be the same'
              'size')
        return np.nan, np.nan
    dx = np.zeros(da.size)
    dy = np.zeros(dd.size)

    """ Convert CD matrix to arcsec/pix """
    cd = cdmatrix.copy() * 3600.

    """ Invert the CD matrix """
    from numpy import linalg
    invcd = linalg.inv(cd)

    """ Loop through the arcsec offsets and convert to pixel offsets """
    for i in range(da.size):
        dsky = np.array([da[i], dd[i]])
        dpix = np.dot(invcd, dsky)
        dx[i] = dpix[0]
        dy[i] = dpix[1]

    return dx, dy
