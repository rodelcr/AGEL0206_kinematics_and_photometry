"""
A collection of fairly generic code for handling data
"""

import numpy as np
from scipy import interpolate, optimize
from scipy.ndimage import filters
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.modeling import models, fitting

# ---------------------------------------------------------------------------


def sigclip(indata, nsig=3., mask=None, verbose=False):
    """
    Runs a sigma-clipping on the data.  The code iterates over
    the following steps until it has converged:
     1. Compute mean and rms noise in the (clipped) data set
     2. Reject (clip) data points that are more than nsig * rms from
        the newly computed mean (using the newly computed rms)
     3. Repeat until no new points are rejected
        Once convergence has been reached, the final clipped mean and clipped
        rms are returned.

     Optional inputs:
        nsig    - Number of sigma from the mean beyond which points are
                  rejected.  Default=3.
        mask    - If some of the input data are known to be bad, they can
                   be flagged before the inputs are computed by including
                   a mask.  This mask must be set such that True
                   indicates good data and False indicates bad data
        verbose - If False (the default) no information is printed
    """

    """ Determine what the input data set is """
    if mask is not None:
        data = indata[mask]
    else:
        data = indata.copy()

    """ Report the number of valid data points """
    size = data[np.isfinite(data)].size
    if verbose:
        print(' sigma_clip: Full size of data       = %d' % data.size)
        print(' sigma_clip: Number of finite values = %d' % size)

    """
    Reject the non-finite data points and compute the initial values
    """
    d = data[np.isfinite(data)].flatten()
    mu = d.mean()
    sig = d.std()
    mu0 = d.mean()
    sig0 = d.std()
    if verbose:
        print('')
        print('npix = %11d. mean = %f. sigma = %f' % (size, mu, sig))

    """ Iterate until convergence """
    delta = 1
    clipError = False
    while delta:
        size = d.size
        if sig == 0.:
            clipError = True
            break
        d = d[abs(d - mu) < nsig * sig]
        mu = d.mean()
        sig = d.std()
        if verbose:
            print('npix = %11d. mean = %f. sigma = %f' % (size, mu, sig))
        delta = size-d.size

    """ Error handling """
    if clipError:
        print('')
        print('ERROR: sigma clipping failed, perhaps with rms=0.')
        print('Setting mu and sig to their original, unclipped, values')
        print('')
        mu = mu0
        sig = sig0

    """ Clean up and return the results and clean up """
    del data, d
    return mu, sig

# ===========================================================================
#
# Start of DataGen class
#
# ===========================================================================


class DataGen(object):
    """

    The most generic data-handling methods.  Right now this is empty since
    the sigclip method, which had been in this class, has been moved out
    to be stand-alone since it was not working well as part of this class.

    """

    # -----------------------------------------------------------------------

    def __init__(self, data):
        self.data = np.asarray(data)

    # -----------------------------------------------------------------------

# ===========================================================================
#
# Start of Data1d class
#
# ===========================================================================


class Data1d(Table):
    """

    Code to perform actions on 1-dimesional data sets such as light curves or
    1d spectra.

    """

    # -----------------------------------------------------------------------

    def __init__(self, x, y, var=None, names=None, debug=False):
        """

        Reads in the data, which can be expressed as y = y(x).

        """

        """ Set up the default names """
        if names is None:
            if var is None:
                names = ['x', 'y']
            else:
                names = ['x', 'y', 'var']
        if(debug):
            print(names)
            print('Length of x vector: %d' % x.size)
            print('Length of y vector: %d' % y.size)

        """ Link to the inherited class """
        if var is None:
            Table.__init__(self, [x, y], names=names)
        else:
            Table.__init__(self, [x, y, var], names=names)

        """ Assign simple names to the columns """
        self.x = self[self.colnames[0]]
        self.y = self[self.colnames[1]]
        if var is not None:
            self.var = self[self.colnames[2]]
        else:
            self.var = None

    # -----------------------------------------------------------------------

    def resamp(self, xout=None, verbose=True):
        """
        Resample the data vector (y) onto a new spacing in x.
        There are two possibilities for the output x vector that sets where
         the interpolation happens.  They are:
         1. xout = None [default]
            A linearized set of spacings between the minimum and maximum
            values in the input x vector
         2. xout is set to an array
            A user-defined x array that has been passed through the xout
            parameter
        """

        if xout is None:
            x0 = self.x[0]
            x1 = self.x.max()
            xout = np.linspace(x0, x1, self.x.size)

        ymod = interpolate.splrep(self.x, self.y)
        yout = interpolate.splev(xout, ymod)

        """ Return the resampled vectors """
        if verbose:
            print('resample: replacing input spectrum with resampled'
                  ' version')
            print('resample: for now not resampling the variance')
        return xout, yout

    # -----------------------------------------------------------------------

    def smooth_boxcar(self, filtwidth, verbose=True):
        """
        Does a boxcar smooth of the spectrum.
        The default is to do inverse variance weighting, using the variance
         spectrum if it exists.
        The other default is not to write out an output file.  This can be
        changed by setting the outfile parameter.
        """

        """ Set the weighting """
        if self.var is not None:
            if verbose:
                print('Weighting by the inverse variance')
            wht = 1.0 / self.var
        else:
            if verbose:
                print('Uniform weighting')
            wht = 0.0 * self.y + 1.0

        """ Smooth the spectrum and variance spectrum """
        yin = wht * self.y
        smowht = filters.uniform_filter(wht, filtwidth)
        ysmooth = filters.uniform_filter(yin, filtwidth)
        ysmooth /= smowht
        if self.var is not None:
            varsmooth = 1.0 / (filtwidth * smowht)
        else:
            varsmooth = None

        return ysmooth, varsmooth

    # -----------------------------------------------------------------------

    def fit_poly(self, fitorder, mod0=None, y0=None, fitrange=None, nsig=3.0,
                 doplot=True, markformat='bo', xlabel='x', ylabel='y',
                 title=None):
        """

        This method is essentially just a pair of calls to numpy's polyfit
        and polyval.  However, there are a few add-ons, including the clipping
        of outliers (set by nsig), an optional limitation on the x range of
        the data to be fit (set by fitrange), and optional plotting of the
        resulting fit.

        """

        """ First a sigma clipping to reject clear outliers """
        if fitrange is None:
            tmpfitdat = self.y.data.copy()
        else:
            fitmask = np.logical_and(self.x >= fitrange[0],
                                     self.x < fitrange[1])
            tmpfitdat = self.y.data[fitmask]
        dmu, dsig = sigclip(tmpfitdat, nsig=nsig)
        goodmask = np.absolute(self.y.data - dmu) < nsig * dsig
        badmask = np.absolute(self.y.data - dmu) >= nsig * dsig
        xgood = self.x.data[goodmask]
        dgood = self.y.data[goodmask]
        xbad = self.x.data[badmask]
        dbad = self.y.data[badmask]

        """ Limit the range of fitted points, if requested """
        if fitrange is None:
            xpoly = xgood
            dpoly = dgood
        else:
            fitmask = np.logical_and(xgood >= fitrange[0], xgood < fitrange[1])
            xpoly = xgood[fitmask]
            dpoly = dgood[fitmask]

        """ Fit a polynomial to the trace """
        if fitorder > -1:
            dpoly = np.polyfit(xpoly, dpoly, fitorder)
        elif y0 is not None:
            dpoly = np.array([y0, ])
        else:
            dpoly = None

        """
        Calculate the fitted function
        """
        fitx = np.arange(self.x.size)
        if dpoly is not None:
            fity = np.polyval(dpoly, fitx)
        else:
            fity is None

        """ Plot the results """
        ymin = dmu - 4.5*dsig
        ymax = dmu + 4.5*dsig
        if doplot and fity is not None:
            plt.plot(self.x, self.y.data, markformat)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)

            """
            Show an input constant value (e.g., the mean y value)
            This input value, if it is passed to the method, would have been
             generated before calling the function
            """
            if y0 is not None:
                plt.axhline(y0, color='k', linestyle='--')

            """ Mark the bad points that were not included in the fit """
            plt.plot(xbad, dbad, "rx", markersize=10, markeredgewidth=2)

            """ Show the fitted function """
            plt.plot(fitx, fity, "r")

            """
            Show the range of points included in the fit, if fitrange was set
            """
            if fitrange is not None:
                plt.axvline(fitrange[0], color='k', linestyle=':')
                plt.axvline(fitrange[1], color='k', linestyle=':')
                xtmp = 0.5 * (fitrange[1] + fitrange[0])
                xerr = xtmp - fitrange[0]
                ytmp = fity.min() - 0.2 * fity.min()
                plt.errorbar(xtmp, ytmp, xerr=xerr, ecolor="g", capsize=10)
            plt.xlim(0, self.x.max())
            plt.ylim(ymin, ymax)

        """
        Return the parameters produced by the fit and the fitted function
        """
        print(dpoly)
        return dpoly, fity

    # -----------------------------------------------------------------------

    def fit_mod(self, mod0, usevar=False, fitrange=None, verbose=True):
        """

        Fits a pre-defined model (given by mod0) to the data.  The model
         must have been set by some other method before calling this function.

        Required inputs:
         mod0     - The input model, in the form of an astropy.modeling object

        Optional inputs:
         usevar   - Use the variance array to define the uncertainties to be
                     used in the fitting.  If set to False (the default) then
                     uniform weighting is used.
        fitrange - x-axis range to use for the fitting.  If fitrange is set
                    to None (the default), then the entire data set is used
                    to constrain the model
        """

        """ Summarize the initial guess """
        if verbose:
            print('')
            print('Initial model')
            print('-------------')
            print(mod0)
            print('')
            print('-------------------------------------------')
            print('')

        """ 
        If variance-weighted fitting is requested, make sure that there
        is a variance vector
        """
        if usevar:
            if self.var is None:
                raise KeyError('*** ERROR:  fit_mod\n'
                               'Fitting requested variance weighting but'
                               ' no variance array found.')
            else:
                rms0 = np.sqrt(self.var)
        else:
            rms0 = np.ones(self.x.size)
               
                
        """ Set the range over which to fit the data """
        if fitrange is not None:
            mask = (self.x > fitrange[0]) & (self.x < fitrange[1])
            x = self.x[mask]
            y = self.y[mask]
            rms = rms0[mask]
        else:
            x = self.x.copy()
            y = self.y.copy()
            rms = rms0

        """
        Do the fitting.
        NOTE: The 'weights' for the fitter, if requested, are 1/RMS and NOT
          1/var because that is what the astropy fitters expect.  
          This must be because their figure of merit is set to
            (y_data - y_mod)*weight
          which is later squared somewhere, giving a real figure of merit
            of  (y_data - y_mod)**2 / sigma**2   since weight = 1/sigma
        """
        fit = fitting.LevMarLSQFitter()
        
        """ 
        look into fitinfo, which has an 'additional info' thing which is
        a dictionary that has as one of its keys param_cov (for covariance)
        NB: this only works (for now) for LevMar
        """
        mod = fit(mod0, x, y, weights=1.0/rms)

        """
        Clean up and return best-fit model
        Also return the 'fitinfo' dictionary.  This dictionary contains
         (among other things):
          'param_cov' - covariance matrix for the parameters
          'ier'       - integer indicating fit quality: 1, 2, 3, 4 are
                        successful fits
        """
        del rms, x, y
        return mod, fit.fit_info

    # -----------------------------------------------------------------------

    def fit_gauss(self, bgorder=0, smo=5, gtype='em', usevar=False,
                  mod0=None, bounds=None, fitrange=None, verbose=True):
        """
        Fits a Gaussian plus a background to the data.  The background
         is represented by a polynomial of degree bgorder.  The default value,
         bgorder=0, gives a constant background.
        The data are modeled using the astropy modeling package.  The 
         parameters that are used for the two components are:
           * Background polynomial: c0, c1, ...  [up to bgorder]
           * Gaussian: amplitude, mean, stddev
        """

        """
        If there is no input model, then we need to set up a model and
        give it some input guesses
        """
        if mod0 is not None:
            m_init = mod0
            
        else:
            """
            Do a temporary smoothing of the data to reduce the effect of
            noise on the initial guesses
            """
            tmpsmooth, junk = self.smooth_boxcar(smo, verbose=False)

            """
            The default model for the background is just a constant
            Do a sigma clipping to estimate the base level from the data
            (may or may not get used later)
            """
            base, tmp = sigclip(tmpsmooth)

            """ Set up the background polynomial """
            p = models.Polynomial1D(degree=bgorder, c0=base)

            """
            Get the initial guesses for the Gaussian.
            """
            if gtype == 'abs':
                amp0 = tmpsmooth.min() - base
                mu0 = self.x[np.argmin(tmpsmooth)]
            else:
                amp0 = tmpsmooth.max() - base
                mu0 = self.x[np.argmax(tmpsmooth)]
            sig0 = 3.5

            del tmpsmooth

            """
            Create the initial-guess model
            NOTE: Should probably add bounds
            """
            g = models.Gaussian1D(amplitude=amp0, mean=mu0, stddev=sig0)
            m_init = p + g

        mod, fit_info = self.fit_mod(m_init, usevar=usevar, fitrange=fitrange,
                                     verbose=verbose)
        return mod, fit_info

    # -----------------------------------------------------------------------

    def _make_gauss(self, p):
        """

        Creates a model comprised of one or more Gaussian profiles plus a
        (for now) constant background.
        NOTE: the only oddity is that, if there are more than one Gaussian in
        the profile, then the "mu" term for the subsequent Gaussians
        (i.e., p[4], p[7], ..) are actually _offsets_ between the mean of the
        subsequent Gaussian and the mean of the first.  For example,
        mu_2 = p[0] + p[4]

        Inputs:
          p  - The parameter values.  The length of this vector will be 1+3*n,
               where n>=1, for one constant background parameter (p[0]) plus
               one or more Gaussian parameters, which come in sets of three.
               Thus, p can be decomposed as follows:
                 p[0] - background: required
                 p[1] - mu_1: required
                 p[2] - sigma_1: required
                 p[3] - amplitude_1: required
                 p[4] - offset between mu_2 and mu_1: optional
                 p[5] - sigma_2: optional
                 p[6] - amplitude_2: optional
                 ... etc. for as many Gaussians are used to construct the
                   profile
        """

        """ Calculate the number of Gaussians in the model """
        ngauss = int((p.size-1)/3)
        if p.size - (ngauss*3+1) != 0:
            print('')
            print('ERROR: Gaussian model contains the incorrect number of'
                  'parameters')
            print('')
            return np.nan

        """ Calculate y_mod using current parameter values """
        ymod = np.zeros(self.x.size) + p[0]
        for i in range(ngauss):
            ind = i*3+1
            if i == 0:
                mu = p[ind]
            else:
                mu = p[1] + p[ind]

            ymod += p[ind+2] * np.exp(-0.5 * ((self.x - mu)/p[ind+1])**2)

        return ymod

    # -----------------------------------------------------------------------

    def _checkmod_gauss(self, p, p_init, fitind):
        """

        Compares the data to the model.  The model consists of at least one
        gaussian plus a constant background and is created by a call to
        make_gauss.
        Thus the comparison is between ymod(x) and y, where the latter is the
        measured quantity.

        NOTE: the only oddity in the model is that, if there are more than
        one Gaussian in the profile, then the "mu" term for the subsequent
        Gaussians (i.e., p[4], p[7], ..) are actually _offsets_ between the
        mean of the subsequent Gaussian and the mean of the first.  For
        example, mu_2 = p[0] + p[4]

        Inputs:
          p  - The parameter values.  The length of this vector will be 1+3*n,
               where n>=1, for one constant background parameter (p[0]) plus
               one or more Gaussian parameters, which come in sets of three.
               Thus, p can be decomposed as follows:
                 p[0] - background: required
                 p[1] - mu_1: required
                 p[2] - sigma_1: required
                 p[3] - amplitude_1: required
                 p[4] - offset between mu_2 and mu_1: optional
                 p[5] - sigma_2: optional
                 p[6] - amplitude_2: optional
                 ... etc. for as many Gaussians are used to construct the
                 profile
        """

        """
        Create the full list of model parameters by combining the fitted
        parameters and the fixed parameters
        """
        pfull = p_init.copy()
        pfull[fitind] = p

        """
        Compute the difference between model and real values
        """
        ymod = self._make_gauss(pfull)
        diff = self.y - ymod

        return diff

    # -----------------------------------------------------------------------

    def fit_gauss_old(self, init=None, fix=None, ngauss=1, verbose=True):
        """

        Old routine for fitting a gaussian plus background.  This approach
        uses the scipy.optimize.leastsq routine rather than the astropy
        modeling routines.

        """

        """
        Set up the container for the initial guesses for the parameter values
        """
        nparam = 3*ngauss + 1
        p_init = np.zeros(nparam)

        """
        Put default choices, which may be overridden, into p_init
        """
        p_init[0] = np.median(self.y, axis=None)
        for i in range(ngauss):
            """
            In this loop the parameters are set as follows:
              p_init[ind] is either mu (if i==0) or an offset from p_init[1]
              p_init[ind+1] is sigma
              p_init[ind+2] is amplitude
            """
            ind = 3*i + 1
            if i == 0:
                tmp = self.y.argsort()
                p_init[ind] = 1.0 * tmp[tmp.shape[0]-1]
            else:
                p_init[ind] = 5. * i * (-1.)**(i+1)
            p_init[ind+1] = 3.
            p_init[ind+2] = self.y.max() - p_init[0]

        """
        Override the default values if init has been set.
        NOTE: the init parameter must be an array (or list) of length nparam
         otherwise the method will quit
        NOTE: A value of -999 in init means keep the default value
        """
        if init is not None:
            if len(init) != nparam:
                print('')
                print('ERROR: locate_trace -- init parameter must have'
                      'length %d' % nparam)
                print('  (since ngauss=%d ==> nparam= 3*%d +1)' %
                      (ngauss, ngauss))
                print('')
                return np.nan
            for j in range(nparam):
                if init[j] > -998.:
                    p_init[j] = init[j]

        """
        Set up which parameters are fixed based on the fix parameter
        The default value (fix=None) means that all of the parameters are
         varied in the fitting process
        """
        fixstr = np.zeros(nparam, dtype='S3')
        if fix is None:
            fixvec = np.zeros(nparam)
        else:
            fixvec = np.atleast_1d(fix)
            if fixvec.size != nparam:
                print('')
                print('ERROR: locate_trace - fix parameter must have length %d'
                      % nparam)
                print('  (since ngauss=%d ==> nparam= 3*%d +1)' % 
                      (ngauss, ngauss))
                print('')
                return np.nan
            fixstr[fixvec == 1] = 'Yes'
        fitmask = fixvec == 0
        fitind = np.arange(nparam)[fitmask]

        """ Fit a Gaussian plus a background to the compressed spectrum """
        mf = 100000
        p = p_init[fitmask]
        p_fit, ier = optimize.leastsq(self._checkmod_gauss, p,
                                      (p_init, fitind),
                                      maxfev=mf)
        """
        Create the full parameter list for the fit by combining the fitted
        parameters and the fixed parameters
        """
        p_out = p_init.copy()
        p_out[fitind] = p_fit

        """ Give results """
        if(verbose):
            print('')
            print('Profile fit results')
            print('-------------------')
            print('                                        Held')
            print('Parameter          Init Value  fixed? Final Value')
            print('--------------    ----------  ------ -----------')
            print('background         %9.3f     %3s     %9.3f'
                  % (p_init[0], fixstr[0], p_out[0]))
            for i in range(ngauss):
                ind = 3 * i + 1
                j = i + 1
                if i == 0:
                    mustr = 'mu_1'
                else:
                    mustr = 'offset_%d' % j
                print('%-9s          %9.3f     %3s     %9.3f'
                      % (mustr, p_init[ind], fixstr[ind], p_out[ind]))
                print('sigma_%d            %9.3f    %3s    %9.3f'
                      % (j, p_init[ind+1], fixstr[ind+1], p_out[ind+1]))
                print('amp_%d              %9.3f    %3s    %9.3f'
                      % (j, p_init[ind+2], fixstr[ind+2], p_out[ind+2]))
            print('')

        return p_out
