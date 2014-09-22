"""Functions for fitting distributions. Code mostly taken from Neil
Crighton's old astro package

"""

from __future__ import division

import numpy as np
import inspect, os

from scipy.special import gammainc
from scipy.optimize import leastsq
import matplotlib.pyplot as pl
import matplotlib.transforms as mtran

import lm
from utilities import \
     scoreatpercentile, indexnear, Bunch, between,  Gaussian
from convolve import convolve_psf
from spec import plotlines
from plot import puttext
from io import readtxt
from spline import InterpCubicSpline


Ckms = 299792.458                 # speed of light km/s, exact
wlya = 1215.6701

# read mean continuum
def read_pc_all():
    """ Principle components (eigenvectors of the covariance matrix)
    for Suzuki et al. (2005) QSO spectra from 1020 to 1600 Angstroms.
    """
    path =  os.path.abspath(os.path.dirname(__file__))
    filename = path + '/PCAcont/Suzuki05/tab3.txt'
    names = 'wa,mu,musig,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10'
    p = readtxt(filename, skip=23, names=names)
    pc = np.array([p['e%i' % i] for i in range(1, 11)])
    return p.wa, p.mu, pc


def findchisq(func, x, y, ysig, a):
    """ Calculate chi2 for function and data values

    func: model function, func(x, *a) models y data
    x,y: data points
    ysig: 1 sigma errors on y
    a: model parameters
    """
    resid = (y - func(x, *a)) / ysig
    return np.dot(resid, resid)


def polyfitr(x, y, order=2, clip=6, xlim=None, ylim=None,
             mask=None, debug=False):
    """ Fit a polynomial to data, rejecting outliers.

    Fits a polynomial f(x) to data, x,y.  Finds standard deviation of
    y - f(x) and removes points that differ from f(x) by more than
    clip*stddev, then refits.  This repeats until no points are
    removed.

    Inputs
    ------
    x,y:
        Data points to be fitted.  They must have the same length.
    order: int (2)
        Order of polynomial to be fitted.
    clip: float (6)
        After each iteration data further than this many standard
        deviations away from the fit will be discarded.
    xlim: tuple of maximum and minimum x values, optional
        Data outside these x limits will not be used in the fit.
    ylim: tuple of maximum and minimum y values, optional
        As for xlim, but for y data.
    mask: sequence of pairs, optional
        A list of minimum and maximum x values (e.g. [(3, 4), (8, 9)])
        giving regions to be excluded from the fit.
    debug: boolean, default False
        If True, plots the fit at each iteration in matplotlib.

    Returns
    -------
    coeff, x, y:
        x, y are the data points contributing to the final fit. coeff
        gives the coefficients of the final polynomial fit (use
        np.polyval(coeff,x)).

    Examples
    --------
    >>> x = np.linspace(0,4)
    >>> np.random.seed(13)
    >>> y = x**2 + np.random.randn(50)
    >>> coeff, x1, y1 = polyfitr(x, y)
    >>> np.allclose(coeff, [1.05228393, -0.31855442, 0.4957111])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, order=1, xlim=(0.5,3.5), ylim=(1,10))
    >>> np.allclose(coeff, [3.23959627, -1.81635911])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, mask=[(1, 2), (3, 3.5)])
    >>> np.allclose(coeff, [1.08044631, -0.37032771, 0.42847982])
    True
    """

    good = ~np.isnan(x) & ~np.isnan(y)
    x = np.asanyarray(x[good])
    y = np.asanyarray(y[good])
    isort = x.argsort()
    x, y = x[isort], y[isort]

    keep = np.ones(len(x), bool)
    if xlim is not None:
        keep &= (xlim[0] < x) & (x < xlim[1])
    if ylim is not None:
        keep &= (ylim[0] < y) & (y < ylim[1])
    if mask is not None:
        badpts = np.zeros(len(x), bool)
        for x0,x1 in mask:
            badpts |=  (x0 < x) & (x < x1)
        keep &= ~badpts

    x,y = x[keep], y[keep]
    if debug:
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y,'.')
        ax.set_autoscale_on(0)
        pl.show()

    coeff = np.polyfit(x, y, order)
    if debug:
        pts, = ax.plot(x, y, '.')
        poly, = ax.plot(x, np.polyval(coeff, x), lw=2)
        pl.show()
        raw_input('Enter to continue')
    norm = np.abs(y - np.polyval(coeff, x))
    stdev = np.std(norm)
    condition =  norm < clip * stdev
    y = y[condition]
    x = x[condition]
    while norm.max() > clip * stdev:
        if len(y) < order + 1:
            raise Exception('Too few points left to fit!')
        coeff = np.polyfit(x, y, order)
        if debug:
            pts.set_data(x, y)
            poly.set_data(x, np.polyval(coeff, x))
            pl.show()
            raw_input('Enter to continue')
        norm = np.abs(y - np.polyval(coeff, x))
        stdev = norm.std()
        condition =  norm < clip * stdev
        y = y[condition]
        x = x[condition]

    return coeff,x,y

def wleastsq(x, y, ysig=None):
    """ Calculate the line of best fit with weights.

    Input
    -----
      x : sequence of floats
          input x data
      y : sequence of floats, len(x)
          input y data, f(x).
      ysig : sequence of floats, len(x), optional
         If the y none sigma errors are given, points are weighted
         by their inverse variance.

    Returns
    -------
      (b, a),(b_sigma, a_sigma)
        The fitted parameters and their one sigma errors.  The fitted
        line has equation y = a + b*x.
    """
    if ysig == None:
        ysig = np.ones(len(x))
    yvar = ysig * ysig   # variance

    s = np.sum(1. / yvar)
    sx = np.sum(x / yvar)
    sxx = np.sum(x*x / yvar)
    sy = np.sum(y / yvar)
    sxy = np.sum(x*y / yvar)

    # See NR section 15.2 for a derivation of the below solutions for
    # the best fit values of a and b.
    # 
    # y = a + b*x 

    temp = s*sxx - sx*sx
    a = (sxx*sy - sx*sxy) / temp
    b = (s*sxy - sx*sy) / temp
    sig_a = np.sqrt(sxx / temp)
    sig_b = np.sqrt(s / temp)

    return (b,a),(sig_b,sig_a)

def fitgauss(x, y, sig, sigma, x0, height):
    """ Fits a gaussian to data with errors using chi**2 minimisation.

    Inputs
    ------
    x, y:
        data values.
    sig:
        One sigma errors in y data values.
    sigma, x0, height:
        Initial guess values for the gaussian parameters.

    Returns
    -------
    results:
        Best fitting sigma, x position and height, and parameter
        covariance matrix (diagonals give the variance for each best
        fitting parameter).

    Example
    -------
    >>> from fit import fitgauss
    >>> def g(x, sigma, x0, height): 
    ...     return height * np.exp(-0.5 * ((x-x0) / sigma) ** 2)
    >>> x = np.linspace(0,50, 100)
    >>> ptrue = 5.8, 22.2, 9
    >>> plot(x, g(x, *ptrue), lw=2)
    >>> y = g(x, *ptrue) + np.random.randn(len(x))*5
    >>> plot(x, y)
    >>> pbest,cov = fitgauss(x,y,[5]*len(y), 10, 20, 3)
    >>> plot(x, g(x,*pbest))
    >>> perr = np.sqrt(cov.diagonal())
    >>> plot(x, g(x,pbest[0]+perr[0], *pbest[1:]), 'r')
    >>> plot(x, g(x,pbest[0]-perr[0], *pbest[1:]), 'r')
    >>> plot(x, g(x,pbest[0], pbest[1]-perr[1], pbest[2]), 'm')
    >>> plot(x, g(x,pbest[0], pbest[1]+perr[1], pbest[2]), 'm')
    >>> plot(x, g(x,pbest[0], pbest[1], pbest[2]-perr[2]), 'c')
    >>> plot(x, g(x,pbest[0], pbest[1], pbest[2]+perr[2]), 'c')
    """
    x,y,sig = map(np.asarray, (x,y,sig))
    def g(x, sigma, x0, height):
        """ Gaussian."""
        return height * np.exp(-0.5 * ((x-x0) / sigma) ** 2)

    guesses = sigma, x0, height
    popt, r = nhmc_fit(g, x, y, sig, guesses)
    print 'Reduced chi2 = %s' % (findchisq(g, x, y, sig, popt) / (len(x) - 3))
    return popt, r


def nhmc_fit(func, xdata, ydata, ysig, guess, **kw):
    """
    Use non-linear least squares to fit a function to data.

    Assumes ``ydata = f(xdata, *params) + eps``

    Parameters
    ----------
    func : callable
        The model function, func(x, ...).  It must take the independent
        variable as the first argument and the M parameters to fit as
        separate remaining arguments.
    xdata : An N-length sequence
        The independent variable where the data is measured.
    ydata : N-length sequence
        The dependent data - nominally func(xdata, ...).
    ysig : N-length sequence
        The standard-deviation of ydata.  These will be used as weights
        in the least-squares problem.
    guess : M-length sequence
        Initial guess for the parameters.

    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the
        squared values of ``(func(xdata, *popt) - ydata)/ysig`` is
        minimized.
    r : object
        Lots of information about the model and fit, including
        parameter covariance matrix, correaltion matrix, paramter one
        sigma errors.

    Notes
    -----
    The algorithm uses the Levenburg-Marquardt algorithm:
    scipy.optimize.leastsq. Additional keyword arguments are passed directly
    to that algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from fit import nhmc_fit
    >>> def func(x, a, b, c):
    ...     return a*np.exp(-b*x) + c

    >>> x = np.linspace(0,4,50)
    >>> y = func(x, 2.5, 1.3, 0.5)
    >>> ydata = y + 0.2*np.random.rand(len(x))
    >>> plot(x, y)
    >>> plot(x, ydata)
    >>> popt, r = nhmc_fit(func, x, ydata, [0.2]*len(x), (1,1,1))
    >>> plot(x, func(x, *popt))

    """
    xdata, ydata, ysig = (np.asarray(a) for a in (xdata, ydata, ysig))
    args = (xdata, ydata, func, 1./ysig)
    npar = len(guess)

    def f(par, xdata, ydata, function, weights):
        return weights * (function(xdata, *par) - ydata)

    def findchisq(par):
        resid = (ydata - func(xdata, *par)) / ysig
        return np.dot(resid, resid)

    popt, pcov, infodict, mesg, ier = leastsq(f, guess, args=args, 
                                              full_output=1, **kw)
    
    #if ier != 1:
    #    raise RuntimeError, "Optimal parameters not found: " + mesg

    # not sure why this is necessary...
    if (len(ydata) > npar) and pcov is not None:
        assert np.allclose((f(popt, *args)**2).sum(), findchisq(popt))
        s_sq = (f(popt, *args)**2).sum()/(len(ydata) - npar)
        pcov = pcov * s_sq
        parsig = np.sqrt(pcov.diagonal())
        corr = np.empty((npar, npar), float)
        for i in xrange(npar):
            for j in xrange(npar):
                corr[i,j] = pcov[i,j] / (parsig[i] * parsig[j])
    else:
        pcov = None
        parsig = corr = None

    #fcode = func.__code__
    #parnames = inspect.getargs(fcode).args[1:]

    results = Bunch(par=popt, parsig=parsig,
                    cov=pcov, corr=corr,
                    func=func, x=xdata, y=ydata, ysig=ysig,
                    guess=guess, findchisq=findchisq,
                    normchisq=findchisq(popt)/(len(xdata) - npar),
                    info=dict(infodict=infodict, mesg=mesg, ier=ier),
                    #parnames=parnames,
                    )

    results.__doc__ = """\
parnames  : parameter names from func call signature
par       : best-fitting parameters
parsig    : 1 sigma errors in best-fitting parameter values
cov       : covariance matrix of parameter values
corr      : normalised covariance matric (-1 -> 1)
func      : the input function: func(x, *par) is the best fitting model
x,y       : x, y data
ysig      : 1 sigma y errors
guess     : input parameter guesses
findchisq : function to find chisq from parameters. e.g. findchisq(par)
normchisq : chisq per degree of freedom for best fitting parameters
info      : information from the leastsq fit output
"""
    return popt, results

def _general_function(params, xdata, ydata, function):
    return function(xdata, *params) - ydata

def _weighted_general_function(params, xdata, ydata, function, weights):
    return weights * (function(xdata, *params) - ydata)

def curve_fit(f, xdata, ydata, p0=None, sigma=None, fulloutput=True, **kw):
    """
    Use non-linear least squares to fit a function, f, to data.

    Assumes ``ydata = f(xdata, *params) + eps``

    Parameters
    ----------
    f : callable
        The model function, f(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : An N-length sequence or an (k,N)-shaped array
            for functions with k predictors.
        The independent variable where the data is measured.
    ydata : N-length sequence
        The dependent data --- nominally f(xdata, ...)
    p0 : None, scalar, or M-length sequence
        Initial guess for the parameters.  If None, then the initial
        values will all be 1 (if the number of parameters for the function
        can be determined using introspection, otherwise a ValueError
        is raised).
    sigma : None or N-length sequence
        If not None, it represents the standard-deviation of ydata.
        This vector, if given, will be used as weights in the
        least-squares problem.

    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared error
        of ``f(xdata, *popt) - ydata`` is minimized
    pcov : 2d array
        The estimated covariance of popt.  The diagonals provide the variance
        of the parameter estimates.

    Notes
    -----
    The algorithm uses the Levenburg-Marquardt algorithm:
    scipy.optimize.leastsq. Additional keyword arguments are passed directly
    to that algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import curve_fit
    >>> def func(x, a, b, c):
    ...     return a*np.exp(-b*x) + c

    >>> x = np.linspace(0,4,50)
    >>> y = func(x, 2.5, 1.3, 0.5)
    >>> yn = y + 0.2*np.random.normal(size=len(x))

    >>> popt, pcov = curve_fit(func, x, yn)

    """
    if p0 is None or np.isscalar(p0):
        # determine number of parameters by inspecting the function
        import inspect
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if len(args) < 2:
            raise ValueError, "p0 not given as a sequence and inspection"\
                " cannot determine the number of fit parameters"
        if p0 is None:
            p0 = 1.0
        p0 = [p0]*(len(args)-1)

    args = (xdata, ydata, f)
    if sigma is None:
        func = _general_function
    else:
        func = _weighted_general_function
        args += (1.0/np.asarray(sigma),)
    popt, pcov, infodict, mesg, ier = leastsq(func, p0, args=args, 
                                              full_output=1, **kw)
    
    #if ier != 1:
    #    raise RuntimeError, "Optimal parameters not found: " + mesg

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    if fulloutput:
        return popt, pcov, infodict, mesg, ier
    else:
        return popt, pcov


def median_continuum(flux, error, numsig=1.5, plot=False):
    """ Estimates the continuum using a median and sigma clipping.

    Given the fluxes and one sigma errors for the section, calculates
    the flux median. Then rejects all flux values less than numsig*sig
    lower than the median, where sig is the median one sigma error
    value of the flux values. Repeat this process, only retaining the
    not-rejected flux values each time, until no flux values are
    rejected. Once this condition is reached, take the current median
    value as the continuum.

    Returns the continuum value.
    """

    if plot:
        pl.plot(flux)
    while True:
        medfl = np.median(flux)
        meder = np.median(error)
        if plot:
            l = pl.axhline(medfl)
        cond = (flux > (medfl - meder*numsig))
        badflux = flux[~cond]
        if len(badflux) == 0:
            return medfl
        flux = flux[cond]
        error = error[cond]


def lyafcont(wa, fl):
    """ Using the model of Bernardi et al, guess the continuum for the
    lya forest. Only works for flux-calibrated spectra. Continuum is
    only valid over the lya forest region! The wa array must extend
    over the lya forest and the reat wa region 1330-1360 Ang.

    Inputs
    ------
    wa, fl: arrays, shape (N,)

    Returns
    -------
    continuum: array, shape (N,)
    """
    c0 = 0.0224
    c1 = -1.56
    c2 = 0.0033
    c3 = 1073
    c4 = 9
    c5 = 0.0023
    c6 = 1123
    c7 = 9
    c8 = 0.021
    c9 = 1216
    c10 = 29
    # only use closest region for now, to minimise flux calibration
    # problems
    powlaw_regions = (1330, 1360), #(1440, 1490), (1650, 1680)

    # work out c0 from (hopefully) clean regions outside forest
    guesses = []
    for w1,w2 in powlaw_regions:
        cond = (w1 <= wa) & (wa <= w2)
        if not np.any(cond):
            print 'no good regions found!'
            continue
        print 'mean wa', np.mean(wa[cond])
        print 'median fl', np.median(fl[cond])
        c0guess = np.median(fl[cond]) / (c0 * (np.mean(wa[cond]) / wlya) ** c1)
        guesses.append(c0guess)

    if len(guesses) > 0:
        mult = np.median(guesses)
    else:
        mult = 1

    print 'multiplier=',mult
    return mult * ( c0 * (wa / wlya) ** c1 +
                    c2 * np.exp(-0.5 * ((wa - c3) / c4)**2) +
                    c5 * np.exp(-0.5 * ((wa - c6) / c7)**2) +
                    c8 * np.exp(-0.5 * ((wa - c9) / c10)**2) )
    
def poly_scale_spec(s, sref, mask=None, order=9, clip=None, debug=True):
    """ Remove the large-scale variations between two spectra by
    dividing them by one another and then fitting a low-order
    polynomial.

    Note we assume the reference spectrum covers a larger wavelength
    range than the spectrum to be scaled.

    s           : Spectrum object to be scaled to match reference
    sref        : Reference spectrum object
    order = 9   : Order of polynomial to fit
    mask = None : Input to polyfitr, masks wavelength regions.
    clip = None : Number of sigma to clip after each iteration.

    Returns
    -------
    out: array
        sr * out matches s
    """

    # Assume reference spectrum covers a larger wavelength range.
    s1 = sref.rebin(wstart=s.wa[0], dw=s.dw, dv=s.dv, npts=len(s.wa))

    c0 = (s1.er > 0) & (s.er > 0)
    c1 = ~np.isnan(s1.fl) & ~np.isnan(s.fl) & (s1.fl!=0) & (s.fl!=0)
    good = c0 & c1
    scale = s1.fl[good] / s.fl[good]
    wa = s.wa[good]
    pl.plot(wa, scale, '.')
    if clip is not None:
        coeff,x,y = polyfitr(wa, scale, order=order, mask=mask,
                             debug=debug, clip=clip)
    else:
        coeff,x,y = polyfitr(wa, scale, order=order, mask=mask, debug=debug)

    pl.plot(s.wa, np.polyval(coeff,s.wa))

    return np.polyval(coeff,s.wa)

def poly_scale_spec2(wa, fl, waref, flref, mask=None, order=9,
                    clip=None, debug=False):
    """ Remove the large-scale variations between two spectra by
    dividing them by one another and then fitting a low-order
    polynomial.

    Note we assume the reference spectrum covers a larger wavelength
    range than the spectrum to be scaled.

    wa, fl      : Spectrum to be scaled to match reference
    waref, flref: Reference spectrum
    order = 9   : Order of polynomial to fit
    mask = None : Input to polyfitr, masks wavelength regions.
    clip = None : Number of sigma to clip after each iteration.

    Returns
    -------
    out: array
        fl * out matches flref where they overlap.
    """

    #find overlap
    wmin = max(wa[0], waref[0])
    wmax = min(wa[-1], waref[-1])
    cond = (wmin < waref) & (waref < wmax)
    flref1 = flref[cond]
    waref1 = waref[cond]
    fl1 = np.interp(waref1, wa, fl)

    good = ~np.isnan(flref1) & ~np.isnan(fl1) & (flref1 != 0) & (fl1 != 0)
    scale = flref1[good] / fl1[good]
    wa2 = waref1[good]
    if debug:
        pl.plot(wa2, scale, '.')
    if clip is not None:
        coeff,x,y = polyfitr(wa2, scale, order=order, mask=mask,
                             debug=debug, clip=clip)
    else:
        coeff,x,y = polyfitr(wa2, scale, order=order, mask=mask, debug=debug)

    if debug:
        pl.plot(s.wa, np.polyval(coeff, wa))

    return np.polyval(coeff, wa)

def scale_by_median(wa, fl, waref, flref, mask=None):
    """ Find multiplier necessary to make median of two fluxes match."""

    #assert len(wa) == len(fl)
    #assert len(waref) == len(flref)
    wmin = max(wa[0], waref[0])
    wmax = min(wa[-1], waref[-1])
    cond = (wmin < waref) & (waref < wmax)
    flref1 = flref[cond]
    waref1 = waref[cond]
    fl1 = np.interp(waref1, wa, fl)

    masked = np.zeros(len(waref1), dtype=bool)
    for wmin, wmax in mask:
        i,j = waref1.searchsorted([wmin,wmax])
        masked[i:j] = True

    c0 = ~np.isnan(flref1) & ~np.isnan(fl1) & (flref1 != 0) & (fl1 != 0)\
         & ~masked
    assert ((len(flref1[c0]) > 0) & (len(fl1[c0]) > 0)), "No overlap!"

    mult = np.median(flref1[c0]) / np.median(fl1[c0])

    return mult

def delta_chi2(ptarg, ndeg):
    """ Find delta chi2 corresponding to the given probability for n
    degrees of freedom (number of data points - number of free
    parameters).

    ptarg is the target probability (ptarg=0.9 correponds to 90%).

    ptarg is the probability that a model fitted to data with the
    given number of degrees of freedom will have a chisq val
    delta_chi2 larger than the minimum chisq value for the
    best-fitting parameters.
    """
    assert 0 <= ptarg <= 1
    d0, d1 = 0, 1e5
    a = 0.5 * ndeg
    while True:
        d = 0.5*(d0 + d1)
        p = gammainc(a, 0.5*d)
        #print ptarg, p, d, d0, d1
        if p > ptarg:
            d1 = d
        else:
            d0 = d
        if abs(p - ptarg) < 1e-4:
            break

    return d


