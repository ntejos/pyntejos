""" Module with various useful functions that don't really belong in any
of the other modules. From Neil Crighton's old astro"""

import cPickle as pickle
import os, pdb
import math
import sys
from pprint import pformat
from textwrap import wrap
if float(sys.version[:3]) < 2.4:
    from sets import ImmutableSet as set
import numpy as np

Ckms = 299792.458         # speed of light km/s, exact
pi = math.pi

class Bunch(object):
    """Bunch class from the python cookbook with __str__ and __repr__
    methods. Similar to an IDL structure.

    >>> s = Bunch(a=1, b=2, c=['bar', 99])
    >>> s.a
    1
    >>> s.c
    ['bar', 99]
    >>> s
    Bunch(a, b, c)
    >>> print s
    Bunch(
    a = 1
    b = 2
    c = ['bar', 99])
    """
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
    def __repr__(self):
        temp = ', '.join(sorted(str(attr) for attr in self.__dict__ if not str(attr).startswith('_')))
        return 'Bunch(%s)' % '\n      '.join(wrap(temp, width=69))

def nan2num(a, replace=0):
    """ Replace nan or inf entries with the `replace` keyword value.

    If replace is 'mean', use the mean of the array to replace
    values. If it's 'interp', then intepolate from the nearest values.
    """
    b = np.array(a, copy=True)
    bad = np.isnan(b) | np.isinf(b)
    if replace == 'mean':
        replace = b[~bad].mean().astype(b.dtype)
    elif replace == 'interp':
        x = np.arange(len(a))
        replace = np.interp(x[bad], x[~bad], b[~bad]).astype(b.dtype)
        
    b[bad] = replace
    return b

def indgroupby(a, name):
    """ Find the indices giving rows in `a` that have common values
    for the field `name`.

    Inputs
    ------
    a: structured array
    name: field of `a` (str)

    Returns
    -------
    unique_vals, indices:
        unique sorted values of a[name], and a list of indices into
        `a` giving rows with each unique value.
    """
    b = a[name]
    isort = b.argsort()
    # find indices in the sorted array where we change to a new value
    # of name.
    temp = np.nonzero(b[isort][:-1] != b[isort][1:])[0]
    jbreaks = np.concatenate([[0], temp + 1, [len(b)]])
    # make a list of indices into a for each unique value
    for j0,j1 in zip(jbreaks[:-1], jbreaks[1:]):
        ind = isort[j0:j1].view(np.ndarray)
        # yield the unique value, and the indices into `a` giving rows
        # with that value
        yield b[ind[0]], ind

def between(a, vmin, vmax):
    """ Return a boolean array True where vmin <= a < vmax.

    (Careful of floating point issues when dealing with equalities
    though)
    """
    a = np.asarray(a)
    c = a < vmax
    c &= a >= vmin
    return c

# algorithm from scipy
def scoreatpercentile(a, perc):
    """Calculate the score at the given 'perc' percentile of the
    sequence a.  For example, the score at perc=50 is the median.

    'perc' can be a sequence of percentile values.

    If the desired quantile lies between two data points, we linearly
    interpolate between them.
    """
    # TODO: this should be a simple wrapper around a well-written quantile
    # function.  GNU R provides 9 quantile algorithms (!), with differing
    # behaviour at, for example, discontinuities.
    vals = np.sort(a, axis=0)

    if not hasattr(perc, '__iter__'):
        perc = [perc]

    out = []
    for p in perc:
        i = 0.01 * p * (vals.shape[0] - 1)
        j = int(i)
        if (i % 1 == 0):
            out.append(vals[j])
        else:
            # linearly interpolate
            out.append(vals[j] + (vals[j+1] - vals[j]) * (i % 1))

    return np.array(out).squeeze()

    
def poisson_noise(flux, sigma, seed=None):
    """ Adds poisson noise to a normalised flux array.

    sigma: One sigma error in the flux at the continuum level (where
    normalised flux=1).

    flux: Array of normalised flux values (i.e. flux values divided by
    the continuum).

    If `seed` is given, it is used to seed the random number
    generator.

    Returns:  flux with noise added, one sigma error array.

    Tests
    -----
    >>> fl = np.linspace(0,1)
    >>> fl0,er0 = poisson_noise(fl, 0.1, seed=114)
    >>> fl1,er1 = np.loadtxt('data/noisepoisson.txt.gz', unpack=1)
    >>> print np.allclose(fl0,fl1), np.allclose(er0,er1)
    True True
    """

    if seed is not None:  np.random.seed(seed)

    flux = np.asarray(flux)
    sigma = float(sigma)
    if np.any(flux < 0):  raise Exception('Flux values must be >= 0!')
    lamb = flux / (sigma * sigma)              # variance per pixel
    flnew = np.random.poisson(lamb)
    flnew = np.where(lamb > 0, flnew/lamb * flux, 0)
    sig = np.where(lamb > 0, flux/np.sqrt(lamb), 0)
    return flnew,sig

def addnoise(flux, sigma, minsig=None, poisson=True, seed=None):
    """ Add noise to a normalised flux array.

    Either gaussian, poisson, or a combination of both noise types is
    added to the flux array, depending on the keyword minsig.

    Parameters
    ----------
    flux: array_like
      Array of normalised flux values.
    sigma: float
      Total desired noise at the continuum (flux=1). Note the
      SNR = 1 / sigma.
    minsig: float, optional
      By default minsig is `None`, which means gaussian noise with
      standard deviation `sigma` is added to the flux. If minsig is
      given, a combination of poisson and gaussian noise is added to
      the flux to give an error of `sigma` at the continuum. In this
      case the gaussian noise component has st. dev. of `minsig`,
      which must be less than `sigma`.

    seed: int, optional
      If seed is given, it is used to seed the random number generator.
      By default the seed is not reset.

    Returns
    -------
    flux with noise added, one sigma error array.

    Tests
    -----
    >>> fl = np.linspace(0,1)
    >>> fl0,er0 = addnoise(fl, 0.2, seed=113)
    >>> fl1,er1 = np.loadtxt('data/noisegauss.txt.gz', unpack=1)
    >>> print np.allclose(fl0,fl1), np.allclose(er0,er1)
    True True
    >>> fl0,er0 = addnoise(fl, 0.2, minsig=0.05, seed=116)
    >>> fl1,er1 = np.loadtxt('data/noiseboth.txt.gz', unpack=1)
    >>> print np.allclose(fl0,fl1), np.allclose(er0,er1)
    True True
    """
    sigma = abs(float(sigma))

    if minsig is None:
        flux = np.asarray(flux)
        if seed is not None:  np.random.seed(seed)
        dev = sigma * np.random.randn(len(flux))
        er = np.empty_like(flux)
        er.fill(sigma)
        return flux + dev, er
    else:
        minsig = abs(float(minsig))
        if minsig > sigma:
            raise Exception('Noise at continuum must be bigger than minsig!')
        # gaussian variance
        var_g = minsig*minsig
        # normalised sigma of poisson noise at the continuum
        sig_p_cont = np.sqrt(sigma*sigma - var_g)
        if poisson:
            flnew,sig_p = poisson_noise(flux, sig_p_cont, seed=seed)
            # gaussian error
            flnew += minsig * np.random.randn(len(flux))
            # total sigma
            er = np.sqrt(sig_p * sig_p + var_g)
        else:
            er = minsig + flux * sig_p_cont
            flnew = np.asarray(flux).copy()
            flnew += er * np.random.randn(len(flux))
        
        return flnew, er


def wmean(val, sig):
    """ Return the weighted mean and error.

    Uses inverse variances as weights.

    val: array with shape (N,)
      Array of values and

    sig: array with shape (N,)
      One sigma errors (sqrt(variance)) of the array values.

    Returns
    -------
    wmean, wsigma: floats
      The weighted mean and error on the weighted mean.

    Tests
    -----
    >>> val = np.concatenate([np.ones(100),np.ones(100)*2])
    >>> sig = range(1,201)
    >>> mean,err = wmean(val,sig)
    >>> np.allclose((mean, err), (1.003026102, 0.78088153553))
    True
    """
    val = np.asarray(val)
    sig = np.asarray(sig)

    # remove any values with bad errors
    condition = (sig > 0.) & ~np.isnan(val) & ~np.isnan(sig)
    if not condition.any():
        raise ValueError('No good values!')
    val = val[condition]
    sig = sig[condition]

    # normalisation
    inverse_variance = 1. / (sig*sig)
    norm = np.sum(inverse_variance)

    wmean = np.sum(inverse_variance*val) / norm
    wsigma = 1. / np.sqrt(norm)

    return wmean,wsigma


def indexnear(ar, val):
    """ Find the element in an array closest to a given value.

    The input array must be sorted lowest to highest.  Returns the
    index of the element with a value closest to the given value.

    Parameters
    ----------
    ar: array_like
      Input array. It must be sorted smallest to largest.
    val: float
      Find the element of `ar` that is closest to `val`.

    Returns
    -------
    index: int
      Index of the `ar` element with the closest value to `val`.

    Examples
    --------
    >>> wa = np.arange(4000, 4500, 0.051)
    >>> i = indexnear(wa, 4302.5)
    >>> print i, wa[i]
    5931 4302.481
    >>> i = indexnear(wa, 4600.0)
    >>> print i, wa[i]
    9803 4499.953
    >>> i = indexnear(wa, 3000.0)
    >>> print i, wa[i]
    0 4000.0
    """
    # TODO: change to be ufunc-like, using np.where?
    ar = np.asarray(ar)
    i = ar.searchsorted(val)
    # needed because searchsort rounds up
    if i == 0:
        return i
    # note if i == len(ar) then ar[i] is invalid, but won't get tested.
    elif i == len(ar) or val - ar[i-1] < ar[i] - val:
        return i-1
    else:
        return i


def calc_Mstar_b(z):
    """ Find the Schechter parameter M* in the rest frame b band at
    redshift z, by interpolating over the Faber at al. 2007 DEEP2
    averaged values, and assuming M*_b = -20.0 at z=0 (rough average
    of the z=0.07-0.1 points in Faber 2007) .
    """
    zvals = 0.0, 0.3, 0.5, 0.7, 0.9, 1.1
    Mvals = -20.00, -21.07, -21.15, -21.51, -21.36, -21.54
    return np.interp(z, zvals, Mvals)


def combinations(items, n):
    """ A generator for the number of ways you can take n items (order
    unimportant) from a list of items."""
    if n == 0:
        yield []
    else:
        for i in xrange(len(items)):
            for c in combinations(items[i+1:], n-1):
                yield [items[i]] + c

def permutations(items):
    """ Permutations are just a special case of combinations."""
    return combinations(items, len(items))

def find_edges_true_regions(condition):
    """ Finds the indices for the edges of contiguous regions where
    condition is True.

    Examples
    --------
    >>> a = np.array([3,0,1,4,6,7,8,6,3,2,0,3,4,5,6,4,2,0,2,5,0,3])
    >>> ileft, iright = find_edges_true_regions(a > 2)
    >>> zip(ileft, iright)
    [(0, 0), (3, 8), (11, 15), (19, 19), (21, 21)]

    """
    indices, = condition.nonzero()
    if not len(indices):
        return None, None
    iright, = (indices[1:] - indices[:-1] > 1).nonzero()
    ileft = iright + 1
    iright = np.concatenate( (iright, [len(indices)-1]) )
    ileft = np.concatenate( ([0], ileft) )
    return indices[ileft], indices[iright]

def stats(arr):
    """ Show the minimum, maximum and mean of an array

    Also show the number of NaN entries (if any).
    """
    
    arr = np.asarray(arr)
    shape = arr.shape
    arr = arr.ravel()
    size = len(arr)
    bad = np.isnan(arr)
    nbad = bad.sum()
    if nbad == size:
        return '#NaN %i of %i' % (nbad, size)
    elif nbad == 0:
        arr = np.sort(arr)
    else:
        arr = np.sort(arr[~bad])
    if len(arr) % 2 == 0:
        i = len(arr) // 2
        median = 0.5 * (arr[i-1] + arr[i])
    else:
        median = arr[len(arr) // 2]

    return 'min %.3g max %.3g median %.3g mean %.3g shape %s #NaN %i of %i' % (
        arr[0], arr[-1], median, arr.mean(), shape, nbad, size)

def Gaussian(x, x0, sigma, height):
    """ Gaussian."""
    return height * np.exp(-0.5 * ((x-x0)/sigma)**2)

