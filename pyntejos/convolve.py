"""From Neil Crighton's old astro package"""

import numpy as np
from utilities import nan2num

def convolve_psf(a, fwhm, edge='invert', replace_nan=True, debug=False):
    """ Given an array of values 'a' and a gaussian full width at half
    maximum 'fwhm' in pixel units, returns the convolution of the
    array with the normalised gaussian.

    Gaussian is calculated for as many pixels required until it drops
    to 1% of peak value. Note that the data will be spoiled at
    distances n/2 (rounding down) from the edges, where n is width of
    the gaussian in pixels.

    The FWHM given should be > 2 pixels to sample the gaussian PSF
    properly.
    """
    const2   = 2.354820046             # 2*sqrt(2*ln(2))
    const100 = 3.034854259             # sqrt(2*ln(100))
    sigma = fwhm / const2
    # gaussian drops to 1/100 of maximum value at x =
    # sqrt(2*ln(100))*sigma, so number of pixels to include from
    # centre of gaussian is:
    n = np.ceil(const100 * sigma)
    if replace_nan:
        a = nan2num(a, replace='interp')
    if debug:
        print "First and last %s pixels of output array will be invalid" % n
    x = np.linspace(-n, n, 2*n + 1)        # total no. of pixels = 2n+1
    gauss = np.exp(-0.5 * (x / sigma) ** 2 )

    return convolve_window(a, gauss, edge=edge)

def convolve_window(a, window, edge='invert'):
    """ Convolve a with a window array. The window array should have
    an odd number of elements.

    The edges of `a` are extended by reflection (and, if requested,
    inversion) to reduce edge effects.

    edge = {'extend', 'reflect', 'invert'} or an integer

    The window is normalised.
    """
    npts = len(window)
    if not npts % 2:
        raise ValueError('`window` must have an odd number of elements!')

    n = npts // 2

    # normalise the window
    window /= window.sum()

    # Add edges to either end of the array to reduce edge effects in
    # the convolution.
    if len(a) < 2*n:
        raise ValueError('Window is too big for the array! %i %i' % (len(a), n))
    if edge == 'invert':
        temp1 = 2*a[0] - a[n:0:-1], a, 2*a[-1] - a[-2:-n-2:-1]
    elif edge == 'reflect':
        temp1 =  a[n:0:-1], a, a[-2:-n-2:-1]
    elif edge == 'extend':
        temp1 =  a[0] * np.ones(n) , a, a[-1] * np.ones(n)
    else:
        try:
            int(edge)
        except TypeError:
            raise ValueError('Unknown value for edge keyword: %s' % edge)
        med1 = np.median(a[:edge])
        med2 = np.median(a[-edge:])
        temp1 =  med1 * np.ones(n) , a, med2 * np.ones(n)

    temp2 = np.convolve(np.concatenate(temp1), window, mode='same')

    return temp2[n:-n]
