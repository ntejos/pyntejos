import numpy as np
from astropy import units as u
"""Module to handle DS9 regions utilities"""


def get_boxregion_string(x, y, a, b, angle):
    """Returns a ds9 string associated to a box with position (x,y) and PA
    Parameters
    ----------
    x : float
        Position of the box in x-coordinates
    y : float
        Position of the box in y-coordinates
    a : float
        Size of the box in the x-coordinates at PA=0 (pixel units)
    b : float
        Size of the box in the y-coordinates at PA=0 (pixel units)
    angle : float
        angle of the box for ds9

    """
    s = 'box({:.8f},{:.8f},{},{},{:.8f}) ||\n'.format(x, y, a, b, angle)
    return s


def create_regfile_nboxes_from_slit(output, xy, PA, l, w, n, pixscale):
    """

    Based on F. Corro script

    Parameters
    ----------
    output : str
        Name of the output file
    xy : (float, float)
        Central position of the slit in (x,y) pixel coordinates
        Where x-coordinates grows to west and y-coordinate grows to north
    PA : float
        Position angle in degrees of the slit (from north to east). Assuming convention that
        x-coordinates grows to west and y-coordinate grows to north
    l : float
        Length of the slit in arcsecs
    w : float
        Width of the slit in arsecs
    n : int
        Odd integer
    pixscale : float
        Pixel scale in arcsec

    Returns
    -------
    A ds9 region file writen in `output`, that has n aligned boxes within the slit

    """

    if (n%2 == 0) or (n < 0):
        raise ValueError("n must be odd positive integer!")

    f = open(output, "w+")
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('physical\n')

    nside = (n - 1)/2
    a = w / pixscale
    b = l / n / pixscale
    angle =  PA  # angle in x-y coordinates is offset from PA
    angle_rad = angle * np.pi / 180.
    for ii in range(-1*nside, nside+1):
        xnew = xy[0] - ii * b * np.sin(angle_rad)
        ynew = xy[1] + ii * b * np.cos(angle_rad)
        s = get_boxregion_string(xnew, ynew, a, b, angle)
        f.write(s)
    f.close()

def create_regfile_nboxes_(output, xy, PA, n, w=1*u.arcsec):
    """ Creates a region file with n boxes within a mage slit


    Parameters
    ----------
    output
    xy
    PA
    n
    w

    Returns
    -------

    """
    pixscale = 0.3*u.arcsec
    l = 10*u.arcsec
    lp = l.to('arcsec').value / pixscale.to('arcsec').value
    wp = w.to('arcsec').value / pixscale.to('arcsec').value
    create_regfile_nboxes_from_slit(output, xy, PA, lp, wp, n)


