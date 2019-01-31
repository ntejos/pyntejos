import numpy as np
from astropy import units as u
"""Module to handle DS9 regions utilities"""


def get_boxregion_string(x, y, a, b, angle, text=None):
    """Returns a ds9 string associated to a box with position (x,y) and angle
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
    text : str
        A comment for the region

    """
    s = 'box({:.8f},{:.8f},{},{},{:.8f}) # text={{{}}}\n'.format(x, y, a, b, angle, text)
    return s


def create_regfile_nboxes_from_slit(output, xy, PA, l, w, n, pixscale, add_label=True):
    """
    Creates a ds9 region file with n boxes within a rectangular
    region representing a slit of width w and length l in arcsecs, centred in pixel xy and
    with a position angle PA. [Note: n must be odd integer]

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
        Odd positive integer
    pixscale : float
        Pixel scale in arcsec
    add_label : bool
        Whether to add a label to the region, (1 to n) where 1 is the northern one when PA=0.

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
    angle =  PA  # angle from x to y coordinates is the same than PA from North to East in the usual convention
    angle_rad = angle * np.pi / 180.
    for ii in range(-1*nside, nside+1):
        xnew = xy[0] + ii * b * np.sin(angle_rad)
        ynew = xy[1] - ii * b * np.cos(angle_rad)
        text = None
        if add_label:
            text = str(ii+nside+1)
        s = get_boxregion_string(xnew, ynew, a, b, angle, text=text)
        f.write(s)
    f.close()
