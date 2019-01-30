import numpy as np
"""Module to handle DS9 regions utilities"""


def get_boxregion_string(x, y, a, b, PA):
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
    PA : float
        Position angle in degrees of the box.

    """

    s = 'box({:.8f},{:.8f},{},{},{:.8f}) ||\n'.format(x, y, a, b, PA)
    return s


def create_regfile_nboxes_from_slit(output, xy, PA, l, w, n):
    """

    Based on F. Corro script

    Parameters
    ----------
    output
    PA
    xy
    l
    w
    n : int
        Odd integer

    Returns
    -------

    """
    if (n%2 == 0) or (n < 0):
        raise ValueError("n must be odd positive integer!")

    f = open(output, "w+")
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('physical\n')

    nside = (n - 1)/2
    a = w
    b = l / n
    pa_rad = PA * np.pi / 180.
    for ii in range(-1*nside, nside+1):
        xnew = xy[0] + ii * b * np.sin(pa_rad)
        ynew = xy[1] + ii * b * np.cos(pa_rad)
        s = get_boxregion_string(xnew, ynew, a, b, PA)
        f.write(s)
    f.close()

