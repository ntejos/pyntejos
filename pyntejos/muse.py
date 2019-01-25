"""This module is meant for analysis of VLT/MUSE data"""

import numpy as np
import astropy.units as u
from mpdaf.obj import Image, Cube, WCS, WaveCoord


def make_empty_cube(radec_center, pixscale, nside, wave_coord):
    """Makes a new datacube with different WCS and 2*nside+1
    pixels per side

    Parameters
    ----------
    radec_center : SkyCoord
        Coordinate of the central spaxel
    pixscale : Angle
        Pixel scale of the cube
    nside : int
        Number of pixels at each side of the central one
    wave_coord : mpdaf.obj.WaveCoord
        The WaveCoord object of the new datacube

    Returns
    -------
     cube : mpdaf.obj.Cube
        Cube object filled with zeros

    """
    ntot = 2 * nside + 1 # total side always odd
    crpix = (nside + 1, nside + 1)
    crval = radec_center[1].to('deg').value, radec_center[0].to('deg').value
    cdelt = pixscale.to('deg').value
    # cd_matrix = np.zeros((2,2))
    # cd_matrix[0,0] = cdelt
    # cd_matrix[1,1] = cdelt
    deg_bool = True
    shape = (ntot, ntot)
    wcs = WCS(crpix=crpix, crval=crval, cdelt=cdelt, deg=deg_bool, shape=shape)
    nw = wave_coord.shape  #
    data = np.zeros((nw, ntot, ntot))
    cube = Cube(wcs=wcs, data=data, wave=wave_coord)
    return cube
