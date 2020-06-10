"""This module is meant for analysis of VLT/MUSE data"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from mpdaf.obj import Image, Cube, WCS, WaveCoord


def get_subcube(inp, xy_center, nside, wv_range=None):
    """

    Parameters
    ----------
    inp : str or Cube
        Input filename or Cube object
    xy_center : (float, float)
        Coordinates of the center in pixels
    nside : int
        Number of pixels to each side of the central one
        for the new cube
    wv_range: (float, float)
        Wavelength range of the subcube (wvmin,wvmax), must be in angstroms
        If None it gets the full spectral range

    Returns
    -------
    cube : Cube
        Subcube

    """

    if not isinstance(inp, Cube): # assume str
        cube = Cube(inp)
    else:
        cube = inp
    yx_center = (xy_center[1], xy_center[0]) # note that center in mpdaf is given as (y,x)
    unit_center = None # to be interpreted as pixels
    subcube = cube.subcube(yx_center, 2*nside+1, lbda=wv_range,
                           unit_center=unit_center,
                           unit_size=unit_center,
                           unit_wave=u.Unit("Angstrom"))

    return subcube

def masked_cube(inp, mask, value=0., method='value'):
    """Returns a masked version of the input cube

    inp : str or Cube
        Input filename or Cube object
    mask : boolean np.array 2D of same shape (ny, nx) than cube spatial
        Pixels to mask
    value : float
        Values to replaces masked spaxels
        if method = 'value'
    method : str
        Method to perform the masking
        value : replaces by the given value
        madian : replaces by the median of all contiguous pixels including the value of the pixel itself

    """
    if not isinstance(inp, Cube): # assume str
        cube = Cube(inp)
    else:
        cube = inp
    newcube = cube.copy()
    nw, ny, nx = cube.shape
    for ii in range(nw):
        image = cube.data[ii,:,:]
        if method == 'value':
            aux = np.where(mask, value, image)
        elif method == 'median':
            inds_y, inds_x = np.where(mask)
            aux = image.data
            for yx in zip(inds_y,inds_x):
                # import pdb;pdb.set_trace()
                median = np.median(image[yx[0]-1:yx[0]+1, yx[1]-1:yx[1]+1])
                aux[yx[0],yx[1]] = median
        newcube[ii,:,:] = aux
        # import pdb;pdb.set_trace()
    return newcube



def make_empty_cube(radec_center, pixscale, nside, wave_coord, wcs=None):
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
    wcs : mpdaf.obj.WCS ; optional
        If given, it will overwrite radec_center, pixscale, nside

    Returns
    -------
     cube : mpdaf.obj.Cube
        Cube object filled with zeros

    """
    if wcs is None:
        ntot = 2 * nside + 1 # total side always odd
        ny, nx = ntot, ntot
        crpix = (nside + 1, nside + 1)
        crval = radec_center[1].to('deg').value, radec_center[0].to('deg').value
        cdelt = pixscale.to('deg').value
        dy, dx = cdelt, -1*cdelt  # dx is negative for convention East goes to negative x's
        # cd_matrix = np.zeros((2,2))
        # cd_matrix[0,0] = cdelt
        # cd_matrix[1,1] = cdelt
        deg_bool = True
        shape = (ny, nx)
        wcs_new = WCS(crpix=crpix, crval=crval, cdelt=(dy, dx), deg=deg_bool, shape=shape)
    else:
        wcs_new = wcs
        nx = wcs.naxis1
        ny = wcs.naxis2

    nw = wave_coord.shape  #
    data = np.zeros((nw, ny, nx))
    cube = Cube(wcs=wcs_new, data=data, var=data, wave=wave_coord)
    return cube



def extend_cube(cube, floor=0., stddev=0.1):
    """Extends a given cube by adding pixels around the original cube filled with zeros or
    random values around a floor.

    Parameters
    ----------
    cube : mpdaf.obj.Cube

    floor : float
        Value for the floor
    stddev : float
        Standard deviation around the floor to generate values

    Returns
    -------
    newcube : mpdaf.obj.Cube

    """
    import random

    # shape
    nw, ny, nx = cube.shape
    wcs = cube.wcs

    # create new cubes for data and var
    new_data = np.zeros((nw, 2 * ny, 2 * nx))
    new_var = np.zeros((nw, 2 * ny, 2 * nx)) + stddev**2

    # add noise to the new cubes
    for i in range(new_data.shape[0]):
        for j in range(new_data.shape[1]):
            for k in range(new_data.shape[2]):
                new_data[i, j, k] = random.gauss(floor, stddev)

    # add original data to the new cubes
    new_data[:, int(0.5 * ny):int(1.5 * ny), int(0.5 * nx):int(1.5 * nx)] = cube._data
    new_var[:, int(0.5 * ny):int(1.5 * ny), int(0.5 * nx):int(1.5 * nx)] = cube._var

    # copy original cube and modify
    newcube = cube.copy()
    newcube._data = new_data
    newcube._var = new_var
    newcube._mask = ~new_data.astype(bool)
    newcube.wcs.set_crpix1(2 * nx)
    newcube.wcs.set_crpix2(2 * ny)
    newcube.wcs.naxis1 = 2 * nx
    newcube.wcs.naxis2 = 2 * ny
    return newcube

