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


def get_nocont_cube(cube, order=1, nsig=(-2.0,2.0), inspect=False, verbose=False):
    """
    Substracts spectral continuum to a cube spaxels per spaxel
    Only for unmasked spaxels

    Parameters
    ----------
    cube : mpdaj.obj.Cube
        Cube object
    order : int
        Order of the polinomyal to fit the continuum
        This is fed into mpdaf.obj.Spectrum.poly_spec() method
    nsig : (float, float)
        The low and high rejection factor in std units (-3.0,3.0)
        This is fed into mpdaf.obj.Spectrum.poly_spec()
    inspect : bool
        Whether to inspect the continuum fits
    verbose : bool


    Returns
    -------
    cube_nocont : Cube
        Copy of Cube with continuum substracted

    """
    cube_new = cube.copy()
    nw, ny, nx = cube.shape
    # subtract continuum on every spectrum
    q = 0
    ntot = nx*ny
    for i in range(nx):
        for j in range(ny):
            q += 1
            if np.alltrue(cube.mask[:,j,i]) == True:
                continue
            spec = cube[:,j,i]
            cond = spec.data == 0.
            if np.sum(cond) == len(cond):  # empty data, just continue
                cube_new[:, j, i] = spec
                continue
            try:
                cont = spec.poly_spec(order, nsig=nsig)
            except:
                import pdb; pdb.set_trace()
            s = 'Spaxel ({},{}) [{}/{}]'.format(i,j,q,ntot)
            if verbose:
                print(s)
            if inspect:
                spec.plot(title=s)
                cont.plot(color='r', linestyle='solid')
                answer = input("Do you want to continue? [y/n]: ")
                if answer in ["y", "Y", 'yes']:
                    pass
                else:
                    raise RuntimeError('You decided not to continue.')

                plt.show()
                plt.clf()
            spec_n = spec.copy()
            spec_n.data = spec.data - cont.data.data
            cube_new[:,j,i] = spec_n

    return cube_new


def cube_ima2abs(cube_imag, pixelscale=0.2*u.arcsec, arc_name='PSZ1GA311_G1', verbose=False):
    """

    Parameters
    ----------
    cube_imag : Cube
        Cube in the image plane

    arc_name: str
        Name of arc

    Returns
    -------
    cube_abs : Cube
        Cube resampled into the absorber plane

    """
    from arctomo.lensing import cube_ima2abs
    return cube_ima2abs(cube_imag, pixelscale=pixelscale, arc_name=arc_name, verbose=verbose)


    # first loop for establishing the new radec range
    nw, ny, nx = cube_imag.shape
    ntot = nx * ny

    # Now determine the physical side as the original image to
    # use in the new one. De-lensing should make the real image smaller
    dy, dx = cube_imag.wcs.get_axis_increments() # note order (y,x)
    # import pdb; pdb.set_trace()
    assert np.allclose(np.fabs(dx), np.fabs(dy), rtol=1e-05, atol=1e-08, equal_nan=False),\
        'The cube has different increments for x and y. Not implemented for this.'
    orig_pixscale = (np.fabs(dx)*u.deg).to('arcsec')  # usually dx is negative
    factor = orig_pixscale.value / pixelscale.to('arcsec').value
    nside = int(factor*nx/2.)

    ra_abs = []
    dec_abs = []
    specs = [] # will store mpdaf Spec objects
    q = 1
    for x in range(nx):
        for y in range(ny):
            dec, ra = cube_imag.wcs.pix2sky((y, x))[0]
            ra_new, dec_new = ima2abs(ra, dec, arc_name=arc_name)
            if verbose:
                pos1 = SkyCoord(ra,dec, unit='deg')
                pos2 = SkyCoord(ra_new, dec_new, unit='deg')
                # sep = pos1.separation(pos2)
                sep_ra, sep_dec = pos1.spherical_offsets_to(pos2)

                print('Spaxel ({},{}) [{}/{}]'.format(x, y, q, ntot))
                print("  ra:  {:.6f} --> {:.6f} (Delta={:.1f}arcsec)".format(ra, ra_new, sep_ra.arcsec))
                print("  dec: {:.6f} --> {:.6f} (Delta={:.1f}arcsec)".format(dec, dec_new, sep_dec.arcsec))
            ra_abs += [ra_new]
            dec_abs += [dec_new]
            spec_aux = cube_imag[:, y, x]  # check order (y,x)
            specs += [spec_aux]
            q += 1

    # esablish the (ra,dec) center of new cube
    ra_center = 0.5 * (np.min(ra_abs) + np.max(ra_abs))
    dec_center = 0.5 * (np.min(dec_abs) + np.max(dec_abs))
    radec_center = (ra_center, dec_center) * u.deg
    cube_abs = make_empty_cube(radec_center, pixelscale, nside=nside, wave_coord=cube_imag.wave)
    # mask out all values
    cube_abs.mask[:,:,:] = True

    # fill in the specs
    counter = np.zeros_like(cube_abs.data[0].data)
    for ii, spec in enumerate(specs):
        # ynew, xnew = cube_abs.wcs.sky2pix((dec_abs[ii],ra_abs[ii]))[0]
        ynew, xnew = cube_abs.wcs.sky2pix((dec_abs[ii], ra_abs[ii]), nearest=True)[0]
        # unmask first, otherwise values are not assigned
        cube_abs.mask[:,ynew, xnew] = False
        cube_abs.data[:, ynew, xnew] += specs[ii].data
        cube_abs.var[:, ynew, xnew] += specs[ii].var
        # counter
        counter[ynew, xnew] += 1
    # are we loosing information ?
    if np.max(counter) > 1:
        cond = counter > 1
        nbad = np.sum(cond)
        print(
        "Warning(NT): there are {}/{} pixels in the abs-plane that have contributions from two or more pixels in the image plane. You should consider reducing the image scale (current = {} arcsec) of the new datacube.".format(
            nbad, ntot, pixelscale.to('arcsec').value))
    # crop cube
    cube_abs.crop()
    return cube_abs, counter


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

