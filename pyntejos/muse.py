"""This module is meant for analysis of VLT/MUSE data"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
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
    yx_center = (xy_center[1], xy_center[0]) # note that center in pmdaf is given as (y,x)
    unit_center = None # to be interpreted as pixels
    subcube = cube.subcube(yx_center, 2*nside+1, lbda=wv_range,
                           unit_center=unit_center,
                           unit_size=unit_center,
                           unit_wave=u.Unit("Angstrom"))
    return subcube


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
    dy, dx = cdelt, -1*cdelt  # dx is negative for convention East goes to negative x's
    # cd_matrix = np.zeros((2,2))
    # cd_matrix[0,0] = cdelt
    # cd_matrix[1,1] = cdelt
    deg_bool = True
    shape = (ntot, ntot)
    wcs = WCS(crpix=crpix, crval=crval, cdelt=(dy, dx), deg=deg_bool, shape=shape)
    nw = wave_coord.shape  #
    data = np.zeros((nw, ntot, ntot))
    cube = Cube(wcs=wcs, data=data, var=data, wave=wave_coord)
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
            cont = spec.poly_spec(order, nsig=nsig)
            s = 'Spaxel ({},{}) [{}/{}]'.format(i,j,q,ntot)
            if verbose:
                print(s)
            if inspect:
                spec.plot(title=s)
                cont.plot(color='r', drawstyle='line')
                plt.show()
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
    from pyntejos.lensing import ima2abs

    # first loop for establishing the new radec range
    nw, ny, nx = cube_imag.shape
    ntot = nx * ny

    # Now determine the physical side as the original image to
    # use in the new one. De-lensing should make the real image smaller
    dy, dx = cube_imag.wcs.get_axis_increments() # note order (y,x)
    # import pdb; pdb.set_trace()
    assert np.fabs(dx) == np.fabs(dy), 'The cube has different increments for x and y. Not implemented for this.'
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
                print('Spaxel ({},{}) [{}/{}]'.format(x, y, q, ntot))
                print("  ra:  {:.6f} --> {:.6f}".format(ra, ra_new))
                print("  dec: {:.6f} --> {:.6f}".format(dec, dec_new))
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




