"""Module for dealing with Magellan/Mage data"""

"""This scripts creates a "cube" version for a single MagE slit.
(c) NT 2018"""

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
from astropy.io.fits.hdu.table import BinTableHDU
from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU
from astropy.utils import isiterable
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra.io import readspec
from linetools.spectra.utils import collate
import copy
import json
from pyntejos import utils as ntu
from pyntejos import regions as ntr
from mpdaf.obj import Cube, WCS, WaveCoord

usage = """\
Usage: make_mage_cube.py config_file.json
"""

def read_single_mage_file(filename, **kwargs):
    """Reads single mage file produced by MIDAS (SLopez format) and returns
    a XSpectrum1D object"""

    # hd = hdulist[0].header
    # uwave = setwave(hd) * u.Angstrom
    spec = ascii.read(filename, data_start=3, data_end=-1, header_start=None)
    # import pdb;pdb.set_trace()
    wv = spec['col2'].data * u.Angstrom
    fl = spec['col3'].data
    try:
        sig = spec['col4'].data
    except:
        sig = np.ones_like(fl)

    xspec1d = XSpectrum1D.from_tuple((wv, fl, sig), **kwargs)
    return xspec1d


def get_CD_matrix(pixscale_x, pixscale_y, pos_angle):
    """gets CD matrix based on pixscale_x and pixelscale_y and position angle in degrees
    pixscale_x : size of pixel x in degrees
    pixscale_y : size of pixel y in degrees
    pos_angle : position angle in degrees
    """
    ang_rad = pos_angle * np.pi / 180.
    # val1 is defined to_east
    # val2 is defined to_north
    cd1_1 = pixscale_x * -1 * np.cos(ang_rad)  # proj. x on val1
    cd1_2 = pixscale_y * 1 * np.sin(ang_rad) # proj. y on val1
    cd2_1 = pixscale_x * 1 * np.sin(ang_rad) # proj. x on val2
    cd2_2 = pixscale_y * 1 * np.cos(ang_rad) # proj. y on val2
    # with this rotation the (x,y)=(1,11) is the northern position of the slit for pos_angles 0<PA<180
    # return values
    return cd1_1, cd1_2, cd2_1, cd2_2


def create_ref_hdulist(reference_muse_file, reference_mage_file, params):
    """Create the HDUList by taking a muse file as reference and updating the header
    by keeping/removing relevant/irrelevant information. Uses the MUSE one as the
    template"""

    # load the reference HDULists
    hdulist_muse = fits.open(reference_muse_file)
    hdulist_mage = fits.open(reference_mage_file)

    # clean/update MUSE Headers
    # primary header
    pr_hdr = hdulist_muse[0].header
    # update
    keys_mage = ['ORIGIN', 'TELESCOP','INSTRUME','OBJECT', 'DATE-OBS',   'AIRMASS', 'SLITNAME']
    for key in keys_mage:
        pr_hdr[key] = hdulist_mage[0].header[key]
    # new keys
    pr_hdr['DATE'] = Time.now().iso  # current time
    pr_hdr['EXPTIME'] = params['REDUCED_EXPTIME']
    pr_hdr['RA'] = hdulist_mage[0].header['RA-D']
    pr_hdr['DEC'] = hdulist_mage[0].header['DEC-D']
    pr_hdr['ARCFILE'] = "N/A"
    pr_hdr['DATAMD5'] = "N/A"
    pr_hdr['LST'] = 'XXX'
    pr_hdr['MJD-OBS'] =  'XXX'
    if params['POS_ANGLE'] != "None":
        pr_hdr['POS_ANGLE'] = params['POS_ANGLE']
    else: # get it from mage header
        pr_hdr['POS_ANGLE'] = hdulist_mage[0].header['ROTANGLE'] - 44.5  # this is the current offset in angle
    pr_hdr['PIXSCALE_X'] = params['PIXSCALE_X']
    pr_hdr['PIXSCALE_Y'] = params['PIXSCALE_Y']

    # compute rotation matrix
    cd1_1, cd1_2, cd2_1, cd2_2 = get_CD_matrix(pr_hdr['PIXSCALE_X'], pr_hdr['PIXSCALE_Y'], pr_hdr['POS_ANGLE'])

    for key in pr_hdr.keys():
        if "ESO " in key:  # remove ESO Headers
            del pr_hdr[key]

    # extensions
    hdu_1 = hdulist_muse[1].header
    hdu_2 = hdulist_muse[2].header
    for hdu in [hdu_1, hdu_2]:
        hdu['CD3_3'] = params['CD3_3']
        hdu['CRPIX3'] = params['CRPIX3']
        hdu['CRVAL3'] = params['CRVAL3']
        # rotation
        hdu['CD1_1'] = cd1_1
        hdu['CD1_2'] = cd1_2
        hdu['CD2_1'] = cd2_1
        hdu['CD2_2'] = cd2_2
        # reference pixel
        hdu['CRPIX1'] = params['CRPIX1']
        hdu['CRPIX2'] = params['CRPIX2']
        hdu['CRVAL1'] = params['CRVAL1']
        hdu['CRVAL2'] = params['CRVAL2']

        # delete irrelevant keys
        for key in hdu.keys():
            if "ZAP" in key:
                del hdu[key]
            if key == 'COMMENT':
                if "ZAP" in hdu['COMMENT']:
                    del hdu['COMMENT']
        # update units
        hdu['BUNIT'] = params['BUNIT']
    # return
    return hdulist_muse

def write_config_make_MagE_cube_dummy(filename):
    """Dumps a dummy configuration file for mage.make_MagE_cube()"""
    text ="""
{
    "directory_mage" : "Directory of the MagE data for a single slit with format and naming convention of S. Lopez",
    "reference_mage" : "Filename of a MagE .fits image of reference of the observations. Headers will be updated from this",
    "reference_muse" : "Filename of a MUSE .fits cube of reference, to define the relevant keywords of the MagE cube.",
    "output_cube" : "Name of output .fits file of the MagE datacube (e.g. cube_mage.fits)",
    "CRPIX1" :  1,
    "CRPIX2" :  1,
    "CRVAL1" : "RA coordinate in degrees of cube spaxel (1,1). Must be float!",
    "CRVAL2" : "Dec coordinate in degrees of cube spaxel (1,1). Must be float!",
    "PIXSCALE_Y" : "Pix scale along the slit in degrees/pixel. Must be float!",
    "PIXSCALE_X" : "Pix scale transverse to the slit direction in degrees/pixel. Must be float!",
    "POS_ANGLE" : "Position angle of the slit in degrees. Use string None if the angle is read from the header of the reference MagE file",
    "CD3_3"  :  "Delta lambda. Must be float!",
    "CRPIX3" :  1,
    "CRVAL3" :  "Lambda_0. Must be float!",
    "BUNIT" : "Unit of flux, as string: e.g. XXX *erg/s/cm**2/Angstrom",
    "REDUCED_BY" : "Name of person who performed the reduction.",
    "REDUCED_DATE" : "Date of the reduction",
    "REDUCED_EXPTIME" : "Total exposure time of the cube, from the reduction. Must be int or float!"
}"""
    f = open(filename, 'w')
    f.write(text)
    f.close()
    # with open(filename, 'w') as fp:
    #     json.dump(params, fp, indent=4)


def make_MagE_cube_old(config_file):
    """Creates a IFU cube for MagE data of a single slit MagE

    Old version copies the MUSE header, new version starts from scratch using MPDAF classes

    Parameters
    ----------
    config_file : str
        Configuration file with important parameters. See mage.dump_config_make_MagE_cube()

    """

    # read config_filename and add filename to params
    f = open(config_file)
    params = json.load(f)
    params['config_filename'] = config_file

    # reference files
    reference_mage_file = params['reference_mage']
    reference_muse_file = params['reference_muse']

    # Add comment on how the header was created
    comment = 'Cube created by pyntejos.mage.make_MagE_cube.py based on headers from {} and {}. Astrometry ' \
              'and wavelength axis given by input configuration ' \
              'file {}.'.format(reference_mage_file.split('/')[-1],reference_muse_file.split('/')[-1], params['config_filename'])

    # Print message
    print('')
    print(comment.replace("Cube created", "MagE cube will be created"))

    dirname = params['directory_mage']
    filenames = glob.glob(dirname+"/*_??.fits")
    # sort them
    filenames.sort()

    # create the new datacube structure
    len_wv = readspec(filenames[0]).npix
    cube = np.zeros_like(np.ndarray(shape=(len_wv, len(filenames), 1)))
    stat = np.zeros_like(cube)
    # model = np.zeros_like(cube)

    # read the files and fill the data cubes (cube and stat)
    print("Reading files from directory {} ordered as:".format(dirname))
    ny = len(filenames)
    for ii,fname in enumerate(filenames):
        fn = fname.split('/')[-1]
        print("\t {}: {}".format(ii + 1, fn))
        # MAKE SURE THE SPECTRA IS PROPERLY SORTED
        nfile = fn.split('.')[0].split('_')[-1]
        nfile = int(nfile)
        assert nfile == ii + 1, "The files in the directory are not sorted properly. Please check."
        try:
            spec = readspec(fname)
            if 0: # dummy
                spec.flux = 1
                spec.sig = 0
            cube[:,0,ny-ii-1] = spec.flux.value
            if np.isnan(spec.sig):
                try:
                    spec_sig = readspec(fname.replace('.fits','_sig.fits'))
                    spec.sig = spec_sig.flux.value
                except:
                    spec.sig = 0
            stat[:,0,ny-ii-1] = (spec.sig.value)**2
            # model[:,ii,0] = model_sp.flux.value

            if 0: # dummy
                spec.flux = 1
                spec.sig = 0
            cube[:,0,ny-ii-1] = spec.flux.value
            stat[:,0,ny-ii-1] = (spec.sig.value)**2
            # model[:,ii,0] = model_sp.flux.value
        except:
            print("Something is wrong with spectrum {}".format(fname.split('/')[-1]))
            import pdb; pdb.set_trace()
            raise ValueError("Something is wrong with spectrum {}".format(fname.split('/')[-1]))

    # get the reference hdulist
    hdulist_new = create_ref_hdulist(params['reference_muse'], params['reference_mage'], params)

    # update the data
    hdulist_new[1].data = cube
    hdulist_new[2].data = stat

    # add model HDU
    # hdu_model = copy.deepcopy(hdulist_new[2])
    # hdu_model.data = model
    # hdu_model.header['EXTNAME'] = 'MODEL'
    # hdulist_new.append(hdu_model)

    # write the cube
    hdulist_new.writeto(params['output_cube'], clobber=True)
    print('Wrote file: {}'.format(params['output_cube']))


def make_MagE_cube_v2(config_file):
    """A new version to create the MagE cubes, it uses MPDAF objects.
    It works for .fits 1-d spectra files stored in a directory
    Errors are looked for in the same directory with *_sig.fits extension


    Parameters
    ----------
    config_file : str
        Configuration file with important parameters. See mage.dump_config_make_MagE_cube()

    """

    # read config_filename and add filename to params
    f = open(config_file)
    params = json.load(f)
    params['config_filename'] = config_file

    # reference files
    reference_mage_file = params['reference_mage']
    reference_muse_file = params['reference_muse']

    # Add comment on how the header was created
    comment = 'Cube created by pyntejos.mage.make_MagE_cube.py based on headers from {} and {}. Astrometry ' \
              'and wavelength axis given by input configuration ' \
              'file {}.'.format(reference_mage_file.split('/')[-1],reference_muse_file.split('/')[-1], params['config_filename'])

    # Print message
    print('')
    print(comment.replace("Cube created", "MagE cube will be created"))

    dirname = params['directory_mage']
    # filenames = glob.glob(dirname+"/*syn_??.fits")
    filenames = glob.glob(dirname+"/*_bin_flux_??.fits")

    # sort them
    filenames.sort()

    # create the new datacube structure
    spec_ref = readspec(filenames[0])
    nw, ny, nx = spec_ref.npix, len(filenames), 1

    # get wavecoord
    hdr = spec_ref.header
    hdr['CUNIT1'] = 'Angstrom'
    wave = WaveCoord(hdr=hdr)
    # get WCS
    crpix_yx = params['CRPIX2'], params['CRPIX1']  # note y,x order
    crval_decra = params['CRVAL2'], params['CRVAL1']  # ditto
    cdelt_yx = params['PIXSCALE_Y'], -1*params['PIXSCALE_X']  # ditto (negative to decrease towards east)
    PA = params['POS_ANGLE']
    if PA == "None": # get angle from MagE reference header
        print("No PA given in parameter file. Reding PA from MagE reference .fits header")
        hdulist_mage = fits.open(params['reference_mage'])
        PA = hdulist_mage[0].header['ROTANGLE'] - 44.5  # this is the current offset in angle
        print("PA={}deg".format(PA))
    shape = (ny, nx)
    wcs = WCS(crpix=crpix_yx, crval=crval_decra, cdelt=cdelt_yx, deg=True, shape=shape, rot=-1*PA)  # for some reason negative PA works
    # redefine wcs to have CD matrix rather than PC matrix
    if 1: #rot != 0:
        hdr = wcs.to_header()
        hdr.rename_keyword('PC1_1','CD1_1')
        hdr.rename_keyword('PC2_1','CD2_1')
        hdr.rename_keyword('PC1_2','CD1_2')
        hdr.rename_keyword('PC2_2','CD2_2')
        # hdr['CD1_1'] = hdr['PC1_1']
        # hdr['CD2_1'] = hdr['PC2_1']
        # hdr['CD1_2'] = hdr['PC1_2']
        # hdr['CD2_2'] = hdr['PC2_2']
        wcs = WCS(hdr=hdr)

    # create data structures
    data = np.zeros_like(np.ndarray(shape=(nw,ny,nx)))
    var = np.zeros_like(data)

    # read the files and fill the data cubes (cube and variance)
    print("Reading files from directory {} ordered as:".format(dirname))
    for ii,fname in enumerate(filenames):
        fn = fname.split('/')[-1]
        yspaxel = ny - ii  # larger y will be the northern one
        print("\t {}: {} --goes to--> spaxel (1,{})".format(ii + 1, fn, yspaxel))
        # MAKE SURE THE SPECTRA IS PROPERLY SORTED
        nfile = fn.split('.')[0].split('_')[-1]
        nfile = int(nfile)
        assert nfile == ii + 1, "The files in the directory are not sorted properly. Please check."
        try:
            spec = readspec(fname)
            if 0: # dummy
                spec.flux = 1
                spec.sig = 0
            data[:,yspaxel-1,0] = spec.flux.value  # position "01" is the most north, y increase to north
            if np.isnan(spec.sig):
                try:
                    spec_sig = readspec(fname.replace('flux','sigm'))
                    spec.sig = spec_sig.flux.value
                except:
                    spec.sig = 0
            var[:,yspaxel-1,0] = (spec.sig.value)**2
        except:
            print("Something is wrong with spectrum {}".format(fname.split('/')[-1]))
            import pdb; pdb.set_trace()
            raise ValueError("Something is wrong with spectrum {}".format(fname.split('/')[-1]))

    # create the Cube object
    cube = Cube(data=data, var=var, wcs=wcs, wave=wave)
    cube.write(params['output_cube'])
    print('Wrote file: {}'.format(params['output_cube']))

    # create white image
    white = cube.sum(axis=0)
    white.write(params['output_cube'].replace('cube', 'white'))
    print('Wrote file: {}'.format(params['output_cube'].replace('cube', 'white')))
    # import pdb; pdb.set_trace()
    return cube


def plot_specs_from_magecube(magecube, only_plot=None, **kwargs):
    """Plots the n spectra from a magecube

    magecube : mpdaf Cube
        original cube
    only_plot : iterable or None
        If not None, iterable with positions to plot. e.g. [1,3,5]
        Convention is that positions start from y=1 upwards

    """

    nw, ny, nx = magecube.shape
    assert nx == 1, "Your magecube does not have the conventional astrometry, where the slit is aligned in the y-axis"
    if only_plot is not None:
        pos_array = only_plot
    else:
        pos_array = range(1, ny+1)
    for ii in pos_array:
        sp = magecube[:,ii-1,0]
        spec = ntu.xspectrum1d_from_mpdaf_spec(sp)
        plt.plot(spec.wavelength, spec.flux, drawstyle='steps-mid', label='(x,y)=({},{})'.format(1,ii), **kwargs)
    plt.xlabel('Wavelength (AA)')
    plt.ylabel('Relative flux')
    plt.legend()
    plt.show()


def write_magecube_as_xspectrum1d(magecube_filename):
    magecube = Cube(magecube_filename)
    nw, ny, nx = magecube.shape
    assert nx == 1, "Your magecube does not have the conventional astrometry, where the slit is aligned in the y-axis"

    spec_list = []
    for ii in range(ny):
        sp = magecube[:,ii-1,0]
        spec = ntu.xspectrum1d_from_mpdaf_spec(sp)
        spec_list += [spec]
    specs = collate(spec_list)
    new_name = magecube_filename.replace('.fits','_xspec.fits')
    specs.write(new_name)
    return specs


def compute_chi2_magecubes(magecube1, magecube2, chi2_wvrange, renorm_wvrange, plot_wvrange, plot=False,
                           text1='musecube1', text2='musecube2'):
    """
    Computes a (renormalized) spectral Chi2 analysis spaxel per spaxel betweein magecube1 and magecube2.
    Computes a (normalized) total flux Chi2 analysis spaxel per spaxel between magecube1 and magecube2.

    magecubes must have the same spatial shape (in xy pixels), but can have different WaveCoord
    magecube1 is rebinned spectrally to match the WaveCoord of magecube2 within the chi2_wvrange


    Parameters
    ----------
    magecube1 : mpdaf.obj.Cube
        Magecube 1
    magecube2 : mpdaf.obj.Cube
        Magecube 2, must be of same shape as magecube1
    chi2_wvrange : (float, float)
        Wavelength range for performing the chi2 comparison. Ideally it has to be a small
        region centred in a spectral feature (e.g. emission line)
    renorm_wvrange : (float, float)
        Wavelength range for performing the flux renormalization of MagE to MUSE. Ideally it
        has to be a small region close to `chi2_range`.
    plot_wvrange : (float, float)
        Wavelength range for plotting. Must be larger than either renorm_range or chi2_range.
    plot : bool
        Whether to plot the iterations for visual inspection

    chi2_range : (float, float)
        Wavelenght range to perform chi2 computation

    Returns
    -------
    chi2_spec, chi2_flux, flux_mage1, flux_mage2 : float, float, np.array(), np.array()

    """

    assert magecube2.shape[1:] == magecube1.shape[1:], 'magecube2 must have the same spatial shape as magecube1'
    nw, ny, nx = magecube1.shape
    assert nx == 1, 'Magecubes must have only 1 spaxel in the x-direction'

    chi2_spec = []
    fl_mage1 = []
    fl_mage2 = []
    if plot:
        # prepare the plot layout (1 panel top, ny panels bottom)
        fig = plt.figure(figsize=(10,10))
        axes = [fig.add_subplot(411)] # top panel
        ax_2r = [fig.add_subplot(4,4,4+ii+1) for ii in range(4)]
        ax_3r = [fig.add_subplot(4,4,8+ii+1) for ii in range(4)]
        ax_4r = [fig.add_subplot(4,4,12+ii+1) for ii in range(3)]
        axes = axes + ax_2r + ax_3r + ax_4r

    for ii in range(ny):
        sp1 = magecube1[:, ii, 0]  # Spectrum
        sp2 = magecube2[:, ii, 0]  # Spectrum
        # keep only region of interest
        for sp in [sp1,sp2]:
            sp.mask_region(lmin=plot_wvrange[0], lmax=plot_wvrange[1], inside=False)
        # convert to XSpectrum1D objects
        spec1 = ntu.xspectrum1d_from_mpdaf_spec(sp1)
        spec2 = ntu.xspectrum1d_from_mpdaf_spec(sp2)
        # rebin and renorm mage to muse
        if spec1.wavelength != spec2.wavelength:
            spec3 = spec1.rebin(spec2.wavelength, do_sig=True)  # rebin in case is necessary
        else:
            spec3 = spec1
        renorm2 = ntu.renorm2_factor_from_wvrange(spec2, spec3, wvrange=renorm_wvrange)
        spec3.data['flux'] = spec3.flux * renorm2
        spec3.data['sig'] = spec3.sig * renorm2

        # plot original MagE data in grey
        if plot:
            ax = axes[ii+1]
            ax.plot(spec1.wavelength, renorm2*spec1.flux, drawstyle='steps-mid', c='gray', label='MagE-1 original',
                    lw=3, alpha=0.25)
            # plot the data to compare
            ax.vlines(chi2_wvrange, ymin=0, ymax=10 * np.max(spec2.flux), color='r')  # plot the chi2 limits ranges
            ntu.plot_two_spec(spec2, spec3, text1=text1, text2=text2, ax=ax)

        # compute chi2
        chi2_aux = ntu.chi2_from_two_spec(spec2, spec3, wvrange=chi2_wvrange)
        chi2_spec += [chi2_aux]

        # also keep the total fluxes per spaxel
        cond = (spec2.wavelength.to('AA').value >= chi2_wvrange[0]) & (spec2.wavelength.to('AA').value <= chi2_wvrange[1])
        fl_mage2 += [np.nansum(spec2.flux[cond])]
        cond = (spec1.wavelength.to('AA').value >= chi2_wvrange[0]) & (spec1.wavelength.to('AA').value <= chi2_wvrange[1])
        fl_mage1 += [np.nansum(spec1.flux[cond])]

        if plot:
            # ax.set_title('Pos#{}, Chi2_spec={:.1f}'.format(ii + 1, chi2_aux))
            ax.set_xlim(plot_wvrange[0], plot_wvrange[1])

    chi2_spec = np.array(chi2_spec)
    chi2_spec = np.sum(chi2_spec) / ny
    fl_mage1 = np.array(fl_mage1)
    fl_mage2 = np.array(fl_mage2)
    renorm2 = fl_mage1[0] / fl_mage2[0]
    fl_mage1 = fl_mage1
    fl_mage2 = fl_mage2 * renorm2
    chi2_flux = np.sum((fl_mage1 - fl_mage2) ** 2 / fl_mage2**2)
    # import pdb; pdb.set_trace()
    if plot:
        legend = axes[-1].get_legend()
        ax = axes[0]
        ax.plot(fl_mage2, 'ko-', drawstyle='steps-mid')
        ax.plot(fl_mage1, 'bo-', drawstyle='steps-mid')
        ax.set_title('Chi2_flux={:.1f}, Chi2_spec={:.1f}'.format(chi2_flux, chi2_spec))
        ax.legend_ = legend
        ax.legend()

    return chi2_spec, chi2_flux, fl_mage1, fl_mage2


def determine_best_astrometry(magecube_filename, musecube_filename, xc_array, yc_array, PA_array,
                              chi2_wvrange, renorm_wvrange, plot_wvrange, plot=True):
    """Utility for determining the best astrometry for a MagE Cube based on a comparisom
    with a reference MUSE Cube. It will use different (xc, yc, PA) positions for a virtual MagE slit
    on a MUSE reference datacube and will determine the set that has the minimum Chi2 within a
    spectral region. It will return a astropy.table.Table object with the results.

    [Warning: it may be time consuming scales as len(xc_array)*len(yc_array)*len(pa_array)]

    Parameters
    ----------
    magecube : str
        Filename of the MagE cube
    musecube_filename : str
        Filename of the MUSE reference cube
    xc_array : 1-d np.array
        Array with the different xc position to try for the MagE slit
        In units of MUSE pixels
    yc_array : 1-d np.array
        Array with the different yc position to try for the MagE slit
        In units of MUSE pixels
    PA_array : 1-d np.array
        Array with the different position angles (PA) to try for the MagE slit
    chi2_wvrange : (float, float)
        Wavelength range for performing the chi2 comparison. Ideally it has to be a small
        region centred in a spectral feature (e.g. emission line)
    renorm_wvrange : (float, float)
        Wavelength range for performing the flux renormalization of MagE to MUSE. Ideally it
        has to be a small region close to `chi2_range`.
    plot_wvrange : (float, float)
        Wavelength range for plotting. Must be larger than either renorm_range or chi2_range.
    plot : bool
        Whether to plot the iterations for visual inspection

    Returns
    -------
    results : astropy.table.Table
        Table with the results
    """
    from PyMUSE import musecube as mc
    import os

    # load the cubes
    musecube_pymuse = mc.MuseCube(musecube_filename)
    musecube = Cube(musecube_filename)
    magecube_orig = Cube(magecube_filename)

    # create master subdir
    master_dirname = magecube_filename.replace('.fits', '_astrometry_runs')
    if not os.path.exists(master_dirname):
        os.makedirs(master_dirname)

    # define mage slit geometry
    l=0.3*3*11
    w=1.
    n=11
    pixscale = 0.2  # of MUSE datacube

    # create a table to store the results
    tab = Table()

    # will create auxiliary cubes from MUSE to match geometry of mage
    names = []
    for xc in xc_array:
        for yc in yc_array:
            for pac in PA_array:
                rootname = '{:.1f}-{:.1f}-{:.1f}'.format(xc,yc,pac)
                rootname = rootname.replace('.','p')
                names += [rootname]

                if not os.path.exists(master_dirname+'/'+rootname):
                    os.makedirs(master_dirname+'/'+rootname)

                output = master_dirname + '/' + rootname + '/mage_slit_in_muse_{}.reg'.format(rootname)
                ntr.create_regfile_nboxes_from_slit(output, (xc,yc), pac, l, w, n, pixscale)

                # define name
                newcube_name = master_dirname + '/' + rootname + '/magecube_from_muse_{}.fits'.format(rootname)
                # check whether this cube already exists
                if os.path.isfile(newcube_name):
                    print('File already exists, skiping: {}'.format(newcube_name.split('/')[-1]))
                    continue
                else:
                    print('Writing file: {}'.format(newcube_name))
                # create the data structure
                nw, ny, nx = magecube_orig.shape
                nw = musecube.shape[0] # replace mw with that of muse datacube
                data = np.zeros(shape=(nw, ny, nx))
                var = np.zeros_like(data)

                # get the spectra from MUSE in the given region
                for ii in range(n):
                    plt.figure(musecube_pymuse.n)
                    spec = musecube_pymuse.get_spec_from_ds9regfile(output, mode='sum', i=ii, frac=0.1, npix=0,
                                                          empirical_std=False, n_figure=None,
                                                            save=False, save_mask=False, plot=False)
                    data[:, ny-ii-1, 0] = spec.flux
                    var[:, ny-ii-1, 0] = spec.sig**2
                    # import pdb; pdb.set_trace()
                # redefine the structure
                magecube_new = Cube(wcs=magecube_orig.wcs, wave=musecube.wave, data=data, var=var)
                newcube_name = master_dirname + '/' + rootname + '/magecube_from_muse_{}.fits'.format(rootname)
                magecube_new.write(newcube_name)
    # plt.show()
    tab['name'] = names
    tab['xc'] = [float(name.split('-')[0].replace('p','.')) for name in tab['name']]
    tab['yc'] = [float(name.split('-')[1].replace('p','.')) for name in tab['name']]
    tab['PA'] = [float(name.split('-')[2].replace('p','.')) for name in tab['name']]

    # now read the cubes and perform the chi2
    chi2_spec = []
    chi2_flux = []
    fl_mage_11 = []
    fl_muse_11 = []
    for jj, name in enumerate(tab['name']):
        newcube_name = master_dirname + '/' + name + '/magecube_from_muse_{}.fits'.format(name)
        magecube1 = magecube_orig  # MagE original
        magecube2 = Cube(newcube_name)  # from MUSE
        chi2_s, chi2_f, fl_mage, fl_muse = compute_chi2_magecubes(magecube1, magecube2,
                                                    chi2_wvrange=chi2_wvrange,
                                                    renorm_wvrange=renorm_wvrange,
                                                    plot_wvrange=plot_wvrange, plot=plot,
                                                    text1= 'MagE', text2='MUSE')

        print("{}  [{}/{}]".format(name, jj + 1, len(tab)))
        print(" Total Chi2_spec/DOF = {:.1f}".format(chi2_s))
        print(" Total Chi2 flux/DOF = {:.1f}".format(chi2_f))
        chi2_spec += [chi2_s]
        chi2_flux += [chi2_f]
        fl_mage_11 += [fl_mage]
        fl_muse_11 += [fl_muse]
        if plot:
            fig = plt.gcf()
            fig.suptitle(name)

    tab['fl_muse_11'] = fl_muse_11
    tab['fl_mage_11'] = fl_mage_11
    tab['chi2_spec'] = chi2_spec
    tab['chi2_flux'] = chi2_flux
    return tab


def overwrite_magecube(cube, fl_arrays, sig_arrays):
    """For a given magecube cube, it will rewrite the content of
    flux and sigma arrays.

    Parameters
    ----------
    cube: Cube
        a cube with MagE data
        11-positions, where position 1 is Southernmost   one when PA=0
    fl_arrays : list of np.arrays()
        List of 11 arrays containing the fluxes in each position
        First element corresponds to the Northernmost one
    sig_arrays : list of np.arrays()
        List of 11 arrays containing the flux uncertainties in each position
        First element corresponds to the Northernmost one

    Returns
    -------
    newcube : Cube
        The new overwritten cube

    """
    newcube = cube.copy()
    assert len(fl_arrays) == 11, "The fl_arrays list must be of 11 elements, each one being a np.array()"
    assert len(sig_arrays) == 11, "The sig_arrays list must be of 11 elements, each one being a np.array()"

    nw, ny, nw = cube.shape
    for ii in range(len(fl_arrays)):
        yspaxel = ny - ii  # larger y will be the northern one
        newcube.data.data[:,yspaxel - 1, 0] = np.array(fl_arrays[ii])
        newcube.var.data[:, yspaxel - 1, 0] = np.array(sig_arrays[ii])**2  # this is variance
    return newcube
