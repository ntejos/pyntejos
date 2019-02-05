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
    filenames = glob.glob(dirname+"/*_??.fits")
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
            data[:,ny-ii-1,0] = spec.flux.value  # position "01" is the most north, y increase to north
            if np.isnan(spec.sig):
                try:
                    spec_sig = readspec(fname.replace('.fits','_sig.fits'))
                    spec.sig = spec_sig.flux.value
                except:
                    spec.sig = 0
            var[:,ny-ii-1,0] = (spec.sig.value)**2
            # model[:,ii,0] = model_sp.flux.value
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


def determine_best_astrometry(magecube_filename, musecube_filename, xc_array, yc_array, PA_array,
                              chi2_range, renorm_range, plot_range, plot=True):
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
    renorm_range : (float, float)
        Wavelength range for performing the flux renormalization of MagE to MUSE. Ideally it
        has to be a small region close to `chi2_range`.
    plot_range : (float, float)
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
    chi2 = []
    chi2w = []
    chi2_11 = []
    for name in tab['name']:
        newcube_name = master_dirname + '/' + name + '/magecube_from_muse_{}.fits'.format(name)
        # load the magecube
        # import pdb; pdb.set_trace()
        magecube_new = Cube(newcube_name)
        nw, ny, nx = magecube_new.shape

        chi2_tot = []
        s2n_tot = []
        for ii in range(ny):
            sp1 = magecube_orig[:, ii, 0]  # MagE data
            sp2 = magecube_new[:, ii, 0]  # MUSE data
            # import pdb; pdb.set_trace()
            # mask region of interest
            sp1.mask_region(lmin=plot_range[0], lmax=plot_range[1], inside=False)
            sp2.mask_region(lmin=plot_range[0], lmax=plot_range[1], inside=False)
            # xspectrum1d objects
            spec1 = ntu.xspectrum1d_from_mpdaf_spec(sp1)
            spec2 = ntu.xspectrum1d_from_mpdaf_spec(sp2)
            # rebin and renorm mage to muse
            spec3 = spec1.rebin(spec2.wavelength, do_sig=True)  # MagE rebined to MUSE
            renorm2 = ntu.renorm2_factor_from_wvrange(spec2, spec3, wvrange=renorm_range)
            spec3.data['flux'] = spec3.flux * renorm2
            spec3.data['sig'] = spec3.sig * renorm2

            # plot original MagE data in grey
            if plot:
                plt.plot(spec1.wavelength, renorm2*spec1.flux, drawstyle='steps-mid', c='gray', label='MagE original', lw=3, alpha=0.25)
                # plot the data to compare
                plt.vlines(chi2_range, ymin=0, ymax=10*np.max(spec2.flux), color='r') # plot the chi2 limits ranges
                ntu.plot_two_spec(spec2, spec3, text1='MUSE', text2='MagE')

            # compute chi2
            chi2_aux = ntu.chi2_from_two_spec(spec2,spec3, wvrange=chi2_range)
            s2n_aux = ntu.get_s2n(spec2, wvrange=chi2_range)
            chi2_tot += [chi2_aux]
            s2n_tot += [s2n_aux]
            if plot:
                plt.title('Run:{}\nPos#{}, Chi2={:.1f}, s2n={:.0f}'.format(name,ii+1, chi2_aux, s2n_aux))
                plt.legend()
                plt.xlim(plot_range[0],plot_range[1])
                plt.show()

        chi2_tot = np.array(chi2_tot)
        s2n_tot = np.array(s2n_tot)
        chi2_w = chi2_tot * s2n_tot / np.sum(s2n_tot)
        print(chi2_tot)
        print(chi2_w)
        print("Total Chi2/DOF = {:.1f}".format(np.sum(chi2_tot)/n))
        print("Weighted Chi2/DOF = {:.1f}".format(np.sum(chi2_w)))
        chi2 += [np.sum(chi2_tot)/n]
        chi2w += [np.sum(chi2_w)]
        chi2_11 += [chi2_tot]
    tab['chi2'] = chi2
    tab['chi2w'] = chi2w
    tab.write(master_dirname + '/astrometry_results.dat', format='ascii.fixed_width', overwrite=True)
    tab['chi2_11'] = chi2_11
    tab.write(master_dirname + '/astrometry_results.fits', overwrite=True)
    return tab

