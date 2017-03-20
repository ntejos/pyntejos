import numpy as np
import pylab as pl
from astropy.io import fits
from astropy.table import Table
from linetools.spectra.io import readspec
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra.utils import collate
import numpy as np
from pypit import arcoadd as arco
from astropy import units as u

"""Main module for co-addition of 1-d spectra"""


def coadd_stis_from_x1dfiles(filenames, wv_array=None, rebin=None):
    """

    Parameters
    ----------
    filenames : list
        List of filenames with x1d STIS data
        Must be of the same object and same
        configuration
    wv_array : Quantity array
        Wavelength array to perform the co-add
    rebin : int, optional
        If given, it rebins the current sampling by
        rebin number of pixels

    Returns
    -------
    spec1d : XSpectrum1D
        Co-added version of all the spectra
    """
    spec_list = []
    for filename in filenames:
        aux = load_single_x1d_stis(filename)
        for sp in aux:
            spec_list += [sp]
    # spec_list contains all echelle orders from different files and multi-extensions
    specs = collate(spec_list)  # now all in a single XSpectrum1D object

    if wv_array is None:
        # bring them to a unique native wavelength grid using PYPIT
        cat_wave = arco.new_wave_grid(specs.data['wave'], wave_method='velocity')
    else:
        cat_wave = wv_array.to('AA').value
    if rebin is not None:
        rebin = int(rebin)
        cat_wave = cat_wave[::rebin]
    specs = specs.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none',grow_bad_sig=True)

    # estimate weights for coaddition (PYPYT)
    sn2, weights = arco.sn_weight(specs)

    # coaddition
    spec1d = arco.one_d_coadd(specs, weights)
    return spec1d


def load_single_x1d_stis(filename, debug=False):
    """

    Parameters
    ----------
    filename : str
        Filename of the fits x1d STIS file
        Could me multiextension

    Returns
    -------
    spec_list : list of XSpectrum1D objects, one for each echelle order
        of the single STIS x1d file

    """
    # get number of extensions
    head = fits.getheader(filename, ext=0)
    numext = head['NEXTEND']
    spec_list = []  # store XSpectrum1D here.

    for ext in range(1, numext + 1):
        sp = fits.getdata(filename, ext=ext)
        print("Loading echelle orders from file {}, ext={}".format(filename, ext))
        for ii in range(len(sp.SPORDER)):
            # chop pixels at edges of orders (i.e. poor sensitivity)
            fl = sp.FLUX[ii][5:-20]
            er = sp.ERROR[ii][5:-20]
            wv = sp.WAVELENGTH[ii][5:-20]
            spec = XSpectrum1D.from_tuple((wv,fl,er))
            spec_list += [spec]
            if debug:
                pl.plot(sp.WAVELENGTH[ii], sp.FLUX[ii], drawstyle='steps-mid')
                pl.plot(sp.WAVELENGTH[ii], sp.ERROR[ii], ":")
    return spec_list
