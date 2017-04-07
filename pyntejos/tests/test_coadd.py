from pyntejos import coadd
from matplotlib import pyplot as plt
import glob
from linetools.spectra.io import readspec


def test_coadd_stis_from_x1dfiles(plot=False):
    spec_old = readspec("./data/m5zng1_stis/m5zng1_stis.fits")

    # co-add independently of old
    filenames = glob.glob("./data/m5zng1_stis/o6n*_x1d.fits")
    spec_new = coadd.coadd_stis_from_x1dfiles(filenames, wv_array=None)

    spec_new_rebinned = spec_new.rebin(spec_old.wavelength, do_sig=True, grow_bad_sig=True)
    if plot:
        plt.plot(spec_old.wavelength, spec_old.flux, drawstyle='steps-mid', color='k')
        plt.plot(spec_old.wavelength, spec_old.sig, drawstyle='steps-mid', color='g')
        plt.plot(spec_new.wavelength, spec_new.flux, drawstyle='steps-mid', color='b')
        plt.plot(spec_new.wavelength, spec_new.sig, drawstyle='steps-mid', color='y')

    spec_new.write_to_fits("new.fits")
