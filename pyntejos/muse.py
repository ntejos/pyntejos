"""This module is meant for analysis of VLT/MUSE data"""

import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.nddata import NDData, StdDevUncertainty
import astropy.units as u

filename = "/media/ntejos/disk2/data/MUSE/094.A-0575C/reduced/DATACUBE_FINALSKY_ALL_1.fits"

class IFUCube(object):
    """A Class for dealing with IFU datacubes. Currently
    implemented for VLT/MUSE.

    Parameters
    ----------
    filename : str
        Name of the file containing the IFU datacube
    instrument : str
        Name of the instrument from where the IFU datacube come from.
        e.g. ('MUSE')

    """
    
    def __init__(self,filename,instrument):
        good_inst = ('MUSE')
        if instrument not in good_inst:
            raise NotImplementedError('Not ready for {} data yet'.format(instrument))
        if instrument == 'MUSE':
            self.fl,self.er,self.wa = self.read_muse_data(filename)

    def read_muse_data(self,filename):
        """Reads MUSE datacubes

        Returns
        -------
        fl, er : astropy.NDData arrays
        wa : np.array 

        """
        fl = fits.getdata(filename,ext=1)
        er = fits.getdata(filename,ext=2) #variance

        h_fl = fits.getheader(filename,ext=1)
        h_er = fits.getheader(filename,ext=2)

        #make sure fl and er are correct extensions
        if h_fl['extname'] != 'DATA':
            raise ValueError('Data array is not in extension 1 of {}!'.format(filename))
        if h_er['extname'] != 'STAT':
            raise ValueError('Uncertainty array is not in extension 2 of {}!'.format(filename))
        #read units
        if h_fl['bunit']=='10**(-20)*erg/s/cm**2/Angstrom':
            fl_unit = 10**(-20.) * u.erg / u.s / (u.cm**2) / u.AA
        else:
            raise NotImplementedError('The units of the datacube were not recognized')
        #read units
        if h_er['bunit']=='(10**(-20)*erg/s/cm**2/Angstrom)**2':
            er_unit = 10**(-20.) * u.erg / u.s / (u.cm**2) / u.AA
            #from variance to stddev
            er = np.sqrt(er)
        else:
            raise NotImplementedError('The units of the datacube were not recognized')

        #relevant quantities
        crval1  = h_fl['crval1']
        crpix1  = h_fl['crpix1']
        cd1_1   = h_fl['cd1_1']
        cd1_2   = h_fl['cd1_2']
        naxis1  = h_fl['naxis1']
        crval2  = h_fl['crval2']
        crpix2  = h_fl['crpix2']
        cd2_1   = h_fl['cd2_1']
        cd2_2   = h_fl['cd2_2']
        naxis2  = h_fl['naxis2']
        crval3  = h_fl['crval3']
        crpix3  = h_fl['crpix3']
        cd3_3   = h_fl['cd3_3']
        naxis3  = h_fl['naxis3']

        #wavelength array
        wa = (np.arange(naxis3)+1-crpix3)*cd3_3+crval3

        return fl, er, wa

    def add_wa_mask(self,wa_min,wa_max):
        cond = (self.wa>wa_min) & (self.wa<wa_max)
        self.wa_mask = self.wa_mask | cond

    def reset_wa_mask(self):
        cond = self.wa < np.min(self.wa)
        self.wa_mask = cond

    def reset_xy_mask(self):
        pass