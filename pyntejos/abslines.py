import numpy as np
import matplotlib.pyplot as pl
from astropy.constants import c as C
from scipy.ndimage import uniform_filter as uf
from astro.sampledist import RanDist
from pyntejos.utils import find_edges,clean_array
from astropy import constants as const
import astropy.units as u
from barak.absorb import readatom
from linetools.analysis import absline as ltaa

#Constant
e2_me_c = ((const.e.esu)**2/(const.c.to('cm/s')*const.m_e.to('g'))).value # from Draine (eq. 9.8 and 9.9)

#read atomic data
atomdat = readatom('/home/ntejos/software/vpfit10/atom.dat')


def get_tau0_peak(wa0,fosc,logN,b):
    """Get the value of the optical depth at the line center,
    tau0. From Draine 2011 (see Chapter 9). It neglects estimulated
    emission which is fine for IGM or ISM except for radio-frecuency
    transitions.
    
    Inputs
    ------
    wa0:   rest-frame wavelenght of the transition in A
    fosc:  oscillator strenght for the transition
    logN:  log10 of the column density in cm^-2
    b:     Doppler parameter in km/s
    
    Returns
    -------
    tau0:  optical depth at the line center
    """
    # give units
    wa0 = wa0 * u.AA
    b = b * u.km/u.s
    return ltaa.get_tau0(wa0, fosc, logN, b)

def logN_b_to_Wr(logN,b,wa0,fosc,gamma):
    """I will use the approximation given by Draine book (eq. 9.27),
    whic comes from atomic physics considerations + Rodgers & Williams
    1974 (I couldn't find the reference though)
    
    Inputs
    ------
    logN:  log10 of the column density in cm^-2.
    b:     Doppler parameter in km/s.
    wa0:   rest-frame wavelenght of the transition in A.
    fosc:  oscillator strenght for the transition.
    gamma: Gamma parameter for the transition in s^-1.
    
    Returns
    -------
    Wr:  rest-frame equivalent width in A
    
    """
    # give units
    wa0 = wa0 * u.AA
    N = 10**logN * (1/u.cm*u.cm)
    b = b * u.km/u.s
    gamma = gamma * (1/u.s)

    # call linetools function
    return ltaa.Wr_from_N_b(N, b, wa0, fosc, gamma)


def logN_b_to_Wr_ion(logN,b,ion='HI',wa0=1215.67):
    """I will use the approximation given by Draine book (eq. 9.27),
    whic comes from atomic physics considerations + Rodgers & Williams
    1974 (I couldn't find the reference though)
    
    Inputs
    ------
    logN:  log10 of the column density in cm^-2.
    b:     Doppler parameter in km/s.
    ion:   Name of ion.
    wa0:   rest-frame wavelenght in A.

    Returns
    -------
    Wr:  rest-frame equivalent width in A
    """
    
    logN = np.array(logN)
    b    = np.array(b)
    wa0  = np.array(wa0)   

    #read atomic data
    try:
        atom_info = atomdat[ion]
    except:
        assert False,\
            'The ion name provided does not match any of the ones in the current database.'
    
    #find the transition parameters by matching the closest wa0
    atom_wa = np.array(atom_info['wa'])
    ind = np.where(np.fabs(wa0-atom_wa)==np.min(np.fabs(wa0-atom_wa)))[0][0]
    wa0   = atom_info[ind]['wa']  #redifine wa0
    fosc  = atom_info[ind]['osc'] #oscillator strenght
    gamma = atom_info[ind]['gam'] #Gamma parameter
    
    #find rest-frame equivalent width
    Wr = logN_b_to_Wr(logN,b,wa0,fosc,gamma) #in Angstroms
    return Wr

def compute_Wmin(wa,fl,er,sl=3.,R=20000,FWHM=10,wa0=1215.67,mask_Gal= True, fl_th = 0):
    """For a given spectrum and transition, it computes the minimun
    rest-frame equivalent width for that transition to be observed. It
    return a tuple of redshift and Wmin (z,Wmin).

    We calculate the Wmin = sl * wa * / (1+z) / R / (fl/er), where z =
    wa/wa0 - 1 (wa0 is the rest frame wavelenght of the transition) and
    R is the resolution of the spectrograph. We then smooth Wmin with a
    boxcar along FWHM pixels.

    Inputs
    ------
    wa : array
      observed wavelenght array
    fl : array
      observed flux array
    er : array
      uncertainty in the observed flux array
    sl : float
      significance level for the detection of line
    R : float
      resolution of the spectrograph
    FWHM : float
      Size (in pixels) of the smoothing window to calculate Wmin; the
      smoothing kernel is boxcar
    wa0 : float
      rest-frame wavelenght of the transition of interest
    mask_Gal : boolean
      If true, mask out Galactic features within 200 km/s
    fl_th : float
      Flux threshold, mask out spectral regions where fl<fl_th
      
    Returns
    -------
    z : array
      redshift, array of same size as wa with wa/wa0 - 1
    Wmin : array
      Minimum rest-frame equivalent width in A, array of same size as wa.
    """
   
    #masked regions (potential Galactic absorption)
    if mask_Gal:
        galactic=np.array([1334.5323 , # CII
                           1238.821  , # NV
                           1242.804  , # NV
                           1302.1685 , # OI
                           1304.8576 , # OI*
                           1306.0286 , # OI**
                           1304.3702 , # SiII
                           1260.4221 , # SiII
                           1334.8132 , # PIII
                           1259.519  , # SII
                           1253.811  , # SII
                           1250.584  , # SII
                           1260.533])  # FeII
        zgal=galactic/wa0 - 1. #Position of Galactic features.
        dzgal=(zgal+1.) * 200 / C.to('km/s').value # this is to mask
                                                   # regions within 200
                                                   # km/s of Galactic
                                                   # absorption.
    
    wa = np.array(wa)
    fl = np.array(fl)
    er = np.array(er)
    
    z = wa/wa0 - 1.   # spectrum in z coordinates
    

    SN = (fl/er)
    Wmin = sl * wa0 / R / SN #equivalent to sl*wa / (1. + z) / R / (S/N)
    Wmin = np.where(Wmin<=0,1e10,Wmin)
    Wmin = np.where(np.isnan(Wmin),1e10,Wmin)
    Wmin = np.where(np.isinf(Wmin),1e10,Wmin)
    Wmin = uf(Wmin.astype(float),FWHM) # smoothed version (uniform prefered over gaussian) 
    Wmin = np.where(fl<fl_th,1e10,Wmin) #mask out spectral regions where fl<tl_th
    if mask_Gal:
        for zi,dzi in zip(zgal,dzgal): #mask out Galactic regions
            cond = (z>zi-dzi)&(z<zi+dzi)
            Wmin = np.where(cond,1e10,Wmin)
    return z, Wmin

def compute_local_SN(wa0,wa,fl,er, npix=100, fl_th = 0.9, debug=0):
    """It computes an average local signal to noise ratio (S/N) around
    the position wa0, over npix pixels, in a given spectrum.
    
    Inputs
    ------
    wa0 : float
      central wavelenght to compute the S/N ratio.
    wa : array
      array of wavelenght values for the spectrum.
    fl : array
      array of fluxes values for the spectrum.
    er : array
      array for uncertainties in the flux spectrum.
    npix : float (forced to be int)
      number of pixels to compute the S/N estimation.
    fl_th : float or array of same dimension than fl
      flux threshold. If fl < fl_th those pixels are masked out of the
      calculation
    
    Returns
    -------
    sn : float
      average signal to noise under the above constraints.
    std : float
      standard deviation of the sn
    """

    wa = np.array(wa)
    fl = np.array(fl)
    er = np.array(er)
    
    fl = clean_array(fl)
    er = clean_array(er)
    
    assert (len(wa)==len(fl)) and (len(wa)==len(er)),\
        'wa,fl and er do not have same dimensions!'
    
    npix = int(npix)
    
    #if fl_th is float , then use it in units of fl
    if type(fl_th) is float:
        fl_th = fl * fl_th
    
    #assert whether wa0 is in wa
    cond = (wa0>=wa[0])&(wa0<=wa[-1])
    assert cond, '[compute_local_SN] wa0 is out of range'
    
    #find pixel index for wa0
    cond = np.fabs(wa0-wa) == np.min(np.fabs(wa0-wa))
    ind = np.where(cond)[0][0]
    
    #define the wa window with at least npix pixels
    i = 0
    while True:
        imin = int(ind-0.5*npix-i)
        imax = int(ind+0.5*npix+i)
        if imax > len(wa) -1:
            imax = len(wa)-1
        if imin < 0:
            imin = 0
        if debug:
            print imin,imax   
        
        cond = (wa >= wa[imin]) & (wa < wa[imax]) & (fl > fl_th)
        if np.sum(cond) >= npix:
            if debug:
                print np.sum(cond)
            break
        else:
            i += 1
        assert i < len(wa), \
            'The spectrum does not satisfy the given conditions, try other values!'
    
    #define the chunk of spectra to look at
    fl_aux = fl[cond]
    er_aux = er[cond]
    sn = fl_aux/er_aux
    sn = clean_array(sn,value=np.median(sn))
    return np.mean(sn), np.std(sn)    



def random_abs(zobs,Wobs,Nrand,wa,fl,er,sl=3.,R=20000,FWHM=10.,wa0=1215.67,zmax=None,zmin=None):
    """From a list of observed redshifts (zobs) and observed rest-frame
    equivalent widths (Wobs), it creates a random catalog.  For a given
    real absorber with Wobs it creates Nrand random ones at the new
    zrand, defined by where the line could have been observed. It
    returns those zrand.
    
    Inputs
    ------
    zobs:    observed redshifts.
    Wobs:    observed rest-frame equivalent widths (in A).
    Nrand:   number of randWr = logN_b_to_Wr()om lines per real one generated (integer).
    wa:      numpy array of wavelenght covered by the spectrum.
    fl:      numpy array of normalized flux.
    er:      numpy array of error in the normalized flux of the spectrum for 
             a given wavelenght.
    sl:      significance level for the detection of the absorption line.
    R:       resolution of the spectrograph, assumed constant
    FWHM:    Full-width at half maximum in pixels (assumed constant). This 
             parameter defines the smoothing scale for Wmin. 
    wa0:     rest-frame wavelenght (in A).
    zmax:    if given, is the maximum redshift allowed for the randoms
    zmin:    if given, is the minimum redshift allowed for the randoms
    
    We first compute Wmin (see abslines.Wmin). We compute the redshifts
    where Wobs could have been observed according to the given Wmin,
    and place Nrand new absorbers with the same properties as the given
    one accordingly.
    """
    zobs = np.array(zobs)
    Wobs = np.array(Wobs)
    wa = np.array(wa)
    fl = np.array(fl)
    er = np.array(er)

    assert len(zobs)==len(Wobs), \
        'zobs and Wobs do not have same dimensions!'
    assert (len(wa)==len(fl)) and (len(wa)==len(er)),\
        'wa,fl and er do not have same dimensions!'

    Nrand   = int(Nrand)
    zrand   = zobs.repeat(Nrand)
    
    z, Wmin = compute_Wmin(wa,fl,er,sl=sl,R=R,FWHM=FWHM,wa0=wa0)
    
    if zmax is None:
        zmax = 1000. ## Here, sort out zlims
    if zmin is None:
        zmin = 0.

    for i in xrange(len(Wobs)):
        
        zgood  = (Wobs[i] > Wmin) & (z >= zmin) & ( z < zmax)
        
        if np.sum(zgood)==0: #This is not good so plot situation.
            pl.plot(z,Wmin,drawstyle='steps-mid')
            pl.axis([z[0],z[-1],0,0.1])
            pl.show()
            print Wmin
        assert np.sum(zgood)>0, \
            'There are not regions in the spectrum with Wmin<{} A. The minimum is {} A. Adjust significance.'.format(Wobs,np.min(Wmin))
        
        rand_z = RanDist(z, zgood*1.)
        aux    = rand_z.random(Nrand)
        zrand[i*Nrand:(i+1)*Nrand] = aux

    return zrand



def compute_redshift_path_per_line(Wobs, wa, fl,er,sl=3.,R=20000,FWHM=10,wa0=1215.67,fl_th=0,zmin=None,zmax=None):
    """Compute the total redshift path for a given line with Wobs
    (could be a np.array), in a given spectrum fl (with error er)."""
    
    #define N, number of lines
    try:
        N = len(Wobs)
    except:
        N = 1

    #redshift path to return
    Dz = []

    #compute Wmin for given spectrum, and dz
    z, Wmin = compute_Wmin(wa,fl,er,sl=sl,R=R,FWHM=FWHM,wa0=wa0,fl_th=fl_th)
    
    ## Here, sort out zlims
    if zmax is None:
        zmax = 1000. 
    if zmin is None:
        zmin = 0.
    
    if N==1:
        #Obtain zgood (1 or 0 depending if Wobs>Wmin)
        zgood = (Wobs >= Wmin) & (z >= zmin) & ( z < zmax)
        zgood = zgood*1
        if np.sum(zgood)==0:
            Dz = 0
            
        #find edges where zgood=1
        lower_edges, upper_edges = find_edges(zgood)
        
        #Final redshift path is given by the difference between edges
        Dz = np.sum(z[upper_edges]-z[lower_edges])
    else:
        for i in xrange(N):
            #Obtain zgood (1 or 0 depending if Wobs>Wmin)
            zgood = (Wobs[i] >= Wmin) & (z >= zmin) & ( z < zmax)
            zgood = zgood*1
            if np.sum(zgood)==0:
                Dz += [0]
                continue
                
            #find edges where zgood=1
            lower_edges, upper_edges = find_edges(zgood)
        
            #Final redshift path is given by the difference between edges
            Dz_aux = np.sum(z[upper_edges]-z[lower_edges])
            Dz += [Dz_aux]
        Dz = np.array(Dz)
            
    return Dz

def compute_absorption_distance_per_line(cosmo, Wobs, wa, fl,er,sl=3.,R=20000,FWHM=10,wa0=1215.67,fl_th=0,zmin=None,zmax=None):
    """Compute the total absorption distance for a given line with Wobs
    (could be a np.array), in a given spectrum fl (with error er). A
    given cosmology has to be given, assuming form from
    astropy.cosmology"""
    
    #define N
    try:
        N = len(Wobs)
    except:
        N = 1
        
    #absorption distance to return
    DX = []

    #compute Wmin for given spectrum, and dz
    z, Wmin = compute_Wmin(wa,fl,er,sl=sl,R=R,FWHM=FWHM,wa0=wa0,fl_th=fl_th)
    
    if zmax is None:
        zmax = 1000. ## Here, sort out zlims
    if zmin is None:
        zmin = 0.
    
    if N==1:
        #Obtain zgood (1 or 0 depending if Wobs>Wmin)
        zgood = (Wobs >= Wmin) & (z >= zmin) & ( z < zmax)
        zgood = zgood*1
        if np.sum(zgood)==0:
            DX = 0
            
        #find edges where zgood==1
        lower_edges, upper_edges = find_edges(zgood)
        
        #Final absorption distance is given by the redshift difference
        #between edges (absorption distance continuum function)
        DX = np.sum(cosmo.absorption_distance(z[upper_edges])-cosmo.absorption_distance(z[lower_edges]))
    else:
        for i in xrange(len(Wobs)):
            #Obtain zgood (1 or 0 depending if Wobs>Wmin)
            zgood = (Wobs[i] >= Wmin) & (z >= zmin) & ( z < zmax)
            zgood = zgood*1
            if np.sum(zgood)==0:
                DX += [0]
                continue
        
            #find edges where zgood=1
            lower_edges, upper_edges = find_edges(zgood)
        
            #Final absorption distance is given by the redshift difference
            #between edges (absorption distance continuum function)
            DX_aux = np.sum(cosmo.absorption_distance(z[upper_edges])-cosmo.absorption_distance(z[lower_edges]))
            DX += [DX_aux]
        DX = np.array(DX)
    return DX
