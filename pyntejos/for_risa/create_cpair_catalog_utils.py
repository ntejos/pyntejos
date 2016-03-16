import numpy as np
import astropy.units as u
from astropy.table import Table
import pyntejos.coord as coord
from pyntejos.utils import give_dv

"""Utilities for create_cpair_catalog.py"""


def create_pairs(ra,dec,z,z_err,cosmo,max_tsep=20, z_err_max1=0.0005,z_err_max2=0.03,objid=None,dvmax=3000,verbose=False):
    """From a set of astronomical objects with ra, dec, z, z_err, it
    creates a 'pair catalog' with the following elements:
    objid1,objid2,ra1,ra2,dec1,dec2,z1,z2,redshift,delta_z,sep_mpc,sep_deg,sen(a),cos(a),quality
    
    Inputs:
    ra, dec, z, z_err: celestial coordinates of objects and redshifts (ra and dec in degrees!)
    z_err_max1: maximum redshift error for at least one member of the pairs
    z_err_max2: maximum allowed redshift uncertanty for both members
    max_tsep: maximum transverse separation allowed for the members of the pairs
    objid: an array with the obj unique ids (optional). If None, objid1
           and objid2 are given by the position in the given arrays.
    dvmax: maximum velocity difference between pairs allowed in km/s, using nominal redshifts.
    """
    
    ra    = np.array(ra)
    dec   = np.array(dec)
    z     = np.array(z)
    z_err = np.array(z_err)
    assert (len(ra) == len(dec)) and (len(ra)==len(z)) and (len(ra)==len(z_err)), \
        'ra,dec,z,z_err arrays are not of the same size!'
    
    pairs = [] 
    
    ra_rad  = ra  * np.pi/180. #ra in radians
    dec_rad = dec * np.pi/180. #dec in radians
    
    for i in np.arange(len(ra)): 
        if verbose and (i%1000 == 0):
            print '{}/{}'.format(i,len(ra))
        #maximum radius to perform the search in radians
        max_radius = max_tsep * u.Mpc / cosmo.comoving_distance(z[i])
        #angular separation between (ra[i],dec[i]) and all other objects in radians
        obj_obj_sep = coord.ang_sep(ra[i],dec[i],ra,dec) * coord.RAD_PER_DEG
        #condition
        cond = (obj_obj_sep < max_radius)
        
        obj_good = np.where(cond)[0] # has the original position of
                                     # objects satisfying the condition
        for j in obj_good:
            quality = 0  #for quality, switch to 1 when both have spec-z
            if j > i:
                #the redshifts of potential members of pairs
                zi = z[i]
                zj = z[j]
                zi_err = z_err[i]
                zj_err = z_err[j]
                
                #redshift condition here
                z_diff = np.fabs(zi-zj)
                err_z  = np.sqrt(zi_err**2 + zj_err**2) #combined z-error
                condz  = (np.min(zi_err,zj_err) < z_err_max1) & (err_z < z_err_max2)
                
                if err_z < z_err_max1:#if both members have good
                                      #redshifts, only use the velocity condition below
                    condz = True
                
                #velocity condition
                condv = give_dv(np.fabs(zi-zj)+np.mean([zi,zj]),np.mean([zi,zj])) < dvmax
                
                if (condz)&(condv):
                    #use the mean redshift as redshift
                    redshift = (zi/zi_err + zj/zj_err) / (1./zi_err + 1./zj_err)
                    redshift_err = np.sqrt(2*(zi_err**2))
                    
                    if (zi_err < z_err_max1) and (zj_err < z_err_max1): #if both have good measurement, set quality=1
                        quality = 1 #if both members have spec-z set good quality
                    
                    #recompute ang sep in a one-by-obe basis, in radians
                    sep_deg = coord.ang_sep(ra[i], dec[i], ra[j], dec[j])
                    sep_rad = sep_deg * coord.RAD_PER_DEG
                    sep_mpc = sep_rad * (cosmo.comoving_distance(redshift)).value #angular separation in Mpc (co-moving)
                    cond_ang_sep = (sep_mpc <= max_tsep)
                    
                    if (cond_ang_sep):
                        
                        #print i,j,redshift

                        ra1 = ra[i]
                        ra2 = ra[j]
                        dec1 = dec[i]
                        dec2 = dec[j]
                        
                        #rotation angle
                        if (ra1 < ra2):
                            sena = (dec2 - dec1) * coord.RAD_PER_DEG / sep_rad
                            cosa = (ra2 - ra1) * coord.RAD_PER_DEG  / sep_rad
                        else:
                            sena = (dec1 - dec2) * coord.RAD_PER_DEG / sep_rad
                            cosa = (ra1 - ra2) * coord.RAD_PER_DEG  / sep_rad
                        
                        #keep the pair
                        if objid is not None:
                            pair_info = (objid[i],objid[j],ra[i],ra[j],dec[i],dec[j],z[i],z[j],redshift,redshift_err,sep_mpc,sep_deg,sena,cosa,np.fabs(z[i]-z[j]),quality)
                        else:
                            pair_info = (i,j,ra[i],ra[j],dec[i],dec[j],z[i],z[j],redshift,redshift_err,sep_mpc,sep_deg,sena,cosa,np.fabs(z[i]-z[j]),quality)
                        pairs = pairs + [pair_info]
                        
                        #keep track of progress
                        #if i/100 == i/100.:
                        #    print i,j, '\t',sep_mpc, '\t',np.fabs(z[i]-z[j]), redshift, redshift_err,quality
    #re-structure the pair array
    if len(pairs)>0:
        pairs = np.array(pairs)
        pairs = pairs.view(np.recarray)
        pairs = np.rec.fromrecords(pairs,names='objid1,objid2,ra1,ra2,dec1,dec2,z1,z2,redshift,redshift_err,sep_mpc,sep_deg,sena,cosa,sep_z,quality')
        #Table version
        pairs = Table(pairs)
        pairs.sort(['redshift'])
        pairs['objid1'] = [int(objid) for objid in pairs['objid1']] #make sure there are ints
        pairs['objid2'] = [int(objid) for objid in pairs['objid2']] #make sure these are ints
    else:
        pairs = Table()
    return pairs
 