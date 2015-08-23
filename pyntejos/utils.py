from astropy.cosmology import WMAP7 as cosmo
from astropy.constants import c as C
import numpy as np

"""Module for utils treating catalogs"""

def compare_z(z1,z2,dv):
    """Return true if the difference between z1 and z2 is within dv at
    the mean redshift. Otherwise return False. dv has to be much
    smaller and speed of light as the function does not account for
    relativity."""
    z1 = np.array(z1)
    z2 = np.array(z2)
    dc = np.array(dv)

    dz = np.fabs(z1 - z2)
    z  = np.mean([z1,z2])
    if dz / (1. + z) < dv / C.to('km/s').value:
        return True
    else:
        return False
    

def group_z(z,dv=1000):
    """Group redshifts within dv (km/s) at the redshift of the
    group. Returns an array of id_membership, of the same dimension
    than z"""
    z_original = np.array(z)
    z = np.array(z)
    z.sort() 
    
    ids = np.zeros(len(z)) # This is the output, here we store the id_membership of each z
    ids_aux = np.zeros(len(z)) # This is the same, but matching a sorted z array
                               
    q = 0 #counter for groups
    ids[0] = q 
    for i in np.arange(1,len(z)):
        if compare_z(z[i],z[i-1],dv):
            ids_aux[i] = q
        else:
            q += 1
            ids_aux[i] = q

    #remap ids_aux to ids
    for i in np.arange(len(z_original)):
        cond = z_original[i] == z
        if np.sum(cond)>1: #in case there are 2 or more with same z
            ids[i] = ids_aux[cond][0]
        else:
            ids[i] = ids_aux[cond]

    return ids.astype(int)

def give_dv(z,zmean,rel=True):
    """Gives velocity difference in km/s between z and zmean, using
    relativistic approximation for a locally flat space-time."""
    z = np.array(z)
    zmean = np.array(zmean)
    
    if rel:
        dv = ((1+z)**2-(1+zmean)**2) / ((1+z)**2 + (1+zmean)**2) * C.to('km/s').value 
    else:
        dv = (z - zmean) / (1. + zmean) * C.to('km/s').value
    
    return dv

def give_dz(dv,zmean,rel=True):
    """Gives redshift difference of dv in km/s at the redshift zmean,
    using relativistic approximation for a locally flat space-time."""
    
    dv = np.array(dv) 
    zmean = np.array(zmean)
    
    if rel:
        beta = dv/C.to('km/s').value
        aux = np.sqrt((1.+beta)/(1.-beta))
        dz = (1. + zmean) * (aux - 1.)
    else:
        dz = dv * (1. + zmean) / C.to('km/s').value
    return dz


def poisson_err(n):
    """Gets poissonian error analytical approximation from equations of
    Gherels 1986. Returns the upper and lower uncertainties"""
    
    n = np.array(n)
    
    errp = np.sqrt(n+0.75)+1.
    errp = np.where(errp==errp,errp,errp)
    errm = np.where(n>0.25,np.sqrt(n-0.25),0)
    
    return errp,errm




def find_edges(a):
    """Assume a is 1-D array of values, where 0 mean masked out. This
    function will provide the indices of lower and upper edges of
    chunks having values different than 0. Useful for spectra
    analyses"""
    a = np.array(a)
    
    #append 0 before and after the original array to an auxiliary array
    #in this way:
    #lower case limits are always 0,something
    #upper case limits are always something,0
    a_aux = np.append(0,a)
    a_aux = np.append(a_aux,0)
    
    lower = []
    upper = []
    
    
    for i in range(1,len(a_aux)-1): #only for indices with original info
        if (a_aux[i]!=0) and (a_aux[i-1]==0): #lower case
            lower += [i]
        if (a_aux[i]!=0) and (a_aux[i+1]==0): #upper case
            upper += [i]
    lower = np.array(lower)
    upper = np.array(upper)
    #lower and upper have indices of a_aux
    #we substract one to get the indices in the original array a
    lower = lower - 1
    upper = upper - 1
    assert len(lower)==len(upper),'Something is wrong with find_edges function. Debug code!'
    return lower.astype(int),upper.astype(int)
    
def is_local_minima(a):
    """For a given array a, it returns true for local minima"""
    a = np.array(a)
    mask = []
    for i in range(1,len(a)-1):
        cond = (a[i] < a[i-1]) and (a[i] < a[i+1])
        if cond:
            mask += [1]
        else:
            mask += [0]
    mask = np.array(mask)
    mask = np.append(0,mask)
    mask = np.append(mask,0)
    return mask == 1

def is_local_maxima(a):
    """For a given array a, returns true for local maxima"""
    a = np.array(a)
    mask = []
    for i in range(1,len(a)-1):
        cond = (a[i] > a[i-1]) and (a[i] > a[i+1])
        if cond:
            mask += [1]
        else:
            mask += [0]
    mask = np.array(mask)
    mask = np.append(0,mask)
    mask = np.append(mask,0)
    return mask == 1

def associate_redshifts(z1,z2,dv):
    """Returns an array of same lenght as z1, with a 1 if z1 redshift
    lie within dv from any of the redshifts given by the array z2;
    otherwise it has a value of 0"""
    
    z1 = np.array(z1)
    z2 = np.array(z2)
    
    association = np.zeros(len(z1))
    for z in z2:
        dv_aux = np.fabs(give_dv(z1,z))
        association = np.where(dv_aux<dv,1,association)
    return association

def get_dv_closest_z(z1,z2,give_ids=False):
    """Returns an array of same lenght as z1, with the velocity difference
    (in km/s) associates to the closest redshift in z2, at restframe
    given by z2. (Using relativistic approximation for flat-space time;
    see give_dv() function)

    If give_ids is True , it also returns the indices of z2 where the
    difference is minimum.

    """

    z1 = np.array(z1)
    z2 = np.array(z2)
    
    dv = []
    inds = []
    for z in z1:
        dv_aux = give_dv(z,z2)
        #find minimum difference
        cond = np.fabs(dv_aux)==np.min(np.fabs(dv_aux))
        ind = np.where(cond)[0][0]
        #append dv
        dv += [dv_aux[ind]]
        inds += [ind]
    dv = np.array(dv)
    if give_ids:
        return dv, inds
    else:
        return dv

def clean_array(x, value=0):
    """Get rid of nan and inf in the array x, and replace them with the
    given value."""
    x = np.array(x)
    cond = ( np.isnan(x) ) | ( np.isinf(x) ) 
    x = np.where(cond,value,x)
    return x


def get_original_indices(original,new):
    """For a pair of arrays containing the same information but sorted in a
    different way, this function provides the indices associated to the
    original array that make up the new array so
    original[indices]=new.

    Add check to make sore arrays contain exact same information

    """
    original = np.array(original)
    new = np.array(new)
    
    #use np.argsort()
    inds_orig_sorted = np.argsort(original)
    inds_new_sorted = np.argsort(new)

    #find indices such that original[indices] = new
    indices_aux = np.argsort(inds_new_sorted)
    indices = inds_orig_sorted[indices_aux]

    return indices
