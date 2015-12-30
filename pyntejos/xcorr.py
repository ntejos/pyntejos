#Module to deal with the 2D2PCF
import time
import numpy as np
import pylab as pl

def auto_pairs_rt(x, y, z, rbinedges, tbinedges,wrapped=True):
    """ Find the number of pairs in 2d and 3d for galaxies with
    coordinate x, y, z.

    x,y,z: comoving coordinates in Mpc.  x is (roughly) the redshift
           direction.

    rbinedges and tbinedges correspond to the arrays defining the bins
    edges """

    start=time.clock()
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    rbinedges = np.array(rbinedges)
    tbinedges = np.array(tbinedges)
    
    npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)

    for i in xrange(len(x) - 1):
        # radial separation
        if wrapped:
            rsep = np.abs(x[i+1:] - x[i])
        else:
            rsep = x[i+1:] - x[i]
        # transverse separation
        tsep = np.hypot(y[i+1:] - y[i], z[i+1:] - z[i])

        vals, _ = np.histogramdd((rsep, tsep), (rbinedges, tbinedges),
                              range=None, normed=False, weights=None)
        npair_rt += vals

    end=time.clock()
    print '\t Time elapsed = %s seconds.'%(end-start)
    return npair_rt



def cross_pairs_rt(x1, y1, z1, x2, y2, z2, rbinedges, tbinedges,wrapped=True):
    """ Find the number of pairs in 2d and 3d for galaxies with
    coordinate (x1, y1, z1) and (x2,y2,z2).

    x,y,z: comoving coordinates in Mpc.  x is (roughly) the redshift
           direction
    """
    start=time.clock()
    x1 = np.array(x1)
    y1 = np.array(y1)
    z1 = np.array(z1)
    x2 = np.array(x2)
    y2 = np.array(y2)
    z2 = np.array(z2)
    rbinedges = np.array(rbinedges)
    tbinedges = np.array(tbinedges)
    
    npair_rt = np.zeros((len(rbinedges) - 1, len(tbinedges) - 1), float)

    if len(x1)>len(x2):
        auxx=x1
        auxy=y1
        auxz=z1
        x1=x2
        y1=y2
        z1=z2
        x2=auxx
        y2=auxy
        z2=auxz
    
    for i in xrange(len(x1) - 1):
        # radial separation
        if wrapped:
            rsep = np.abs(x1[i] - x2)
        else:
            rsep = x1[i] - x2
        # transverse separation
        tsep = np.hypot(y1[i]-y2, z1[i]-z2)
        
        vals, _ = np.histogramdd((rsep, tsep), (rbinedges, tbinedges),
                              range=None, normed=False, weights=None)
        npair_rt += vals

    end=time.clock()
    print '\t Time elapsed = %s seconds.'%(end-start)
    return npair_rt


def W12(DD,RR,f=None):
    """ It computes,

    (1) the Natural/Peebles estimator: DD/(f*RR) - 1
    (2) its poissonian error
    
    for a given number of pairs. """
    
    Ndd = np.sum(DD) # different than Ngal*(Ngal-1)/2
    Nrr = np.sum(RR) # different than Nrdm*(Nrdm-1)/2

    if f==None:
        f = Ndd/Nrr
    RR = np.where(RR==0,1e-10,RR)
    
    W1 = DD / (f*RR)  - 1.
    W1 = np.where(RR==0,-2,W1)

    t1 = DD / ((f*RR)**2)
    t2 = (DD**2) / (RR**3)
    
    err_W1 = np.sqrt(t1 + t2)
    err_W1 = np.where(err_W1==0,1,err_W1)
    return W1, err_W1



def W3b(DD,RR,DR,RD,Ndd=None,Nrr=None,Ndr=None,Nrd=None):
    """ It computes,

    (1) the Landy&Szalay estimator: (DD - DR - RD + RR) / RR 
    (2) its poissonian error
    
    for a given number of pairs. f1,f2,f3 corresponds to the
    normalization factors such that: W3 = DD/(f1*RR) - DR/(f2*RR) -
    RD/(f3*RR) + 1"""

    if Ndd is None:
        Ndd = np.sum(DD)
    if Nrr is None:
        Nrr = np.sum(RR)
    if Ndr is None:
        Ndr = np.sum(DR)
    if Nrd is None:
        Nrd = np.sum(RD)
    
    #normalize counts
    DD = DD / Ndd 
    RR = RR / Nrr
    DR = DR / Ndr
    RD = RD / Nrd


    #assert np.sum(RR==0)==0, "RR has values equal to zero. Try again."
    RR = np.where(RR==0,1e-10,RR)

    W3 = DD/RR  - DR/RR - RD/RR + 1.
    W3 = np.where(DD==0,-1,W3)

    t1 = DD*Ndd / ((Nrr*RR)**2)
    t2 = DR*Ndr / ((Nrr*RR)**2)
    t3 = RD*Nrd / ((Nrr*RR)**2)
    t4 = ((DD*Ndd - DR*Ndr - RD*Nrd)**2 ) / ((Nrr*RR)**3)
    
    err_W3 = np.sqrt(t1 + t2 + t3 + t4)
    err_W3 = clean_matrix(err_W3,n=0)
    err_W3 = np.where(err_W3<=0,1.9,err_W3)
    
    return W3,err_W3

def W3(DD,RR,DR,RD,Ndd=None,Nrr=None,Ndr=None,Nrd=None):
    """Returns the Landy & Scalay estimator for the cross-correlation
    and its error. The pair counts DD, RR, DR, RD are not required to
    be normalized normalized. Ndd,Nrr,Ndr,Nrd are the normalizations
    factors for each pair count (if None, the normalizations are taken
    from the sum over the whole array."""
    
    if Ndd is None:
        Ndd = np.sum(DD)
    if Nrr is None:
        Nrr = np.sum(RR)
    if Ndr is None:
        Ndr = np.sum(DR)
    if Nrd is None:
        Nrd = np.sum(RD)

    #normalize counts
    DD = DD / Ndd 
    RR = RR / Nrr
    DR = DR / Ndr
    RD = RD / Nrd
    
    RR = np.where(RR==0,1e-20,RR)
    W3 = DD/RR  - DR/RR - RD/RR + 1.
    
    #error from LS93, eqs. (48) (46) (43) and definition of Gp
    #Here Gp = RR when normalized
    err_W3 = np.sqrt( (1 + W3)**2 / (Ndd*RR) )
    #err_W3 = np.sqrt( (1 + W3) / (Ndd*DD) )

    #integral constrain correction, LS93 eq. (48) (46)
    C = np.sum(W3*RR)
    #print 'integral constraint: %s'%C
    W3 = (1 + C)*(1+W3) - 1
    if np.fabs(C)>1:
        print 'integral constraint: %s'%C


    #clean W3 and errors
    #W3     = clean_matrix(W3,n=0)
    W3     = np.where(RR==0,0,W3)
    W3     = np.where(DD==0,-1,W3)
    err_W3 = clean_matrix(err_W3,n=0)
    #err_W3 = np.where(err_W3<=0,1.9,err_W3)
    return W3, err_W3


def W_sig(W,err_W,s=1.):
    
    """ It return the values of W which are s times greater than
    W/err_W, the rest of the values are set to 0."""

    Ws = W
    Ws = np.where(W/err_W < s, -1, Ws)
    return Ws



def clean_matrix(W,n=-99):
    W = np.where(np.isnan(W),n,W)
    W = np.where(np.isinf(W),n,W)
    return W


def plot_NM(tbinedges,rbinedges, thingstoplot,thingsnames,n,m,vmin=None, vmax=None):
    
    """ It plots N rows and M columns of the array of things,
    thingstoplot. The name of each plot is given by thingsnames
    array"""
    
    pl.clf()
    pl.subplots_adjust(wspace=0.6,hspace=0.4)
    for i,t,name in zip(xrange(len(thingstoplot)),thingstoplot,thingsnames):
        
        ax = pl.subplot(n,m,i+1)
        c  = pl.pcolormesh(tbinedges,rbinedges,t,vmin=vmin,vmax=vmax)#,cmap=pl.cm.gray)
        cb = pl.colorbar(c,ax=ax,fraction=0.1)
        #ax.axis('tight')
        ax.axis('scaled')
        ax.set_title(name,fontsize='large')
        ax.set_xlabel(r'$r_{perp}$ (Mpc/h)')
        ax.set_ylabel(r'$r_{LOS}$ (Mpc/h)')
