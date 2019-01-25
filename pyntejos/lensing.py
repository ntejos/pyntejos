import numpy as np
from mpdaf.obj import Image

"""Utilities for lensing de-lensing"""

def ima2abs(ra_deg, dec_deg, arc_name='PSZ1GA311_G1'):
    """ Return de-lensed coordinate for a given arc
    using the corresponding deflection matrix.
    (c) Sebita 2019"""

    # number of positions
    N = 1
    # coordinates in image plane (deg)
    # xdeg=np.array([237.5082875])
    # ydeg=np.array([-78.1846194])
    xdeg=np.array([ra_deg])
    ydeg=np.array([dec_deg])

    # lens stuff, deflexion matrices, factor to normalize
    # lensfactor = d_ls/ds(z=0.73) / d_ls/ds(z=2.37)
    if arc_name == 'PSZ1GA311_G1':
        lensfactor = 0.4958 # depends on redshift
        alpha_scale = 8.33333333333333E-06  # scale in deflection matrix
        prefix = '/media/ntejos/disk1/projects/arc_tomo/PSZ1GA311/data/delensing/'
        alpha_x = lensfactor*Image(prefix + 'dplx1_Fel_wcs.fits')
        alpha_y = lensfactor*Image(prefix + 'dply1_Fel_wcs.fits')
    else:
        raise ValueError('Not implemented for that arc_name.')

    # coordinates in pixels of alpha matrix
    xpix,ypix = [0]*N,[0]*N
    for i in range(1):
        xpix[i]=alpha_x.wcs.sky2pix((ydeg[i],xdeg[i]),nearest=True)[0] [1]
        ypix[i]=alpha_x.wcs.sky2pix((ydeg[i],xdeg[i]),nearest=True)[0] [0]

    # map to absorber plane (using the alpha matrix)
    xdeg_l = xdeg + alpha_x[ypix,xpix].data/3600./np.cos(np.deg2rad(ydeg))
    ydeg_l = ydeg - alpha_y[ypix,xpix].data/3600.

    return (xdeg_l[i],ydeg_l[i])