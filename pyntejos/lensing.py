import os
import numpy as np
from mpdaf.obj import Image

"""Utilities for lensing de-lensing"""

def ima2abs(ra_deg, dec_deg, arc_name='PSZ1GA311_G1'):
    """ Return de-lensed coordinate for a given arc
    using the corresponding deflection matrix.
    (c) Sebita 2019

    Parameters
    ----------

    ra_deg
    dec_deg

    arc_name : str
        currently available ['PSZ1GA311_G1', 'SGASJ1226_G1']

    """

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
        # alpha_scale = 8.33333333333333E-06  # scale in deflection matrix
        path = os.getcwd()
        if path.startswith('/disk03/ntejos/'):
            prefix = '/disk03/ntejos/projects/arc_tomo/PSZ1GA311/data/delensing/'
        else:
            prefix = '/media/ntejos/disk1/projects/arc_tomo/PSZ1GA311/data/delensing/'
        alpha_x = lensfactor*Image(prefix + 'dplx1_Fel_wcs.fits')
        alpha_y = lensfactor*Image(prefix + 'dply1_Fel_wcs.fits')

    elif arc_name == 'SGASJ1226_G1':
        path = os.getcwd()
        if path.startswith('/disk03/ntejos/'):
            prefix = '/disk03/ntejos/projects/arc_tomo/SGASJ1226/data/delensing/'
        else:
            prefix = '/media/ntejos/disk1/projects/arc_tomo/SGASJ1226/data/delensing/'
        lensfactor = 1.  # it is supposed to be included in the matrices in this case
        alpha_x = lensfactor * Image(prefix + 'dplx_0.77_abs.fits')
        alpha_y = lensfactor * Image(prefix + 'dply_0.77_abs.fits')

    else:
        raise ValueError('Not implemented for that arc_name.')


    # coordinates in pixels of alpha matrix
    xpix,ypix = [0]*N,[0]*N
    for i in range(N):
        xpix[i]=alpha_x.wcs.sky2pix((ydeg[i],xdeg[i]),nearest=True)[0] [1]
        ypix[i]=alpha_y.wcs.sky2pix((ydeg[i],xdeg[i]),nearest=True)[0] [0]
        # print("Pixel (x,y)=({},{})".format(xpix[i],ypix[i]))
    # import pdb;    pdb.set_trace()
    # map to absorber plane (using the alpha matrix)
    xdeg_l = xdeg + alpha_x[ypix,xpix].data/3600./np.cos(np.deg2rad(ydeg))
    ydeg_l = ydeg - alpha_y[ypix,xpix].data/3600.

    return (xdeg_l[i], ydeg_l[i])


def ima2abs_new(ra_deg, dec_deg, arc_name='PSZ1GA311_G1'):
    """ Return de-lensed coordinate for a given arc
    using the corresponding deflection matrix.
    (c) NT

    Parameters
    ----------

    ra_deg : np.array
    dec_deg : np.array

    arc_name : str
        currently available ['PSZ1GA311_G1', 'SGASJ1226_G1']

    """

    # number of positions
    N = len(ra_deg)
    # coordinates in image plane (deg)
    # xdeg=np.array([237.5082875])
    # ydeg=np.array([-78.1846194])
    xdeg=np.array([ra_deg])
    ydeg=np.array([dec_deg])

    # lens stuff, deflexion matrices, factor to normalize
    # lensfactor = d_ls/ds(z=0.73) / d_ls/ds(z=2.37)
    if 'PSZ1GA311' in arc_name:
        if arc_name == 'PSZ1GA311_G1':
            lensfactor = 0.4958 # depends on redshift
        elif arc_name == 'PSZ1GA311_G2':
            lensfactor = XXX
        else:
            raise ValueError('Not implemented for that arc_name.')

        path = os.getcwd()
        if path.startswith('/disk03/ntejos/'):
            prefix = '/disk03/ntejos/projects/arc_tomo/PSZ1GA311/data/delensing/'
        else:
            prefix = '/media/ntejos/disk1/projects/arc_tomo/PSZ1GA311/data/delensing/'
        alpha_x = lensfactor*Image(prefix + 'dplx1_Fel_wcs.fits')
        alpha_y = lensfactor*Image(prefix + 'dply1_Fel_wcs.fits')

    elif arc_name == 'SGASJ1226_G1':
        path = os.getcwd()
        if path.startswith('/disk03/ntejos/'):
            prefix = '/disk03/ntejos/projects/arc_tomo/SGASJ1226/data/delensing/'
        else:
            prefix = '/media/ntejos/disk1/projects/arc_tomo/SGASJ1226/data/delensing/'
        lensfactor = 1.  # it is supposed to be included in the matrices in this case
        alpha_x = lensfactor * Image(prefix + 'dplx_0.77_abs.fits')
        alpha_y = lensfactor * Image(prefix + 'dply_0.77_abs.fits')

    elif arc_name == 'SGASJ0033_G1':
        path = os.getcwd()
        if path.startswith('/disk03/ntejos/'):
            prefix = '/disk03/ntejos/projects/arc_tomo/SGASJ0033/data/delensing/'
        else:
            prefix = '/media/ntejos/disk1/projects/arc_tomo/SGASJ0033/data/delensing/'
        lensfactor = 0.53  # from Cedric ... need to be checked it is supposed to be included in the matrices in this case
        alpha_x = lensfactor * Image(prefix + 'SGAS0033_dx.fits')
        alpha_y = lensfactor * Image(prefix + 'SGAS0033_dy.fits')

    else:
        raise ValueError('Not implemented for that arc_name.')
    wcs = alpha_x.wcs

    # coordinates in pixels of alpha matrix
    xpix,ypix = [0]*N,[0]*N
    for i in range(N):
        ypix[i], xpix[i] = wcs.sky2pix((ydeg[i],xdeg[i]),nearest=True)[0]
        import pdb;    pdb.set_trace()
        # print("Pixel (x,y)=({},{})".format(xpix[i],ypix[i]))
    # map to absorber plane (using the alpha matrix)
    xdeg_l = xdeg + alpha_x[ypix,xpix].data/3600./np.cos(np.deg2rad(ydeg))
    ydeg_l = ydeg - alpha_y[ypix,xpix].data/3600.

    return (xdeg_l[i], ydeg_l[i])
















