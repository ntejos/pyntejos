import numpy as np
from astropy.io import ascii
from astropy.table import vstack
import os
"""Module for getting different Exposure Time Calculators"""


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


def s2n_COS(t, FUV, tnorm=1000):
    # load cos etc simulations for flat spectrum 1000s exp, and FUV = 17 mag
    data_path('cos_etc_g130m_v25.1.1.csv')
    cos_g130m = ascii.read(data_path('cos_etc_g130m_v25.1.1.csv'))
    cos_g160m = ascii.read(data_path('cos_etc_g160m_v25.1.1.csv'))
    #separate them at ~1450 A
    cond = cos_g130m['wavelength'] < 1450
    cos_g130m = cos_g130m[cond]
    cond = cos_g160m['wavelength'] >= 1450
    cos_g160m = cos_g160m[cond]
    # merge both
    cos = vstack([cos_g130m, cos_g160m], join_type='exact')

    # Signal
    signal = cos['target_counts'] * t / tnorm * 10.**((FUV- 17.)/(-2.5))

    # Noise terms
    dark = cos['dark_counts'] * t / tnorm
    sky = cos['sky_counts'] * t / tnorm

    # Noise
    var = signal + dark + sky
    sig = np.sqrt(var)

    #append S/N to cos
    sn = signal/sig * np.sqrt(3) # per-resolution element of 3 pixels
    cos['sn'] = sn
    return cos
