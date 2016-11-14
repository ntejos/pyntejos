#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    ustr = unicode
except NameError:
    ustr = str

"""This script reads (ra,dec) coordinates and matches them to
catalogs of known objects.
"""
def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Matches input radec coordinates to a given catalog within a given angular separation.')
    parser.add_argument("radec", nargs='?', type=str, default=None, help="RA,DEC (e.g. 152.25900,7.22885), JXX (e.g. J100902.16+071343.8)")
    parser.add_argument("angsep", type=float, help="Angular separation in arcmin")
    parser.add_argument("--catalog", type=str, default='QSO', help="Catalog to match the coordinate to (default is 'QSO')")
    parser.add_argument("--epoch", default=2000., type=float, help="Epoch [Not functional]")
    parser.add_argument("--redshift", default=None, help="Redshift of the source; if given co-moving distance is calculated.")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from linetools.scripts.utils import coord_arg_to_coord
    from linetools import utils as ltu
    from astropy.io import fits, ascii
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from pyntejos.catalogs import add_radec_deg_columns

    pargs = parser(options=args)
    # RA,DEC
    icoord = coord_arg_to_coord(pargs.radec)
    coord = ltu.radec_to_coord(icoord)

    # define catalog
    print('Reading {} catalog'.format(pargs.catalog))
    if pargs.catalog == 'QSO':
        # read qsos from MILLIQUAS catalog
        col_names = ['ra_d', 'dec_d', 'name', 'description', 'rmag', 'bmag', 'comment', 'psf_r', 'psf_b', 'z', 'cite', 'zcite', 'qso_prob', 'Xname', 'Rname', 'Lobe1', 'Lobe2'] 
        cat = ascii.read('/media/ntejos/disk1/catalogs/qsos/milliquas/milliquas.txt', format='fixed_width', names=col_names)
    elif pargs.catalog == 'GC':
        # read MW globular cluster catalog 
        cat = ascii.read('/media/ntejos/disk1/catalogs/globular_clusters/mwgc10_1.dat', format='fixed_width')
        # add ra_d, dec_d columns
        cat = add_radec_deg_columns(cat)
    else:
        print(' Not implemented for such catalog.')
        return

    # cross-match
    print('Cross-matching...')
    cat_coords = SkyCoord(cat['ra_d'], cat['dec_d'], unit='deg')
    seplim = pargs.angsep * u.arcmin
    sep2d = coord.separation(cat_coords)
    cond = sep2d <= seplim
    cat = cat[cond]
    if len(cat) < 1:
        print("No matches found.")
    else:
        cat['sep2d'] = sep2d[cond]
        if pargs.redshift is not None:
            from astropy.cosmology import Planck15 as cosmo
            import pdb; pdb.set_trace()
            sep = (cosmo.kpc_comoving_per_arcmin(float(pargs.redshift)) * cat['sep2d']).to('Mpc')
            cat['sep_mpc'] = sep.value
        cat.sort('sep2d')
        print(cat)

if __name__ == '__main__':
    main()