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
        description='Coadds STIS echelle x1d.fits files. It gathers all the individual orders among all '
                    'the different exposures (including multi-extension files) and rebins to native sampling'
                    'only once.')
    parser.add_argument("filenames", nargs='+', type=str, default=None, help="Filenames to use for co-addition "
                                                                             "(e.g. o6n403010_x1d.fits, o6n403020_x1d.fits, etc.")
    parser.add_argument("-rebin", type=int, default='3', help="Number of native pixels to rebin (default is 3).")
    parser.add_argument("-o", type=str, default='out.fits', help="Output filename (default out.fits)")
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from pyntejos import coadd

    pargs = parser(options=args)
    filenames = pargs.filenames
    output = pargs.o # output
    rebin = pargs.rebin
    if rebin < 2:
        rebin = None

    # co-add
    spec = coadd.coadd_stis_from_x1dfiles(filenames, wv_array=None, rebin = rebin)
    spec.write_to_fits(output)


if __name__ == '__main__':
    main()