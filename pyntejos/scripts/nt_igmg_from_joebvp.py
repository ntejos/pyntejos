#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)


try:
    ustr = unicode
except NameError:
    ustr = str

"""This script creates a IGMGuesses JSON file from a JOEBVP input file.
"""
def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Coadds COS x1d.fits files.')
    parser.add_argument("filename", type=str, default=None, help="Input JOEBVP file.")
    parser.add_argument("specfile", type=str, default=None, help="Name of the spectrum file.")
    parser.add_argument("-o", type=str, default='out_igmg.json', help="Output filename (default out_igmg.json)")
    parser.add_argument("-fwhm", type=int, default=3, help="FWHM of the spectrum (default is 3)")
    parser.add_argument("-radec", type=str, default="J000000+000000", help="Coordinates of the object in J2000 (default is J000000+000000)")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from pyntejos.utils import igmgjson_from_joebvp
    from astropy.coordinates import SkyCoord
    from linetools.utils import radec_to_coord

    pargs = parser(options=args)
    filename = pargs.filename
    specfile = pargs.specfile
    output = pargs.o # output
    fwhm = pargs.fwhm
    radec = radec_to_coord(parser.radec)

    igmgjson_from_joebvp(filename, radec, specfile, fwhm, output)


if __name__ == '__main__':
    main()