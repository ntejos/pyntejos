#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    ustr = unicode
except NameError:
    ustr = str

"""This scripts creates a "cube" version for a single MagE slit.
(c) NT 2018.
"""

def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Makes a MagE datacube based on individual MagE spectra of a single MagE slit.')
    parser.add_argument("input", type=str, default=None, help="Input Mage cube .fits")
    parser.add_argument("-airvac", type=str, default='air', help="Whether thee file is calibrated in air (air) or vaccumm (vac). Default is 'air'")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    from pyntejos import mage
    pargs = parser(options=args)
    mage.write_magecube_as_xspectrum1d(pargs.input, airvac=pargs.airvac)

if __name__ == '__main__':
    main()
