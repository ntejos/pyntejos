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
    parser.add_argument("configfile", type=str, default=None, help="Configuration file. To dump an "
                                                                   "example version use -d option.")
    parser.add_argument("-d", type=str, default=None, help="Use -d to dump a dummy configuration file "
                                                           "config_make_MagE_cube_dummy.json.")
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    from pyntejos import mage

    pargs = parser(options=args)
    config_file = pargs.configfile
    # asking for dummy config file?
    if pargs.d is not None:
        mage.write_config_make_MagE_cube_dummy("./config_make_MagE_cube_dummy.json")
        return
    else:
    mage.make_MagE_cube(pargs.configfile)


if __name__ == '__main__':
    main()
