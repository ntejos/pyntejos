#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    ustr = unicode
except NameError:
    ustr = str

"""This scripts prints the offset between position 1 and reference position
(c) NT 2019.
"""

def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Prints the offset in arcsec of a Position 1 w/r to a reference.')
    parser.add_argument("pos1", type=str, default=None, help="Position 1 in hhmmssddmmss")
    parser.add_argument("pos_ref", type=str, default=None, help="Reference position in hhmmssddmmss")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    from pyntejos.coord import get_offsets

    pargs = parser(options=args)
    pos1 = pargs.pos1
    pos_ref = pargs.pos_ref
    get_offsets(pos1, pos_ref)

if __name__ == '__main__':
    main()