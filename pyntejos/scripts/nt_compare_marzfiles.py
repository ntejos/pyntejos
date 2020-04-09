#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    ustr = unicode
except NameError:
    ustr = str

"""This scripts compares two MARZ files, and write info in the terminal
(c) NT 2020.
"""

def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Compares two MARZ files .mz and prints the differences in the terminal')
    parser.add_argument("marzfile1", type=str, default=None, help="Marzfile 1")
    parser.add_argument("marzfile2", type=str, default=None, help="Marzfile 2")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    from pyntejos.utils import compare_2_marzfiles

    pargs = parser(options=args)
    mz1 = pargs.marzfile1
    mz2 = pargs.marzfile2
    compare_2_marzfiles(mz1, mz2)


if __name__ == '__main__':
    main()
