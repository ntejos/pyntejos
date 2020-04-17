import urllib
from IPython import embed
"""This module should contain functions to parse specific formats"""

def write_lsf_linetools_from_stscifile(stsci_file, output):
    """reads a stsci LSF file for COS, and writes the proper format that linetools
    expects.
    It currently works for the table versions at \

    www.stsci.edu/hst/instrumentation/cos/performance/spectral-resolution

    available on 16 april 2020... (they sometimes change the format)
    """

    f = open(stsci_file, 'r')
    fo = open(output, 'w')
    lines = f.readlines()

    # parse first line
    l0 = lines[0]
    l0_new = l0.replace(' ', 'A,')
    l0_new = l0_new.replace('\n', 'A\n')
    l0_new = 'pixel,'+l0_new
    fo.write(l0_new)

    for ii, line in enumerate(lines[1:]):
        line_new = '{}\t'.format(ii+1) + line
        fo.write(line_new)
    fo.close()
    f.close()
