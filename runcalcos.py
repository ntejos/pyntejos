""" Routine for running CalCOS """

import calcos
import pyfits
import glob


def run(rawdir, outdir):

    print "Fetching input files"
    rawtag = glob.glob(rawdir + "*rawtag*")
    asn = glob.glob(rawdir + "*asn*")

    print "Setting RANDSEED keyword value to 1"
    for i in range(len(rawtag)):
        hdulist = pyfits.open(rawtag[i], mode="update")
        prihdr = hdulist[0].header
        prihdr.update("RANDSEED", 1)
        hdulist.flush()

    print "Running CalCOS..."
    for i in range(len(asn)):
        calcos.calcos(asn[i], verbosity=2, outdir=outdir)
