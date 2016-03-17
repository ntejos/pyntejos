import numpy as np
from astropy.io import fits
from astropy.cosmology import Planck13 as cosmo
from pyntejos.catalogs import unify_cluster_catalog_redmapper
from create_cpair_catalog_utils import create_pairs
import sys
"""This script is to create a full cluster-pair catalog from the
redmapper DR8 or DES_SVA1 cluster catalog of Rykoff et al."""

# Global parameters are read from configuration file
from run_all_config import params as param

# read stuff from the system and replace those from the config file
filename = sys.argv[1]
dataset = sys.argv[2]
param['filename'] = filename
param['dataset'] = dataset

# read clusters from redmapper
if param['dataset'] == 'sdss_dr8':
    redmapper = fits.getdata('/home/ntejos/catalogs/clusters/redmapper/SDSS_DR8/redmapper_dr8_public_v5.2_catalog.fits')
    # redmapper = fits.getdata('/home/ntejos/catalogs/clusters/redmapper/SDSS_DR8/dr8_run_redmapper_v5.10_lgt5_catalog.fit')
if param['dataset'] == 'des_sva1':
    redmapper = fits.getdata('/media/ntejos/disk1/catalogs/clusters/redmapper/DES_SVA1/redmapper_sva1_public_v6.3_catalog.fits')
    # relax the spectroscopic criteria for DES?
    # opt.z_err_max1 = 0.005
else:
    redmapper = fits.getdata(param['filename'])

# unify cluster catalog for right columns and naming conventions
clusters = unify_cluster_catalog_redmapper(redmapper, cosmo)

#add dataset
clusters['dataset'] = param['dataset']

# get rid of low and high redshift clusters and low-richness ones
cond = (clusters['redshift'] >= param['zmin']) & (clusters['redshift'] < param['zmax']) & (clusters['richness'] > param['richmin'])
clusters = clusters[cond]

# create unique identifiers
clusters['objid'] = [i+1 for i in range(len(clusters))]

# create cluster-pairs
cpairs = create_pairs(clusters['ra'], clusters['dec'], clusters['redshift'],
                clusters['redshift_err'], cosmo, max_tsep=param['max_tsep'],
                z_err_max1=param['z_err_max1'], z_err_max2=param['z_err_max2'],
                objid=clusters['objid'], dvmax=param['max_dv'], verbose=True)

# write table
name = 'cpairs_richmin{}_l{}_{}'.format(int(param['richmin']), int(param['max_tsep']), param['dataset'])
cpairs.write(name+'.fits', format='fits', overwrite=True)
# cpairs.write(name+'.txt', format='ascii')
