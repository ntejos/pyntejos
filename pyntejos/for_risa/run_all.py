#!/usr/bin/env python

import sys
import os

"""This python script runs both create_cpairs_catalog.py and
cross_match the created catalog with """

filename = sys.argv[1]
dataset = sys.argv[2]

print filename
print dataset

print 'Creating the cluster-pair catalog'
os.system('python create_cpair_catalog.py {} {}'.format(filename, dataset))
print 'Cross-matching the cluster-pair catalog to a sample of bright UV QSOs'
os.system('python cross_match_cpairs_qsos.py {}'.format(dataset))