#!/usr/bin/env python

import numpy as np
import json
import argparse
import sys
from astropy.constants import c
from astropy.table import Table, Column
import astropy.units as u

"""This script reads a .json file given by IGMguesses.py and writes a
component table.
"""

# get arguments

parser = argparse.ArgumentParser(description='Parser for igmguesses_json_to_comptable.py')
parser.add_argument("in_file", type=str, help="IGMguesses .json file")
parser.add_argument("out_file", type=str, help="Output filename (ascii table)")

args = sys.argv[1:]
pargs = parser.parse_args(args)

# Load in json file and get relevant quantities
jf = open(pargs.in_file)
data = json.load(jf)
comps = data['cmps']

# create columns as lists first
name = []
ion = []
comment = []
zcmp = []
zfit = []
vmin_zcmp = []
vmax_zcmp = []
vmin_zfit = []
vmax_zfit = []
wrest = []
logn = []
b = []
quality = []

for comp_name in comps.keys():
    name += [comp_name]
    ion += [comp_name.split('_')[1]]
    comment += [comps[comp_name]['Comment']]
    wrest += [comps[comp_name]['wrest']]
    zfit_aux = comps[comp_name]['zfit']
    zfit += [zfit_aux]  # component redshift from fit
    zcmp_aux = comps[comp_name]['zcomp']
    zcmp += [zcmp_aux]  # middle of component vlim
    vlim_zcmp = np.array((comps[comp_name]['vlim']))  # w/r zcmp
    vmin_zcmp += [vlim_zcmp[0]]
    vmax_zcmp += [vlim_zcmp[1]]

    zlim_zcmp = zcmp_aux + vlim_zcmp * (1 + zcmp_aux) / c.to('km/s').value  # w/r zcmp
    vlim_zfit = (zlim_zcmp - zfit_aux) * c.to('km/s').value / (1 + zfit_aux)  # w/r zfit
    vmin_zfit += [vlim_zfit[0]]
    vmax_zfit += [vlim_zfit[1]]

    # quantities that are independent of redshift
    logn += [comps[comp_name]['Nfit']]
    b += [comps[comp_name]['bfit']]
    quality += [comps[comp_name]['Quality']]

# Create table
comp_table = Table()
comp_table.add_column(Column(name='name', data=name))
comp_table.add_column(Column(name='ion', data=ion))
comp_table.add_column(Column(name='comment', data=comment))
comp_table.add_column(Column(name='zcmp', data=zcmp))
comp_table.add_column(Column(name='vmin_zcmp', data=vmin_zcmp, unit=u.km/u.s))
comp_table.add_column(Column(name='vmax_zcmp', data=vmax_zcmp, unit=u.km/u.s))
comp_table.add_column(Column(name='zfit', data=zfit))
comp_table.add_column(Column(name='vmin_zfit', data=vmin_zfit, unit=u.km/u.s))
comp_table.add_column(Column(name='vmax_zfit', data=vmax_zfit, unit=u.km/u.s))
comp_table.add_column(Column(name='wrest', data=wrest, unit=u.AA))
comp_table.add_column(Column(name='logn', data=logn))
comp_table.add_column(Column(name='b', data=b, unit=u.km/u.s))
comp_table.add_column(Column(name='quality', data=quality))

# Write table
comp_table.write(pargs.out_file, format='ascii', delimiter='|')
