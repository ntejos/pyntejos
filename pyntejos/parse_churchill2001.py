import numpy as np
from astropy.table import Table
from astropy import units as u
from linetools.spectralline import AbsLine


"""This module parses the ASCII tables of
Churchill 2001 paper"""


def replace_upper(s1):
    while '<' in s1:
        ind = s1.find('<')
        aux = s1[ind:ind + 6]
        number = aux[1:]
        s1 = s1.replace(aux, number + '\kern.1666em\pm \kern.1666em-1\t')
        # import pdb;pdb.set_trace()
    return s1


def replace_sim(s1):
    while '\sim' in s1:
        ind = s1.find('\sim')
        aux = s1[ind:ind + 9]
        number = aux[5:]
        s1 = s1.replace(aux, number + '\kern.1666em\pm \kern.1666em-2\t')
    return s1


def parse_tab7(fn_tab7):
    # if 1:
    f = open(fn_tab7)
    lines = f.readlines()
    qso = []
    zabs = []
    vel = []
    mgii_logN = []
    mgii_logN_err = []
    mgii_b = []
    mgii_b_err = []
    feii_logN = []
    feii_logN_err = []
    feii_b = []
    feii_b_err = []
    mgi_logN = []
    mgi_logN_err = []
    mgi_b = []
    mgi_b_err = []

    for ii, line in enumerate(lines):
        if line == '\n':
            continue
        # replace \ldots by a number
        line = line.replace('\ldots', '0.0\kern.1666em\pm \kern.1666em0.0')
        if ii == 169:  # fix Churchill table inconsistency
            line = line.replace('10.70', '<10.70')

        # replace < limits and \sim
        line = replace_upper(line)
        line = replace_sim(line)
        # fo.write(line)
        items = line.split('\kern.1666em')
        # print(len(items), line.count('\ldots'), ii)

        q = items[0].split('\t')[0]
        z = items[0].split('\t')[1]
        if q == '':
            q = qso[ii - 2]
        if z == '':
            z = zabs[ii - 2]
        qso += [q]
        zabs += [z]
        # clean items and remove double \t\t if exists
        for jj in range(len(items)):
            if jj == 0:
                continue
            items[jj] = items[jj].replace('\t\t', '\t')

        vel += [items[0].split('\t')[2]]
        mgii_logN += [items[0].split('\t')[3]]
        mgii_logN_err += [items[2].split('\t')[0]]
        mgii_b += [items[2].split('\t')[1]]
        mgii_b_err += [items[4].split('\t')[0]]
        feii_logN += [items[4].split('\t')[1]]
        feii_logN_err += [items[6].split('\t')[0]]
        # import pdb; pdb.set_trace()
        feii_b += [items[6].split('\t')[1]]
        feii_b_err += [items[8].split('\t')[0]]
        mgi_logN += [items[8].split('\t')[1]]
        mgi_logN_err += [items[10].split('\t')[0]]
        mgi_b += [items[10].split('\t')[1].replace('\n', '')]
        mgi_b_err += [items[12].split('\n')[0].replace('\t', '')]

    # structure the table
    tab = Table()
    tab['ind'] = [ii for ii in range(len(qso))]
    tab['qso'] = qso
    tab['zabs'] = np.array(zabs).astype(float)
    tab['vel'] = np.array(vel).astype(float)
    tab['mgii_logN'] = np.array(mgii_logN).astype(float)
    tab['mgii_logN_err'] = np.array(mgii_logN_err).astype(float)
    tab['mgii_b'] = np.array(mgii_b).astype(float)
    tab['mgii_b_err'] = np.array(mgii_b_err).astype(float)
    tab['feii_logN'] = np.array(feii_logN).astype(float)
    tab['feii_logN_err'] = np.array(feii_logN_err).astype(float)
    tab['feii_b'] = np.array(feii_b).astype(float)
    tab['feii_b_err'] = np.array(feii_b_err).astype(float)
    tab['mgi_logN'] = np.array(mgi_logN).astype(float)
    tab['mgi_logN_err'] = np.array(mgi_logN_err).astype(float)
    tab['mgi_b'] = np.array(mgi_b).astype(float)
    tab['mgi_b_err'] = np.array(mgi_b_err).astype(float)
    return tab
