import numpy as np
from pyntejos import utils as ntu
"""Utilities for CONICYT related stuff. Most noticeable, functions to calculate
the score of publications"""


def add_pos_author(tab, lastname):
    """For a given publication table,
    it adds number within author list where
    string `lastname` appears, and returns a
    new table with that extra column"""

    counter = []
    for row in tab:
        author_list = row['authors']
        q = 0
        for a in author_list:
            if lastname in a:
                counter += [q+1]
                break
            q += 1
        q += 1
        if q > len(author_list):
            print('{}/{}'.format(q, len(author_list)))
            print('Warning: {} does not appear as author in publication {}'.format(lastname, row['name']))
            counter += [0]
            print(row['ind'], counter)
            # import pdb; pdb.set_trace()
    new_tab = tab.copy()
    new_tab['pos_{}'.format(lastname)] = counter
    return new_tab


def get_ci_factor(tab):
    """C_i number of citations divided by number of years since publication"""
    current_year = ntu.get_current_year()
    deno = np.where(current_year==tab['year'], 1., (current_year - tab['year']))
    ci = tab['citations'] / deno.astype(float)
    return ci


def is_published(tab):
    """Returns boolean array on whether a paper has been published"""
    is_pub = []
    nopub = ['prop','yCat', 'soft', 'conf', 'IAUS', '2018arXiv180305005P']
    for row in tab:
        bibcode = row['bibcode']
        q = 0
        for aux in nopub:
            if aux in bibcode:
                is_pub += [False]
                q = 1
                break
        if q == 0:
            is_pub += [True]
    is_pub = np.array(is_pub)
    return is_pub


def get_li_factor(tab, lastname, conc='ini19'):
    """L_i factor defined in FONDECYT Regular o Iniciacion 2019"""
    li = [0.2]*len(tab) # base
    li = np.array(li)
    author_pos = tab['pos_{}'.format(lastname)]
    author_pos = np.array(author_pos)
    if conc == 'ini19':
        li = np.where(author_pos == 5, 0.3, li)
        li = np.where(author_pos == 4, 0.4, li)
        li = np.where(author_pos == 3, 0.6, li)
        li = np.where(author_pos == 2, 0.8, li)
        li = np.where(author_pos == 1, 1.0, li)
        li = np.where(author_pos == 0, 0.0, li)
    elif conc == 'reg19':
        li = np.where(author_pos == 6, 0.3, li)
        li = np.where(author_pos == 5, 0.5, li)
        li = np.where(author_pos == 4, 0.7, li)
        li = np.where(author_pos == 3, 0.9, li)
        li = np.where(author_pos <= 2, 1.0, li)
        li = np.where(author_pos == 0, 0.0, li)
    else:
        raise ValueError('Not implemented for this concurso: {}.'.format(conc))
    return li


def get_si_factor(ci, li, conc='ini19'):
    """S_i factor defined in FONDECYT Regular or Iniciacion"""
    # import pdb; pdb.set_trace()
    if conc in ['ini19', 'reg19']:
        si = li*((1 + ci)**0.5)
    else:
        raise ValueError('Not implemented for this concurso: {}.'.format(conc))
    return si


def get_puntaje(S, NA, conc='ini19'):
    if conc == 'ini19':
        punt = 1. + 2.7*((S/(NA+2))**0.25)
    elif conc == 'reg19':
        punt = 1. + 1.7*(S**0.25)  # does not depend on NA
    else:
        raise ValueError('Not implemented for this concurso: {}.'.format(conc))
    return punt

