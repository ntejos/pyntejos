from barak.io import readtxt
import numpy as np
from astropy.table import Table, Column

def read_Danforth14(filename,keywords=['RA_TARG','DEC_TARG']):
    """Reads the special format of Danforth+14 on IGM ('igm-systems'
    style) catalos from HST/COS.
    
    Inputs:
    filename: name of the file
    keywords: list of keywords to retrieve appart from data. e.g. ['RA_TARG','DEC_TARG']
    
    Return: Astropy Table with original data + new columns having the
    given keywords.

    Example of file content from Danforth+14 IGM catalogs:

#SIMPLE  = T
#XTENSION= 'BINTABLE'
#BITPIX  = 8                       / 8-bit bytes.
#NAXIS   = 2                      / 2-dimensional table.
#NAXIS1 = 64                    / Width of table in bytes.
#NAXIS2  =           62             / Number of rows in table.
#PCOUNT  =   0                 /Random parameter count
#GCOUNT  =  1                  /Group count
#TELESCOP= 'HST'                  / Observatory name.
#INSTRUME= 'COS'                  / Instrument name.
#RADESYS = 'FK5'                  / Astrometric reference system.
#EQUINOX = 2000.                  / Equinox of the coordinates.
#TARGNAME = '1ES1028+511'           / Sight line name.
#TARGROOT = '1es1028'       / 7-character sight line id.
#RA_TARG  = 157.82708333   / RA of the target [deg].
#DEC_TARG =  50.89333333   / DEC of the target [deg].
#TARG_ZEM = 0.36040300             / Nominal emission redshift.
#DATE-OBS = '2011-05-01'    / Approximate date of first observation [YYYY-MM-DD].
#HLSPLEAD = 'Charles Danforth'    / HLSP Lead.
#CITATION = 'Danforth, et al. 2014, ApJ, XXX:XXX' / Reference for HLSP methodology.
#PR_INV_L = 'Danforth'            / PI Last Name.
#PR_INV_F = 'Charles'             / PI First Name.
#AIRORVAC = 'VAC     '        / Are wavelengths in air or vacuum? 
#FILTER01 = 'G130M'           / Filter used in coadd.
#FILTER02 = 'G160M'           / Filter used in coadd.
#EXPTIME1 =        14651          / Max exposure time in Filter 1 [sec].
#EXPTIME2 =        14606          / Max exposure time in Filter 2 [sec].
#SN_RES_1 = 20                / Median S/N per 7-pixel resel in Filter 1.
#SN_RES_2 = 12                / Median S/N per 7-pixel resel in Filter 2.
#TFIELDS = 16                 / Number of fields in each row.
#TTYPE1  = 'Z_SYS'            / 1: system redshift
#TFORM1  = '1E'
#TTYPE2  = 'DELTAV_SYS'       / 2: system velocity half-width
#TFORM2  = '1E'
#TUNIT2  = 'km/s'
#TTYPE3  = 'WAVELENGTH'       / 3: Observed wavelength
#TFORM3  = '1E'
#TUNIT3  = 'Ang'
#TTYPE4  = 'LINE_ID'          / 4: Line ID
#TFORM4  = '10A'
#TTYPE5  = 'z_ABS'            / 5: Line redshift
#TFORM5  = '1E'
#TTYPE6  = 'SIGLEVEL   '      / 6: Line Significance Level
#TFORM6  = '1E'
#TUNIT6  = 'SIGMA'
#TTYPE7  = 'S/N  '            / 7: Local S/N per resel 
#TFORM7  = '1E'
#TTYPE8  = 'EQW  '            / 8: Line Equiv. Width
#TFORM8  = '1J'
#TUNIT8  = 'mAA'
#TTYPE9  = 'EQW_ERR'          / 9: Equiv Width error
#TFORM9  = '1J'
#TUNIT9  = 'mAA'
#TTYPE10 = 'BVALUE'           / 10: Line b-value
#TFORM10 = '1E'
#TUNIT10 = 'km/s'
#TTYPE11 = 'BVALUE_ERR'       / 11: b-value error
#TFORM11 = '1E'
#TUNIT11 = 'km/s'
#TTYPE12 = 'LOGN'             / 12: line log column density
#TFORM12 = '1E'
#TUNIT12 = 'log cm^-2'
#TTYPE13 = 'LOGN_ERR'         / 13: log column density error
#TFORM13 = '1E'
#TUNIT13 = 'log cm^-2'
#TTYPE14 = 'FLAG_FIT'         / 14: line fit quality concerns?
#TFORM14 = '1I'
#TTYPE15 = 'NUM_SYS_LINES'    / 15: number of lines in system
#TFORM15 = '1I'
#TTYPE16 = 'NUM_METAL_LINES'  / 16: number of metal lines in system
#TFORM16 = '1I'
#COMMENT     
#COMMENT 'Generated Fri Jan 17 15:09:22 2014' 
#COMMENT 'Delivered to MAST from the COS Sight Lines HLSP project'
#END
0.002435 40.0  1218.61 "Lya 1215"   0.002435 3.2    3.0    81     2      22.1   1.0    13.27  0.01   1  1  0 
0.003204 65.0  1219.57 "Lya 1215"   0.003225 10.3   3.0    266    3      23.4   1.0    14.18  0.01   1  2  1 
0.003204 65.0  1553.07 "CIV 1548"   0.003141 3.9    6.0    88     14     36.9   10.9   13.39  0.09   1  2  1 
0.021890 39.0  1242.26 "Lya 1215"   0.021890 8.7    10.0   83     23     27.0   4.2    13.26  0.10   0  1  0 
0.022580 39.0  1243.10 "Lya 1215"   0.022580 6.7    10.0   86     11     50.6   8.3    13.24  0.05   0  1  0 
0.044865 39.0  1270.19 "Lya 1215"   0.044865 11.3   10.0   157    36     57.1   5.1    13.52  0.09   1  1  0 
0.050985 48.0  1277.63 "Lya 1215"   0.050985 23.0   9.0    206    13     23.1   1.5    13.89  0.03   0  1  0 
0.051399 38.0  1278.13 "Lya 1215"   0.051399 3.5    9.0    32     17     18.1   8.0    12.82  0.10   1  1  0 
0.080392 37.0  1313.38 "Lya 1215"   0.080392 4.3    7.0    41     23     9.6    23.5   13.00  0.20   0  1  0 
0.083274 37.0  1316.88 "Lya 1215"   0.083274 4.5    7.0    126    68     100.0  1.0    13.40  0.23   1  1  0 
0.088340 37.0  1323.04 "Lya 1215"   0.088340 6.5    7.0    105    46     36.6   6.0    13.36  0.16   0  1  0 
0.094355 37.0  1330.35 "Lya 1215"   0.094355 9.6    8.0    124    28     29.8   3.1    13.47  0.09   1  1  0 
"""
    #open file
    f = open(filename)
    
    #to store keyword values
    values = dict()
    
    #to store column names
    colnames = []
 
    while True:
        line = f.readline()
        words = line.split()
        
        if (words[0][0] != '#') or (len(keywords)==0):
            break
        
        #read keyword values
        for keyword in keywords:
            if words[0][1:]==keyword:
                try:
                    values[keyword] = float(words[2]) #store the value of keyword
                except:
                    values[keyword] = words[2][1:-1] #store the value of keyword
        #read colnames
        if words[0][1:].startswith('TTYPE'):
            name = words[2][1:-1]
            if name == 'S/':
                name = 'SN'
            elif name == 'EQ':
                name = 'EQW'
            elif name == 'SIGLEVE':
                name = 'SIGLEVEL'
            colnames = colnames + [name] #store colname
            if name == 'LINE_ID': #this is to account for the
                                            #fact that "Lya 1215" is
                                            #read as two columns later:
                                            #'"Lya' and '1215"'
                colnames = colnames + ['LINE_ID_WAVELENGTH']
    f.close()
    
    #read content
    data = readtxt(filename,comment='#',names=colnames)
    
    #clean excess of " in data accosiated to 'LINE_ID' and
    #'LINE_ID_WAVELENGTH'
    try:
        aux = data['LINE_ID']
        aux = np.array([aux[i][1:] for i in range(len(aux))])
        data['LINE_ID'] = aux
        aux = data['LINE_ID_WAVELENGTH']
        aux = np.array([aux[i][:-1] for i in range(len(aux))])
        data['LINE_ID_WAVELENGTH'] = aux
    except:
        pass
    #use data in astropy.table Table format
    data = Table(data)

    #add columns for each keyword
    for key in values.keys():
        aux = [values[key] for i in range(len(data)) ]
        aux = Column(name=key,data=aux)
        data.add_column(aux)
        
    return data


def write_DS9reg(x, y, filename=None, coord='IMAGE', ptype='x', size=20,
                 c='green', tag='all', width=1, text=None):
    """Write a region file for ds9 for a  list of coordinates.
    Taken from Neil Crighton's barak.io

    Parameters
    ----------
    x, y : arrays of floats, shape (N,)
        The coordinates. These may be image or WCS.
    filename : str, optional
        A filename to write to.
    coord : str  (`IMAGE` or `J2000`)
        The coordinate type: `IMAGE` (pixel coordinates) or
        `J2000` (celestial coordinates).
    ptype : str or np.array of shape (N,)
        DS9 point type (e.g. `circle`, `box`,  `diamond`,  `cross`, `x`, `arrow`,
        `boxcircle`)
    size : int or np.array of shape (N,)
        DS9 point size.
    c : str or np.array of shape (N,)
        point colour: `cyan` `blue` `magenta` `red` `green` `yellow` `white`
        `black`}.
    tag : str or np.array of shape (N,)
        DS9 tag. e.g. 'all'
    width : int or np.array of shape (N,)
        DS9 width
    text : str or np.array of shape (N,)
        Text
    """
    header = ['global font="helvetica 10 normal" select=1 highlite=1 '
               'edit=0 move=1 delete=1 include=1 fixed=0 source\n']
    header.append(coord + '\n')

    x = np.array(x)
    y = np.array(y)
    if isinstance(ptype, basestring):
        ptype = [ptype] * len(x)
    if isinstance(size, int):
        size = [size] * len(x)
    if isinstance(width, int):
        width = [width] * len(x)
    if isinstance(text, basestring):
        text = [text] * len(x)
    elif text is None:
        text = list(range(len(x)))
    if isinstance(tag, basestring):
        tag = [tag] * len(x)
    if isinstance(c, basestring):
        c = [c] * len(x)
    
    
    regions = []
    # fmt = ('point(%12.8f,%12.8f) # \
    # point=%s %s width=%s text={%s} color=%s tag={%s}\n')
    for i in xrange(len(x)):
        s = 'point({:.8f},{:.8f}) # point={} {} width={} text={{{}}} color={} tag={}\n'\
                .format(x[i], y[i], ptype[i], size[i], width[i], text[i], c[i], tag[i])
        regions.append(s)

    if filename is not None:
        fh = open(filename,'w')
        fh.writelines(header + regions)
        fh.close()
    return header, regions

#def read_ascii(filename,names=None,):
#    """Read ascii file table"""

def read_publication_xml(filename):

    """Produces a table from the ADS xml file

    Returns
    -------
    tab : astropy Table
        Table version of the xml file with most relevant values

    """


    # first read the file
    f = open(filename, 'r')
    lines = f.readlines()
    reset_author = False
    author_aux = []
    authors = []
    bibcode = []
    title = []
    citations = []
    pubdate = []
    journal = []
    citation_exists = False
    for ii, line in enumerate(lines):
        if '<bibcode>' in line:
            aux = line.split('<bibcode>')[1].split('</bibcode>')[0]
            bibcode += [aux]
            if 'arXiv' in line:
                journal += ['arXiv']
            else:
                aux = bibcode[-1][4:].split('.')[0]
                if 'prop' in line:
                    aux = aux + '.prop'
                journal += [aux]
        elif '<title>' in line:
            aux = line.split('<title>')[1].split('</title>')[0]
            title += [aux]
        elif '<citations>' in line:
            aux = line.split('<citations>')[1].split('</citations>')[0]
            citations += [aux]
            citation_exists = True
        elif '<pubdate>' in line:
            aux = line.split('<pubdate>')[1].split('</pubdate>')[0]
            pubdate += [aux]

        # authors
        if reset_author is True:
            author_aux = []
            reset_author = False
        if '<author>' in line:
            aux = line.split('<author>')[1].split('</author>')[0]
            author_aux += [aux]
        elif '</record>' in line: # new article
            authors += [author_aux]
            reset_author = True
            if citation_exists is False:
                citations += [0]
            citation_exists = False

    tab = Table()
    tab['ind'] = range(0,len(bibcode))
    tab['bibcode'] = bibcode
    tab['pubdate'] = pubdate
    tab['year'] = [int(b[:4]) for b in tab['bibcode']]
    tab['citations'] = np.array(citations).astype(int)
    # tab['title'] = title  # something wrong with encoding
    tab['authors'] = authors
    tab['n_authors'] = [len(authors) for authors in tab['authors']]
    tab['1author'] = [a[0].split(',')[0] for a in tab['authors']]
    tab['name'] = ['{}+{}'.format(fa,yr) for fa, yr in zip(tab['1author'], tab['year'])]
    return tab
