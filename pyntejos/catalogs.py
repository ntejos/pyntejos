import numpy as np
from astropy.table import Table, Column, Row, vstack
import astropy.units as u
from pyntejos.coord import dec2s

"""This module is meant to provide catalog utilities,
for instance, unifying different public catalogs
"""

def give_name(ra,dec):
    """Returns a numpy array having a string with name of the target based
    on ra,dec in degrees.

    """
    ra_s,dec_s = dec2s(ra,dec)
    ra_s = ra_s.replace(' ','')[:4]
    dec_s = dec_s.replace(' ','')[:5]
    name = 'J{}{}'.format(ra_s,dec_s)
    return name

def unify_cluster_catalog_redmapper(catalog, dataset, cosmo):
    """For a given catalog of redMaPPer clusters, it unifies its add_columns
    and column names (ra, dec, redshift, richness, mass, etc) for 
    general analysis. 

    Parameters
    ----------
    catalog : Table
        The redmapper catalog as an astropy Table object
    dataset : str
        Choose between: `sdss_dr8` or `des_sva1`
    cosmo : Cosmology object
        Cosmology from astropy.cosmology
    Returns
    -------
    catalog_out : Table
        An astropy Table object with the unified redmapper catalog
    """
    
    clusters = Table(catalog)
    
    #Unify name and dataset
    if 'NAME' not in clusters.keys():
        # This is probably an old format
        try:
            name = [str(catalog['MEM_MATCH_ID'][i]) for i in range(len(clusters))] 
        except:
            raise ValueError('This is probably an old format of redmapper catalog; please use a more recent one.')
    else:
        name = catalog['NAME']
    clusters['NAME_ORIG'] = name
    clusters['dataset'] = dataset
    
    #Unify ID, add column with original indices for redmapper
    objid = np.array([int(i+1) for i in range(len(catalog))])
    clusters['objid_redmapper'] = objid
    
    #Unify RA/DEC
    ra  = catalog['RA']
    dec = catalog['DEC']
    clusters['ra'] = ra 
    clusters['dec'] = dec
    
    #unify best redshifts (spec-z when available, otherwise photo-z)
    if 'Z_SPEC' not in clusters.keys():
        # This is probably an old format
        try:
            clusters.rename_column('BCG_SPEC_Z', 'Z_SPEC')
            clusters.rename_column('Z_LAMBDA_E', 'Z_LAMBDA_ERR')
            clusters.rename_column('LAMBDA_CHISQ', 'LAMBDA') # this is richness
        except:
            import pdb; pdb.set_trace()
            raise ValueError('This is probably an old format of redmapper catalog; please use a more recent one.')
    z = np.where(clusters['Z_SPEC']!=-1,clusters['Z_SPEC'],clusters['Z_LAMBDA'])
    z_err = np.where(clusters['Z_SPEC']!=-1,0.0002,clusters['Z_LAMBDA_ERR'])
    zsp = np.where(clusters['Z_SPEC']!=-1,clusters['Z_SPEC'],-1)
    zph = clusters['Z_LAMBDA']
    zph_err = clusters['Z_LAMBDA_ERR']
    flag_zsp = (clusters['Z_SPEC'] != -1)*1.
    clusters['redshift'] = z
    clusters['redshift_err'] = z_err
    clusters['flag_zsp'] = flag_zsp    
    clusters['zph'] = zph
    clusters['zph_err'] = zph_err
    clusters['zsp'] = zsp
    
    #create unified richness
    richness = clusters['LAMBDA']
    clusters['richness'] = richness 

    # add mass estimate for redmapper
    mass = (10**14 / (cosmo.H0.value/70) ) * np.exp(1.72) * np.power(clusters['richness']/60.,1.08)
    clusters['mass'] = mass
    
    # r200 for clusters
    rho_200 = 200 * cosmo.critical_density(clusters['redshift']).to('Msun/Mpc3') #in Msun/Mpc3
    mass_200 = clusters['mass']*u.Msun # in Msun
    r_200 = np.power(mass_200/rho_200/(4*np.pi/3.), 1/3.).value #in Mpc
    clusters['r200_mpc'] = r_200

    # give unified name
    name = [give_name(ra,dec) for ra,dec in zip(clusters['ra'],clusters['dec'])]
    clusters['name'] = name

    return clusters

def unify_cluster_catalog_GMBCG(catalog, cosmo):
    """For a given catalog of GMBCG clusters, it unifies it to the
    required properties (RA, DEC, redshift, richness, Mass, etc) as
    used in our analysis. It returns a new catalog of only relevant
    information."""
    
    clusters = Table()
    
    #Unify name
    name = Column(name='NAME_ORIG',data=catalog['NAME'])
    clusters.add_column(name)
    
    #Unify ID, add column with original indices for hao
    objid = np.array([int(i+1) for i in range(len(catalog))])
    objid = Column(name = 'objid_gmbcg', data=objid)
    clusters.add_column(objid)
    
    #Unify RA/DEC
    ra  = Column(name='ra' ,data=catalog['RA'])
    dec = Column(name='dec',data=catalog['DEC'])
    clusters.add_columns([ra,dec])
    
    #unify best redshifts (spec-z when available, otherwise photo-z)
    z = np.where(catalog['SPZ']!=0,catalog['SPZ'],catalog['PHOTOZ'])
    z = Column(name = 'redshift', data=z)
    clusters.add_column(z)
    z_err = np.where(catalog['SPZ']!=0,0.0002,catalog['PHOTOZ_ERR'])
    z_err = Column(name='redshift_err',data=z_err)
    clusters.add_column(z_err)
    
    #add photometric and spectroscopic redshifts
    zsp = Column(name='zsp',data=np.where(catalog['SPZ']>0,catalog['SPZ'],-1))
    zph = Column(name='zph',data=catalog['PHOTOZ'])
    zph_err = Column(name='zph_err',data=catalog['PHOTOZ_ERR'])
    clusters.add_columns([zsp,zph,zph_err])
    
    #add flag_spec
    flag_zsp = (catalog['SPZ'] > 0)*1.
    flag_zsp = Column(name='flag_zsp',data=flag_zsp)
    clusters.add_column(flag_zsp)

    #create unified richness for hao
    richness = np.where(catalog['WEIGHTOK'],catalog['GM_NGALS_WEIGHTED'],catalog['GM_SCALED_NGALS'])
    richness = Column(name = 'richness', data=richness)
    clusters.add_column(richness)

    #add mass estimate
    #add column with cluster mass using calibration given by http://gmbcg.blogspot.com/
    #cluster mass = 8.8e13*((0.66*Ngals+9.567)/20.)^1.28 (h^-1 M_sun)
    mass = 8.8e13 * np.power((0.66 * clusters['richness'] + 9.567)/20,1.28)/(cosmo.H0/100) # in M_sun/h
    mass = Column(name = 'mass', data=mass)
    clusters.add_column(mass)
    
    #add column with r200 in Mpc
    rho_200 = 200 * cosmo.critical_density(clusters['redshift']).to('Msun/Mpc3') #in Msun/Mpc3
    mass_200 = clusters['mass']*u.Msun # in Msun
    r_200 = np.power(mass_200/rho_200/(4*np.pi/3.), 1/3.).value #in Mpc
    clusters['r200_mpc'] = r_200

    # give unified name
    name = [give_name(ra,dec) for ra,dec in zip(clusters['ra'],clusters['dec'])]
    clusters['name'] = name

    return clusters



def unify_qso_catalog_xmq(qsos):
    """Unifies the name of columns that are relevant for most analyses"""
    qsos = Table(qsos)
    qsos.rename_column('Z','redshift')
    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('FUV','mag_fuv')
    qsos.add_column(Column(name='objid_xmq',data=np.arange(len(qsos))+1))
    # unify name
    name = [give_name(ra,dec) for ra,dec in zip(qsos['ra'],qsos['dec'])]
    qsos['name'] = name
    return qsos


def unify_qso_catalog_legacy(qsos):
    """Unifies the name of columns that are relevant for most analyses"""
    qsos = Table(qsos)
    # remove those with no medium resolution observation
    cond = qsos['FUV M'] ==  ' . . .'
    qsos = qsos[~cond]

    name_aux = [name.split('">')[1].split('</a')[0] for name in qsos['Target Name']]
    qsos['Redshift'] = np.where(qsos['Redshift'] == '  . .', -1, qsos['Redshift'])
    redshift_aux = [float(z) for z in qsos['Redshift']]

    qsos['redshift'] = redshift_aux
    qsos['NAME_ORIG'] = name_aux
    qsos['mag_fuv'] = -1*np.ones(len(qsos))  # need to get FUV photometry from somewhere

    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('Median S/N','sn_tot')

    # ID
    qsos['objid_legacy'] = np.arange(len(qsos))+1

    # unify name
    name = [give_name(ra,dec) for ra,dec in zip(qsos['ra'],qsos['dec'])]
    qsos['name'] = name

    return qsos

def unify_qso_catalog_uvqs(qsos):
    """Unifies the name of columns that are relevant for most analyses"""
    qsos = Table(qsos)
    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('FUV','mag_fuv')
    qsos.rename_column('Z','redshift')
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    # unify name
    name = [give_name(ra,dec) for ra,dec in zip(qsos['ra'],qsos['dec'])]
    qsos['name'] = name

    return qsos