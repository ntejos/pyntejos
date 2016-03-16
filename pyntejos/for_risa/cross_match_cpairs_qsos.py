import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits,ascii
from astropy.table import Table, Row, Column, vstack
from astropy.cosmology import Planck13 as cosmo
from pyntejos.utils import group_z, give_dz
import pyntejos.coord as coord
from cross_match_cpairs_qsos_utils import *    

"""This script is to cross-match QSOs with cluster-pairs, and find the
best QSO targets for HST based on a given merit. - Nicolas Tejos 2016

"""

#Global parameters are read from configuration file
from cross_match_cpairs_qsos_config import params 

#Read QSO catalog
dataset = params['dataset']
if dataset == 'gabor':
    qsos = '/media/ntejos/disk1/catalogs/UV_QSOs/cluster_fuv19_Gabor.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_gabor(qsos) 
if dataset == 'newz1':
    qsos = '/media/ntejos/disk1/catalogs/UV_QSOs/z1qso_new.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_newz1(qsos) 
if dataset == 'xmq':
    qsos = '/media/ntejos/disk1/catalogs/UV_QSOs/xmq_GALEX_cutz.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_xmq(qsos) 
if dataset == 'lehner':
    qsos = '/media/ntejos/disk2/projects/COS-Web/data/Lehner_20140828/lehner_database.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_lehner(qsos) 
if dataset == 'legacy':
    qsos = '/media/ntejos/disk2/data/COS/legacy-archive/samples/QSOALS_websample.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_legacy(qsos) 
if dataset == 'uvqs':
    qsos = '/media/ntejos/disk1/catalogs/UV_QSOs/hlsp_uvqs_multi_multi_all_multi_v1_redshifts.fits'
    qsos = fits.getdata(qsos)
    qsos = Table(qsos)
    qsos = unify_qso_catalog_uvqs(qsos)
if dataset == 'unobserved':
    # Join uvqs and xmq
    qso_datasets = ['xmq','uvqs']
    qsos = Table()
    for q_d in qso_datasets:
        if q_d == 'xmq':
            # filename = '/media/ntejos/disk1/catalogs/UV_QSOs/xmq_GALEX_cutz.fits'
            filename = './data/xmq_GALEX_cutz.fits'
            qsos_aux = fits.getdata(filename)   
            qsos_aux = unify_qso_catalog_xmq(qsos_aux)
        elif q_d == 'uvqs':
            # filename = '/media/ntejos/disk1/catalogs/UV_QSOs/hlsp_uvqs_multi_multi_all_multi_v1_redshifts.fits'
            filename = './data/hlsp_uvqs_multi_multi_all_multi_v1_redshifts.fits')
            qsos_aux = fits.getdata(filename)   
            qsos_aux = unify_qso_catalog_uvqs(qsos_aux)
        qsos_aux['dataset'] = q_d  # add origin dataset as column
        qsos = vstack([qsos,qsos_aux])
    #get rid of those w/o GALEX mag
    cond = qsos['mag_fuv']>0
    qsos = qsos[cond]

#keep only brighter than mag_fuv_min for those observed ones
if dataset not in ['lehner', 'legacy']:
    cond = qsos['mag_fuv'] < params['mag_fuv_min']
    qsos = qsos[cond]
if dataset not in ['legacy']:
    qsos.add_column(Column(name='flux_fuv',data=fuv_mag_to_flux(qsos['mag_fuv'])))

# make sure that J1410+2304 has the right FUV flux
if dataset == 'lehner':
    cond = qsos['name'] == 'J1410+2304'
    qsos['mag_fuv']=np.where(cond,18.7,qsos['mag_fuv'])

#remove duplicates or bad targets (i.e. bad wavelength, or bright star nearby)
duplicates = ['J1617+0854','J1301+5902','J1124+4201','J1105+3425',\
              'J0235-0402','J1358+0213','J1210+2725','J1017+4702'] #'J1456+2750'
for dup in duplicates:
    if dataset in ['lehner', 'legacy']:
        continue
    cond = qsos['name'] == dup
    qsos = qsos[~cond]


#Set up quality flag for cluster-pairs
"""
if (opt.richmin == 20) and (opt.max_tsep == 20):
    opt['quality'] = 1
elif (opt.richmin == 20) and (opt.max_tsep == 25):
    opt['quality'] = 2
elif (opt.richmin == 15) and (opt.max_tsep == 20):
    opt['quality'] = 3
elif (opt.richmin == 15) and (opt.max_tsep == 25):
    opt['quality'] = 4
elif (opt.richmin == 10) and (opt.max_tsep == 20):
    opt['quality'] = 5
elif (opt.richmin == 10) and (opt.max_tsep == 25):
    opt['quality'] = 6
else:
    opt['quality'] = 0
"""

# redefine ids
qsos['id'] = np.arange(len(qsos))+1

#Read cluster-pair catalog; will read from both SDSS and DES
# cl_datasets = ['sdss_dr8','des_sva1']
cl_datasets = ['des_yr1']
cpairs = Table()
for cl_d in cl_datasets:
    cpairs_aux = fits.getdata('cpairs_richmin{}_l{}_{}.fits'.format(int(params['richmin']),int(params['max_tsep']), cl_d))
    cpairs_aux = Table(cpairs_aux)
    cpairs_aux['dataset'] = cl_d    
    cpairs = vstack([cpairs,cpairs_aux])   

#add objid to cpairs
cpairs.add_column(Column(name='objid',data=np.arange(len(cpairs))+1))
#add comoving line of sight distance for transverse calculation in mpc
losd = cosmo.comoving_transverse_distance(cpairs['redshift']).value
cpairs.add_column(Column(name='losd_mpc',data=losd))


#For a given QSO, store the ids of cluster-pairs within d Mpcs
ids = []
y_mpc = [] #impact parameters
x_mpc = [] #distance to closest cluster along axis

for qso in qsos:
    ids_aux = []
    y_mpc_aux = []
    x_mpc_aux = []
    
    #QSO coordinates
    ra0  = qso['ra']
    dec0 = qso['dec']
    zqso = qso['redshift']

    #Define the closest subsample of cpairs that could have an max_dsep
    #1. get angular separations to each cluster of cpairs
    sep_ang_deg1 = coord.ang_sep(cpairs['ra1'],cpairs['dec1'],ra0,dec0)
    sep_ang_deg2 = coord.ang_sep(cpairs['ra2'],cpairs['dec2'],ra0,dec0)
    sep_ang_deg_max  = np.where(sep_ang_deg1>=sep_ang_deg2,sep_ang_deg1,sep_ang_deg2) #max between the 2
    sep_ang_rad_max = coord.RAD_PER_DEG * sep_ang_deg_max
    
    #2. get comoving transverse separations
    tsep_mpc_max = cpairs['losd_mpc'] * sep_ang_rad_max
    
    #3. Only keep those with cpairs['sep_mpc'] + max_dsep <= tsep_mpc_max  
    cond = tsep_mpc_max <= cpairs['sep_mpc'] + params['max_dsep']
    cpairs_aux = cpairs[cond]

    #get rid of cpairs behind the QSO too
    cond = cpairs_aux['redshift'] <= zqso - give_dz(params['vel_prox'],zqso)
    cpairs_aux = cpairs_aux[cond]

    for i,cpair in enumerate(cpairs_aux):
        min_sep, sep_x_mpc = get_impact_parameters_cpair(cpair,ra0,dec0,cosmo,npoints=20)
        if sep_x_mpc < params['min_xsep']:
            continue 
        if min_sep < params['max_dsep']:
            ids_aux += [cpair['objid']]
            y_mpc_aux += [min_sep]
            x_mpc_aux += [sep_x_mpc]
    print '{} cluster-pairs found for qso {}, {}'.format(len(ids_aux),qso['id'], qso['name'])
    ids += [tuple(ids_aux)]
    y_mpc += [y_mpc_aux]
    x_mpc += [x_mpc_aux]

#append indices of cluster-pairs and relevant info
qsos.add_column(Column(name='cpairs_ids',data=ids))
qsos.add_column(Column(name='cpairs_y_mpc',data=y_mpc))
qsos.add_column(Column(name='cpairs_x_mpc',data=x_mpc))


#number of total cpairs within contraints
n_tot = np.array([len(a) for a in ids])
qsos.add_column(Column(name='n_tot',data=n_tot))

#get number of independent ones
n_ind = []
for idsi in ids:
    if len(idsi)==0:
        n_ind += [0]
    else:
        z = cpairs[list(idsi)]['redshift']
        grouped_z = group_z(z,params['dv_group'])
        n_ind += [np.max(grouped_z)+1]
qsos.add_column(Column(name='n_ind',data=n_ind))       

#get rid of those not having cluster-pairs
cond = qsos['n_ind'] > 0
qsos = qsos[cond]

#create figure of merit: 
if dataset not in ['lehner', 'legacy']:
    add_zmax_pairs(qsos,cpairs) #maximum redshift for intervening cpairs
    bonus = np.where(qsos['zmax_cpairs']<0.15,2,1) # if z<0.15 only G130M needed to cover both HI and OVI, so a bonus is applied
    qsos.add_column(Column(name='bonus',data=bonus))
    merit = qsos['n_ind'] * (qsos['flux_fuv']) * qsos['bonus']

    #estimate SN for 1 HST orbit at 
    add_wa_HI_OVI_and_S2N(qsos, cpairs)

    n_orb_g130 = -1*np.ones(len(qsos))
    n_orb_g160 = -1*np.ones(len(qsos))
    sn_both = [np.min([sn_HI,sn_OVI]) for sn_HI,sn_OVI in zip(qsos['sn_HI'],qsos['sn_OVI'])]
    sn_both = np.array(sn_both)
    #if bonus = 2 only g130m needed
    n_orb_g130 = np.where(qsos['bonus']==2,(params['sn_min']/sn_both)**2,n_orb_g130)
    n_orb_g160 = np.where(qsos['bonus']==2,0,n_orb_g160)
    #if bonus == 1 both gratings are needed
    n_orb_g130 = np.where(qsos['bonus']==1,(params['sn_min']/qsos['sn_OVI'])**2,n_orb_g130)
    n_orb_g160 = np.where(qsos['bonus']==1,(params['sn_min']/qsos['sn_HI'])**2,n_orb_g160)
    #approximate up to integer
    n_orb_g130 = np.where(n_orb_g130>0,np.array([int(v)+1 for v in n_orb_g130]),n_orb_g130)
    n_orb_g160 = np.where(n_orb_g160>0,np.array([int(v)+1 for v in n_orb_g160]),n_orb_g160)
    #total number of orbits
    n_orb_tot = n_orb_g130 + n_orb_g160
    n_orb_tot = np.where(qsos['bonus']==2,n_orb_g130,n_orb_tot)
    #store this new info
    qsos.add_column(Column(name='n_orb_g130m',data=n_orb_g130))
    qsos.add_column(Column(name='n_orb_g160m',data=n_orb_g160))
    qsos.add_column(Column(name='n_orb_tot',data=n_orb_tot))
    
    #New figure of merit
    merit = qsos['n_ind'] / qsos['n_orb_tot']
    qsos.add_column(Column(name='merit',data=merit))

else:
    merit = qsos['sn_tot']
    dumb_data = -1*np.ones(len(qsos))
    qsos.add_column(Column(name='n_orb_g130m',data=dumb_data))
    qsos.add_column(Column(name='n_orb_g160m',data=-dumb_data))
    qsos.add_column(Column(name='n_orb_tot',data=dumb_data))
    qsos.add_column(Column(name='merit',data=merit))

#remove QSOS with merit <= something to speed up
if dataset not in ['lehner', 'legacy']:
    cond = qsos['merit'] > 0 
    qsos = qsos[cond]
    
#sort by merit
qsos.sort(['merit'])
qsos.reverse()

#Estimate orbits for survey
n_ind_sum = [0]
n_orb_sum = [0]
for qso in qsos:
    n_ind_sum += [n_ind_sum[-1] + qso['n_ind']]
    n_orb_sum += [n_orb_sum[-1] + qso['n_orb_tot']]
#add to qsos
qsos.add_column(Column(name='n_ind_sum',data=n_ind_sum[1:]))
qsos.add_column(Column(name='n_orb_sum',data=n_orb_sum[1:]))

#write html table
write_html_table(qsos, cpairs, params)
#write latex table
write_latex_table(qsos,cpairs, max_n_orb=30)

#store fits table
use_cols = ['name','ra','dec','mag_fuv','redshift','n_tot','n_ind']
ascii.write(qsos, output='qsos.txt',include_names=use_cols)

# figures
# make_Nmin_figure()
make_orbits_plot(qsos, params)
