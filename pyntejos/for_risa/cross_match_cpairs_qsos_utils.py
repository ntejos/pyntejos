import numpy as np
import matplotlib.pyplot as pl
from astropy.table import Table, Row, Column, vstack
from astropy import units as u
from astropy.io import fits, ascii
import pyntejos.coord as coord
from pyntejos.utils import give_dz, give_dv
from pyntejos.catalogs import give_name
from pyntejos.catalogs import unify_cluster_catalog_redmapper, \
                                unify_qso_catalog_xmq,\
                                unify_qso_catalog_legacy,\
                                unify_qso_catalog_uvqs

# Colors for plotting
blue = '#007996'
green = '#668D3C'
orange = '#FF8642'
red = '#C0362C'
everything = '#816C5B'
yellow = '#FF8642'
black = 'k'

def fuv_mag_to_flux(mag_fuv):
    """Returns in units of erg/sec/cm2/A/ 1.4* 10**(-15)"""
    return np.power(10**(mag_fuv - 18.82),-1/2.5) #* 1.4* 10**(-15) 

def unify_qso_catalog_gabor(qsos):
    """Unifies the name of columns that are relevant for the analysis"""
    qsos.rename_column('z','redshift')
    return qsos

def unify_qso_catalog_newz1(qsos):
    """Unifies the name of columns that are relevant for the analysis"""
    qsos.rename_column('Z','redshift')
    qsos.rename_column('RAD','ra')
    qsos.rename_column('DECD','dec')
    qsos.rename_column('FUV','mag_fuv')
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    return qsos

def unify_qso_catalog_xmq_old(qsos):
    """Unifies the name of columns that are relevant for the analysis"""
    qsos.rename_column('Z','redshift')
    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('FUV','mag_fuv')
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    #get rid of those w/o GALEX mag
    cond = qsos['mag_fuv']>0
    qsos = qsos[cond]
    return qsos

def unify_qso_catalog_lehner(qsos):
    """Unifies the name of columns that are relevant for the analysis"""
    qsos.rename_column('name','NAME_OLD')
    qsos['mag_fuv'] = np.where(qsos['mag_fuv']<=0,19,qsos['mag_fuv'])
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    return qsos
    
def unify_qso_catalog_legacy_old(qsos):
    """Unifies the name of columns that are relevant for the analysis"""
    
    # remove those with no medium resolution observation
    cond = qsos['FUV M'] ==  ' . . .'
    qsos = qsos[~cond]

    name_aux = [name.split('">')[1].split('</a')[0] for name in qsos['Target Name']]
    qsos['Redshift'] = np.where(qsos['Redshift'] == '  . .', -1, qsos['Redshift'])
    redshift_aux = [float(z) for z in qsos['Redshift']]
    
    dummy = -1*np.ones(len(qsos))
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    qsos.add_column(Column(name='redshift',data=redshift_aux))
    qsos.add_column(Column(name='NAME_OLD',data=name_aux))
    qsos.add_column(Column(name='mag_fuv',data=dummy))

    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('Median S/N','sn_tot')
    return qsos

def unify_qso_catalog_uvqs_old(qsos):
    """Unifies the name of columns that are relevant for the analysis"""

    qsos.rename_column('RA','ra')
    qsos.rename_column('DEC','dec')
    qsos.rename_column('FUV','mag_fuv')
    qsos.rename_column('Z','redshift')
    qsos.add_column(Column(name='id',data=np.arange(len(qsos))+1))
    return qsos

def add_impact_parameters_cpairs(cpairs,ra0,dec0, cosmo):
    """Adds a column with impact parameter between the straight line given
    by a cluster-pair and (ra0,dec0) in Mpc ('d'), as well as the
    distance along the cpair axis to the closest cluster of the pair in
    Mpc ('x').

    """
    
    d = []
    x = []
    for cpair in cpairs:
        d_aux, x_aux = get_impact_parameters_cpair(cpair,ra0,dec0,cosmo,npoints=1000)
        d = d + [d_aux]
        x = x + [x_aux]
    d = np.array(d)
    x = np.array(x)
    d = Column(name='d_mpc',data=d)
    x = Column(name='x_mpc',data=x)
    cpairs.add_columns([d,x])


def get_line_radec_sepdv(ra1,dec1,z1,ra2,dec2,z2,ra0,dec0,zmed,cosmo,npoints=100):
    """Defines a line in terms of (ra,dec) offset around (ra0,dec0) and
    dv around zmed, using a linear sampling. Output:
    ra_mpc,dec_mpc,sep_mpc,dv"""
    
    ra  = np.linspace(ra1,ra2,npoints,endpoint=True)
    dec = np.linspace(dec1,dec2,npoints,endpoint=True)
    z   = np.linspace(z1,z2,npoints,endpoint=True)
    
    #angular separation to the (ra0,dec0) in degrees
    ra_ang_deg  = (ra - ra0) * np.cos(dec0 * coord.RAD_PER_DEG)
    dec_ang_deg = (dec - dec0)
    sep_ang_deg = coord.ang_sep(ra,dec,ra0,dec0)
    
    #angular separation to the (ra0,dec0) in radians
    ra_ang_rad  = coord.RAD_PER_DEG * ra_ang_deg 
    dec_ang_rad = coord.RAD_PER_DEG * dec_ang_deg
    sep_ang_rad = coord.RAD_PER_DEG * sep_ang_deg
    
    #relative co-moving distances to the (ra0,dec0) in Mpc
    td = cosmo.comoving_transverse_distance(z) # in Mpc
    ra_mpc  = td * ra_ang_rad
    dec_mpc = td * dec_ang_rad
    sep_mpc = td * sep_ang_rad
    
    #relative velocity at zmed
    dv = give_dv(z,zmed)

    return ra_mpc,dec_mpc,sep_mpc,dv

def get_impact_parameters_cpair(cpair,ra0,dec0,cosmo,npoints=1000):
    """Returns the impact parameter between the straight line given by a
    cluster-pair and (ra0,dec0), as well as the distance along the
    cpair axis to the closest cluster of the pair.

    """
    
    ra1  = cpair['ra1']
    dec1 = cpair['dec1']
    z1   = cpair['z1']
    ra2  = cpair['ra2']
    dec2 = cpair['dec2']
    z2   = cpair['z2']
    zmed = cpair['redshift']
    
    ra_mpc,dec_mpc,sep_mpc,dv = get_line_radec_sepdv(ra1,dec1,z1,ra2,dec2,z2,ra0,dec0,zmed,cosmo,npoints=npoints)
    
    #get hypotenuse
    sep = np.hypot(ra_mpc,dec_mpc)

    #get minimum distance
    sep_min = np.min(sep)
    
    #lets add the distance along the cpair axis to the closest cluster
    #find the closest cluster distance first
    sep_cl_mpc1 = np.hypot(ra_mpc[0],dec_mpc[0])
    sep_cl_mpc2 = np.hypot(ra_mpc[-1],dec_mpc[-1])
    sep_cl_mpc = np.min([sep_cl_mpc1.value,sep_cl_mpc2.value])
    #then, project along the cluster axis 
    sep_x_mpc = np.sqrt(sep_cl_mpc**2 - sep_min.value**2)
    
    return sep_min.value, sep_x_mpc

def Nmin(e,dz,s,a,a_err):
    e = np.array(e)
    dz = np.array(dz)
    s = np.array(s)
    a = np.array(a)
    a_err = np.array(a_err)
    
    return (e/dz/a)*(s**2)/((e-1.)-s*a_err/a)**2

def add_wa_HI_OVI_and_S2N(qsos,cpairs,verbose=True):
    """It adds expected wavelength for HI and OVI at the cluster-pair
    redshifts, and also adds an estimation of the signal-to-noise at
    those wavelengths in 1 HST orbit exposure time.

    """
    t_orbit = 2700 # in seconds
    
    wa_HI = []
    wa_OVI = []
    sn_HI = []
    sn_OVI = []

    #estimate SN per 1 orbit
    
    for qso in qsos:
        if verbose:
            print 'Estimating SN for qso {}'.format(qso['id'])
        
        wa_HI_aux = []
        wa_OVI_aux = []
        sn_HI_aux = []
        sn_OVI_aux = []
        
        #estimate SN per 1 orbit
        cos_etc = s2n(t_orbit,qso['mag_fuv'])
        
        for cpair_id in list(qso['cpairs_ids']):
            cpair = cpairs[cpair_id-1] #obijds start from 1 (not 0)
            z_pair = cpair['redshift']
            wa0_HI = 1216 * (1+z_pair)
            wa0_OVI = 1031 * (1+z_pair)
            wa_HI_aux  += [wa0_HI]
            wa_OVI_aux += [wa0_OVI]
            #SN now
            cond = np.fabs(cos_etc['wavelength'] -wa0_HI) == np.min(np.fabs(cos_etc['wavelength'] -wa0_HI))
            sn_HI_aux += [cos_etc['sn'][cond][0]]
            cond = np.fabs(cos_etc['wavelength'] -wa0_OVI) == np.min(np.fabs(cos_etc['wavelength'] -wa0_OVI))
            sn_OVI_aux += [cos_etc['sn'][cond][0]]
            
        
        wa_HI += [wa_HI_aux]
        wa_OVI += [wa_OVI_aux]
        sn_HI += [np.min(sn_HI_aux)] # keep the minimum of the list
        sn_OVI += [np.min(sn_OVI_aux)] # keep the minimum of the list
    
    qsos.add_column(Column(name='wa_HI',data=wa_HI))
    qsos.add_column(Column(name='wa_OVI',data=wa_OVI))
    qsos.add_column(Column(name='sn_HI',data=sn_HI))
    qsos.add_column(Column(name='sn_OVI',data=sn_OVI))
        
def add_zmax_pairs(qsos,cpairs):
    zmax_cpairs = []
    for qso in qsos:
        cpairs_ids = list(qso['cpairs_ids'])
        if len(cpairs_ids)==0:
            zmax_aux = 99
        else:
            zmax_aux = np.max(cpairs[cpairs_ids]['redshift'])
        zmax_cpairs += [zmax_aux]
    qsos.add_column(Column(name='zmax_cpairs',data=zmax_cpairs))


def s2n(t,FUV,tnorm=1000):
    # load cos etc simulations for flat spectrum 1000s exp, and FUV = 17 mag
    cos_g130m = ascii.read('cos_etc_g130m_v24.1.csv')
    cos_g160m = ascii.read('cos_etc_g160m_v24.1.csv')
    #separate them at 1400 A
    cond = cos_g130m['wavelength'] < 1405
    cos_g130m = cos_g130m[cond]
    cond = cos_g160m['wavelength'] >= 1405
    cos_g160m = cos_g160m[cond]
    #merge both
    cos = vstack([cos_g130m,cos_g160m], join_type='exact')
    
    # Signal
    signal = cos['target_counts']*t/tnorm * 10.**((FUV-17.)/(-2.5))
    
    # Noise terms
    dark = cos['dark_counts']*t/tnorm
    sky = cos['sky_counts']*t/tnorm
        
    # Noise
    var = signal + dark + sky
    sig = np.sqrt(var)
    
    #append S/N to cos
    sn = signal/sig * np.sqrt(6) # per-resolution element of 3 pixels
    cos.add_column(Column(name='sn', data=sn))
    return cos
    
    

def write_html_table(qsos,cpairs, params):
    """Writes a master html table from the qsos catalog"""
    
    filename = 'qsos_{}_richmin{}_l{}_sn{}_y{}.html'.format(params['dataset'],
        params['richmin'],params['max_tsep'],params['sn_min'], params['max_dsep'])
    f  = open(filename,'w')
    s = """ <html>
    <head>
    <meta charset="utf-8"/>
    <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
    </head>
    <body>
    <table border="1">\n"""
    f.write(s)
    print filename
    #headers    
    s = """<tr><th>QSO</th><th>RA</th><th>DEC</th><th>FUV</th><th>Z</th><th>N_TOT</th><th>N_IND</th><th>N_SUM</th><th>MERIT</th><th>N_ORB_SUM</th><th>CLUSTERS</th><th>Z_PAIRS</th><th>X_MPC</th><th>Y_MPC</th><th>L_MPC</th><th>WA_HI</th><th>WA_OVI</th></tr>\n"""

    f.write(s)
    for qso in qsos:
        name = qso['name']
        ra = '{:.5f}'.format(qso['ra'])
        dec = '{:.5f}'.format(qso['dec'])
        z = '{:.4f}'.format(qso['redshift'])
        fuv = '{:.2f}'.format(qso['mag_fuv'])
        n_tot = qso['n_tot']
        n_ind = qso['n_ind']
        n_ind_sum = qso['n_ind_sum']
        n_orb_sum = qso['n_orb_sum']
        merit = '{:.1f}'.format(qso['merit'])
        qso_url = 'http://skyserver.sdss.org/dr10/en/tools/chart/navi.aspx?ra={}&dec={}&opt='.format(ra,dec)
        qso_chart = '<a href="{}" target="_blank">{}</a>'.format(qso_url,name)
        s_qso = '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td> <td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td></tr> \n'.format(qso_chart,ra,dec,fuv,z,n_tot,n_ind,n_ind_sum,merit,n_orb_sum)
        f.write(s_qso)
        
        for i,cpair_id in enumerate(list(qso['cpairs_ids'])):
            cpair = cpairs[cpair_id-1] #'objid starts from 1 (not 0)'
            ra1,dec1 = cpair['ra1'],cpair['dec1']
            ra2,dec2 = cpair['ra2'],cpair['dec2']
            z_pair = '{:.4f}'.format(cpair['redshift'])
            wa_HI  = int(1216 * (1+cpair['redshift']))
            wa_OVI = int(1031 * (1+cpair['redshift']))
            y_mpc  = '{:.1f}'.format(qso['cpairs_y_mpc'][i])
            x_mpc  = '{:.1f}'.format(qso['cpairs_x_mpc'][i])
            l_mpc  = '{:.1f}'.format(cpair['sep_mpc'])
            
            s1_aux = 'http://skyserver.sdss.org/dr10/en/tools/chart/navi.aspx?ra={}&dec={}&opt='.format(ra1,dec1)
            s2_aux = 'http://skyserver.sdss.org/dr10/en/tools/chart/navi.aspx?ra={}&dec={}&opt='.format(ra2,dec2)
            
            s_pair_i_chart = '<a href="{}" target="_blank"> P{}-CL1</a><a href="{}" target="_blank"> P{}-CL2</a>'.format(s1_aux,i+1,s2_aux,i+1)
            s_pair_i = '<tr><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td> <td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr> \n'.format(s_pair_i_chart,z_pair,x_mpc,y_mpc,l_mpc,wa_HI,wa_OVI)
            f.write(s_pair_i)
            
    
    s = """ </table>\n</body>\n</html>\n"""
    f.write(s)
    f.close()


def write_latex_table(qsos,cpairs,max_n_orb=30):
    """Writes latex table for HST proposal"""

    header = """
    \\begin{table} 
    \scriptsize 
    \centering 
    \\begin{minipage}{170mm}
    
    \caption{ A summary of our proposed QSO observations, in terms of
    known cluster-pairs intersected. In combination with archival data,
    this sample will allow a $99\%$ statistical detection of the WHIM
    in these inter-cluster filaments.}\label{tab:1}
    
    \\begin{tabular}{@{}cccccccccccc@{}}
    \hline
    \hline
    \multicolumn{6}{c}{Targeted QSOs}  & & \multicolumn{5}{c}{Targeted Filaments} \\\ \n
    \cline{1-6}  \cline{8-12} 
    QSO Name & RA  & Dec  & $z_{\\rm QSO}$ & FUV & n$_{\\rm orb}$ &  &$z_{\\rm flmnt}$ & Impact Par. & Both      & HI Lya & OVI  \\\ \n
             & (degrees) & (degrees)    &   & (AB) &  &           &      & (Mpc)       & spec-$z$? & (\AA)  & (\AA)  \\\ \n
    
    \hline
    """
    
    filename = '/media/ntejos/disk1/proposals/HST/COS/c24/filaments/tab1.tex'
    f  = open(filename,'w')
    f.write(header)

    cond = qsos['n_orb_sum'] < max_n_orb
    qsos_aux = qsos[cond]
    
    for qso in qsos_aux:
        name = qso['name']
        ra = '{:.5f}'.format(qso['ra'])
        dec = '{:.5f}'.format(qso['dec'])
        z = '{:.4f}'.format(qso['redshift'])
        fuv = '{:.2f}'.format(qso['mag_fuv'])
        n_orb_sum = '{}+{}'.format(int(qso['n_orb_g130m']),int(qso['n_orb_g160m']))
        #s_qso = '{}&{}&{}&{}&{}&{}& & & & & & \\\ \n'.format(name,ra,dec,z,fuv,n_orb_sum)
        #f.write(s_qso)
        s_qso = '{}&{}&{}&{}&{}&{}'.format(name,ra,dec,z,fuv,n_orb_sum)
        
        
        for i,cpair_id in enumerate(list(qso['cpairs_ids'])):
            cpair = cpairs[cpair_id-1] #'objid starts from 1 (not 0)'
            ra1,dec1 = cpair['ra1'],cpair['dec1']
            ra2,dec2 = cpair['ra2'],cpair['dec2']
            z_pair = '{:.4f}'.format(cpair['redshift'])
            wa_HI  = int(1216 * (1+cpair['redshift']))
            wa_OVI = int(1031 * (1+cpair['redshift']))
            y_mpc  = '{:.1f}'.format(qso['cpairs_y_mpc'][i])
            x_mpc  = '{:.1f}'.format(qso['cpairs_x_mpc'][i])
            l_mpc  = '{:.1f}'.format(cpair['sep_mpc'])
            if cpair['quality']==1:
                both_spec = 'y'
            else:
                both_spec = 'n'
            
            #s_pair = ' & & & & & & & {}&{}&{}&{}&{}\\\ \n'.format(z_pair,y_mpc,both_spec,wa_HI,wa_OVI)
            #f.write(s_pair)
            s_pair = '&&{}&{}&{}&{}&{}\\\ \n'.format(z_pair,y_mpc,both_spec,wa_HI,wa_OVI)
            
        f.write(s_qso+s_pair)
    tail = """      \hline
    \end{tabular}
  \end{minipage}
\end{table}"""
    
    f.write(tail)
    f.close()


def make_Nmin_figure():
    s = np.arange(0,5,0.1)
    Nmin_BLA = Nmin(6,0.008,s,20,5)
    Nmin_BLA_m = Nmin(6-4,0.008,s,20,5)
    Nmin_BLA_p = Nmin(6+4,0.008,s,20,5)
    
    Nmin_HI  = Nmin(2,0.008,s,130,0)
    Nmin_HI_m  = Nmin(2-0.5,0.008,s,130,0)
    Nmin_HI_p  = Nmin(2+0.5,0.008,s,130,0)
    
    Nmin_OVI = Nmin(4,0.008,s,20,1)
    Nmin_OVI_m = Nmin(4-3,0.008,s,20,0)
    Nmin_OVI_p = Nmin(4+4,0.008,s,20,0)
    
    fig = pl.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(s,Nmin_HI,c=blue,lw=2,label='Total HI',ls='-')
    #ax.fill_between(s,Nmin_HI_p,Nmin_HI_m,color=blue,alpha=0.5)
    ax.plot(s,Nmin_BLA,c=red,lw=2,label='BLAs',ls='--')
    #ax.fill_between(s,Nmin_BLA_p,Nmin_BLA_m,color=red,alpha=0.5)
    ax.plot(s,Nmin_OVI,c=green,lw=2,label='OVI',ls=':')
    ax.minorticks_on()
    ax.set_xlabel('Significance level ($\sigma$)',fontsize='x-large')
    ax.set_ylabel('Number of QSO-filament pairs',fontsize='x-large')
    ax.legend(numpoints=1)
    ax.hlines(20,s[0],s[-1],linestyles='dashed',alpha=0.5)
    ax.annotate('In the archive',(0.3,7),fontsize='large',va='bottom')
    ax.fill_between(np.append(s,12),0,8,color=yellow, alpha=0.3)
    ax.set_xlim(0,5)
    ax.set_ylim(0,50)
    ax.annotate('Optimal number',(0.3,20),fontsize='large',va='bottom')
    fig.savefig('Nmin.png',dpi=300)
    pl.show()

def make_orbits_plot(qsos, dataset):
    fig = pl.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    
    ymax = 20
    xmax = 50

    ax.plot(qsos['n_orb_sum'],qsos['n_ind_sum'],'ko-')
    ax.set_xlabel('Number of HST orbits',fontsize='x-large')
    ax.set_ylabel('Number of QSO-filament pairs',fontsize='x-large')
    
    #ax.set_title('Richness={}; L={}; dataset={}; SN={}'.format(opt.richmin,opt.max_tsep,opt.dataset,opt.sn_min))
    ax.minorticks_on()
    ax.hlines(12,0,xmax,linestyles='dashed',alpha=0.5)
    ax.annotate('Our request',(2,12),fontsize='large',va='bottom')
    
    ax.set_ylim(0,ymax)
    ax.set_xlim(0,xmax)
    if params['dataset']=='lehner':
        ax.set_xlim(0,250)
    fig.savefig('Norb_{}.png'.format(dataset),dpi=300)
    pl.show()
    
    
