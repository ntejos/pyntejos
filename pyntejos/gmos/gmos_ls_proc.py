#! /usr/bin/env python
#   (c) 2016-Jun-08  shaw@noao.edu
# Modified by ntejos to make it a more flexible script

import sys
import os
import copy
from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos, onedspec
import fileSelect as fs
from utils import ask_user

def create_database():
    # check whether
    pass

# def gmos_ls_proc():
if 0:
    '''
    GMOS Data Reduction Cookbook companion script to the chapter:
       "Reduction of Longslit Spectra with PyRAF"

    PyRAF script to:
    Process GMOS spectra for AM2306-721, in program GS-2007A-Q-76.

    The names for the relevant header keywords and their expected values are
    described in the DRC chapter entitled "Supplementary Material"

    Perform the following starting in the parent work directory:
        cd /path/to/work_directory

    Place the fileSelect.py module in your work directory. Now execute this
    script from the unix prompt:
        python gmos_ls_proc.py
    '''

    print ("### Begin Processing GMOS/Longslit Images ###")
    print ("###")
    print ("=== Creating MasterCals ===")

    # This whole example depends upon first having built an sqlite3 database of metadata:
    #    cd ./raw
    #    python obslog.py obsLog.sqlite3
    dbFile='./raw/obsLog.sqlite3'

    # From the work_directory:
    # Create the query dictionary of essential parameter=value pairs.
    # Select bias exposures within ~2 months of the target observations:
    qd = {'Full':{'use_me':1,
           'Instrument':'GMOS-S','CcdBin':'2 4','RoI':'Full',
           'Disperser':'B600+_%','CentWave':485.0,'AperMask':'1.0arcsec',
           'Object':'AM2306-72%',
           'DateObs':'2007-06-05:2007-07-07'}
          }
    # Make another copy for the CenterSpec RoI:
    qd['CenSp'] = copy.deepcopy(qd['Full'])
    qd['CenSp'].update({'RoI':'CentSp','Object':'LTT9239'})

    print (" --Creating Bias MasterCal--")

    # Set the task parameters.
    gemtools.gemextn.unlearn()    # Disarm a bug in gbias
    gmos.gbias.unlearn()
    biasFlags = {
        'logfile':'biasLog.txt','rawpath':'./raw/','fl_vardq':'yes',
        'verbose':'yes'
    }
    regions = ['Full','CenSp']
    for r in regions:
        # The following SQL generates the list of full-frame files to process.
        SQL = fs.createQuery('bias', qd[r])
        biasFiles = fs.fileListQuery(dbFile, SQL, qd[r])

        # The str.join() funciton is needed to transform a python list into a 
        # comma-separated string of file names that IRAF can understand.
        if len(biasFiles) > 1:
            gmos.gbias(','.join(str(x) for x in biasFiles), 'MCbias'+r, 
                       **biasFlags)


    # Clean up
    iraf.imdel("gS2007*.fits")

    print (" -- Creating GCAL Spectral Flat-Field MasterCals --")
    # Set the task parameters.
    qd['Full'].update({'DateObs':'*'})
    qd['CenSp'].update({'DateObs':'*'})
    gmos.gireduce.unlearn()
    gmos.gsflat.unlearn()
    # Normalize the spectral flats per CCD.
    # The response fitting should be done interactively.
    flatFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_dark':'no',
        'fl_fixpix':'no','fl_oversize':'no','fl_vardq':'yes','fl_fulldq':'yes',
        'rawpath':'./raw','fl_inter':'no','fl_detec':'yes',
        'function':'spline3','order':'13,11,28',
        'logfile':'gsflatLog.txt','verbose':'yes'
        }
    for r in regions:
        qr = qd[r]
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qr), qr)
        if len(flatFiles) > 0:
            gmos.gsflat (','.join(str(x) for x in flatFiles), 'MCflat'+r,
                    bias='MCbias'+r, **flatFlags)


    iraf.imdel('gS2007*.fits,gsS2007*.fits')

    print ("=== Processing Science Files ===")
    print (" -- Performing Basic Processing --")

    # Set task parameters.
    gmos.gsreduce.unlearn()
    sciFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_gscrrej':'no',
        'fl_dark':'no','fl_flat':'yes','fl_gmosaic':'yes','fl_fixpix':'no',
        'fl_gsappwave':'yes','fl_oversize':'no',
        'fl_vardq':'yes','fl_fulldq':'yes','rawpath':'./raw',
        'fl_inter':'no','logfile':'gsreduceLog.txt','verbose':'no'
    }
    arcFlags = copy.deepcopy(sciFlags)
    arcFlags.update({'fl_flat':'no','fl_vardq':'no','fl_fulldq':'no'})
    stdFlags = copy.deepcopy(sciFlags)
    stdFlags.update({'fl_fixpix':'yes','fl_vardq':'no','fl_fulldq':'no'})

    # Perform basic reductions on all exposures for science targets.
    print ("  - Arc exposures -")
    for r in regions:
        qr = qd[r]
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qr), qr)
        if len(arcFiles) > 0:
            gmos.gsreduce (','.join(str(x) for x in arcFiles), 
                           bias='MCbias'+r, **arcFlags)

    print ("  - Std star exposures -")
    r = 'CenSp'
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qd[r]), qd[r])
    if len(stdFiles) > 0:
        gmos.gsreduce (','.join(str(x) for x in stdFiles), bias='MCbias'+r,
                  flatim='MCflat'+r, **stdFlags)

    print ("  - Science exposures -")
    r = 'Full'
    sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qd[r]), qd[r])
    if len(sciFiles) > 0:
        gmos.gsreduce (','.join(str(x) for x in sciFiles), bias='MCbias'+r,
                  flatim='MCflat'+r, **sciFlags)

    # Clean up
    iraf.imdel('gS2007*.fits')
    
    print (" -- Determine wavelength calibration --")
    # Set task parameters
    gmos.gswavelength.unlearn()
    waveFlags = {
        'coordlist':'gmos$data/CuAr_GMOS.dat','fwidth':6,'nsum':50,
        'function':'chebyshev','order':5,
        'fl_inter':'no','logfile':'gswaveLog.txt','verbose':'no'
        }
    # The fit to the dispersion relation should be performed interactively.
    # Here we will use a previously determined result.
    # Need to select specific wavecals to match science exposures.
    prefix = 'gsS20070623S0'
    for arc in ['071', '081', '091', '109']:
        gmos.gswavelength (prefix+arc, **waveFlags)

    ### End of basic processing. Continue with advanced processing.

    print (" -- Performing Advanced Processing --")
    print (" -- Combine exposures, apply dispersion, subtract sky --")
    # Set task parameters.
    gemtools.gemcombine.unlearn()
    sciCombFlags = {
        'combine':'average','reject':'ccdclip',
        'fl_vardq':'yes','fl_dqprop':'yes',
        'logfile':'gemcombineLog.txt.txt','verbose':'no'
    }
    stdCombFlags = copy.deepcopy(sciCombFlags)
    stdCombFlags.update({'fl_vardq':'no','fl_dqprop':'no'})
    gmos.gstransform.unlearn()
    transFlags = {
        'fl_vardq':'yes','interptype':'linear','fl_flux':'yes',
        'logfile':'gstransLog.txt'
    }
    # The sky regions should be selected with care, using e.g. prows/pcols:
    #   pcols ("tAM2306b.fits[SCI]", 1100, 2040, wy1=40, wy2=320)
    gmos.gsskysub.unlearn()
    skyFlags = {
        'fl_oversize':'no','fl_vardq':'yes','logfile':'gsskysubLog.txt'
    }

    # Process the Standard Star
    prefix = "gs"
    qs = qd['CenSp']
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qs), qs)
    gemtools.gemcombine (','.join(prefix+str(x) for x in stdFiles), 
                         'LTT9239', **stdCombFlags)
    gmos.gstransform ('LTT9239', wavtraname='gsS20070623S0109', **transFlags)
    gmos.gsskysub ('tLTT9239', long_sample='20:70,190:230')

    print (" -- Extract Std spectrum --")
    # Extract the std spectruma using a large aperture.
    # It's important to trace the spectra interactively.
    gmos.gsextract.unlearn()
    extrFlags = {
        'apwidth':3.,'fl_inter':'no','find':'yes',
        'trace':'yes','tfunction':'chebyshev','torder':'6','tnsum':20,
        'background':'fit','bfunction':'chebyshev','border':2,
        'fl_vardq':'no','logfile':'gsextrLog.txt'
    }
    gmos.gsextract ("stLTT9239", **extrFlags)

    print (" -- Derive the Flux calibration --")
    gmos.gsstandard.unlearn()
    sensFlags = {
        'fl_inter':'yes','starname':'l9239','caldir':'onedstds$ctionewcal/',
        'observatory':'Gemini-South','extinction':'onedstds$ctioextinct.dat',
        'function':'chebyshev','order':9,'verbose':'no','logfile':'gsstdLog.txt'
    }
    gmos.gsstandard ('estLTT9239', sfile='std.txt', sfunction='sens', **sensFlags)

    # Process the science targets.
    # Use a dictionary to associate science targets with Arcs and sky regions.
    sciTargets = {
        'AM2306-721_a':{'arc':'gsS20070623S0071','sky':'520:720'}, 
        'AM2306-72_b':{'arc':'gsS20070623S0081','sky':'670:760,920:1020'}, 
        'AM2306-721_c':{'arc':'gsS20070623S0091','sky':'170:380,920:1080'}
    }
    for targ,p in sciTargets.iteritems():
        qs = qd['Full']
        qs['Object'] = targ
        # Fix up the target name for the output file
        sciOut = targ.split('-')[0]+targ[-1]
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qs), qs)
        gemtools.gemcombine(','.join(prefix+str(x) for x in sciFiles), 
                             sciOut, **sciCombFlags)
        gmos.gstransform(sciOut, wavtraname=p['arc'], **transFlags)
        gmos.gsskysub('t'+sciOut, long_sample=p['sky'], **skyFlags)

    # Clean up
    iraf.imdel("gsS2007*.fits")

    ## Apply the sensitivity function.
    gmos.gscalibrate.unlearn()
    calibFlags = {
        'extinction':'onedstds$ctioextinct.dat','fl_ext':'yes','fl_scale':'no', 
        'sfunction':'sens','fl_vardq':'yes','logfile':'gscalibrateLog.txt'
        }
    gmos.gscalibrate('stAM2306*', **calibFlags)
    calibFlags.update({'fl_vardq':'no'})
    gmos.gscalibrate('estLTT9239', **calibFlags)

    print (" -- Extract Target Spectra --")
    onedspec.nsum=4
    onedspec.sarith('cstAM2306b.fits[SCI]', 'copy', '', 'ecstAM2306b.ms',
                          apertures='222-346x4')
    stop
    print ("=== Finished Calibration Processing ===")


def gmos_ls_proc2(
        sciTargets,
        stdTarget,
        dbFile='./raw/obsLog.sqlite3',
        qd_full={'use_me': 1, 'Instrument': 'GMOS-S', 'CcdBin': '2 4', 'RoI': 'Full', 'Disperser': 'B600+_%', 'CentWave': 485.0, 'AperMask': '1.0arcsec', 'Object': 'AM2306-72%','DateObs': '2007-06-05:2007-07-07'},
        qd_censp={'use_me': 1, 'Instrument': 'GMOS-S', 'CcdBin': '2 4', 'RoI': 'CenSp', 'Disperser': 'B600+_%', 'CentWave': 485.0, 'AperMask': '1.0arcsec', 'Object': 'LTT9239','DateObs': '2007-06-05:2007-07-07'},
        biasFlags={'logfile': 'biasLog.txt', 'rawpath': './raw/', 'fl_vardq': 'yes', 'verbose': 'no'},
        flatFlags = {'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias': 'yes', 'fl_dark': 'no', 'fl_fixpix': 'no', 'fl_oversize': 'no', 'fl_vardq': 'yes', 'fl_fulldq': 'yes','rawpath': './raw', 'fl_inter': 'no', 'fl_detec': 'yes', 'function': 'spline3', 'order': '13,11,28', 'logfile': 'gsflatLog.txt', 'verbose': 'no'},
        sciFlags = {'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias': 'yes', 'fl_gscrrej': 'no','fl_dark': 'no', 'fl_flat': 'yes', 'fl_gmosaic': 'yes', 'fl_fixpix': 'no', 'fl_gsappwave': 'yes', 'fl_oversize': 'no', 'fl_vardq': 'yes', 'fl_fulldq': 'yes', 'rawpath': './raw', 'fl_inter': 'no', 'logfile': 'gsreduceLog.txt', 'verbose': 'no'},
        waveFlags = {'coordlist': 'gmos$data/CuAr_GMOS.dat', 'fwidth': 6, 'nsum': 50, 'function': 'chebyshev', 'order': 5, 'fl_inter': 'no', 'logfile': 'gswaveLog.txt', 'verbose': 'no'},
        sciCombFlags = {'combine': 'average', 'reject': 'ccdclip', 'fl_vardq': 'yes', 'fl_dqprop': 'yes', 'logfile': 'gemcombineLog.txt', 'verbose': 'no'},
        transFlags={'fl_vardq': 'yes', 'interptype': 'linear', 'fl_flux': 'yes', 'logfile': 'gstransLog.txt'},
        skyFlags={'fl_oversize': 'no', 'fl_vardq': 'yes', 'logfile': 'gsskysubLog.txt'},
        skip_wavecal=True,
        clean_files=False):

    """
    Parameters
    ----------
    dbFile : str
        Filename containing the SQL sqlite3 database created by obslog.py
        It must be placed in the ./raw/ directory
        Default is `./raw/obsLog.sqlite3`

    sciTargets : dict
        Dictionary with the associations of science targets and its associated ARC for wavelength calibration
        as well as the regions defining the sky along the slit.
        e.g. sciTargetd = {'AM2306-721_a': {'arc': 'gsS20070623S0071', 'sky': '520:720'},
                           'AM2306-72_b': {'arc': 'gsS20070623S0081', 'sky': '670:760,920:1020'}}
        Note that there could be more than one target defined this way.

    stdTarget : dict
        Dictionary with the associations of standard star targets and its associated ARC for wavelength calibration
        as well as the regions defining the sky along the slit.
        e.g. stdTarget = {'LTT1788': {'arc': 'S20180711S0281', 'sky': '170:380,920:1080'}}

    qd_full : dictionary
        Query Dictionary of essential parameter=value pairs for Full RoI. Meant for science object.

    qd_censp : dictionary
        Query Dictionary of essential parameter=value pairs for CenSp RoI. Meant for standard star.

    biasFlags : dict
        Dictionary for the keyword flags of gmos.gbias() function

    flatFlags : dict
        Dictionary for the keyword flags of gmos.gsflat() function

    sciFlags : dict
        Dictionary for the keyword flags of gmos.gsreduce() function
        Based on these flags a set of arcFlags and stdFlags dictionaries will be created
        for basic processing.

    waveFlags : dict
        Dictionary for the keyword flags of gmos.gswavelength() function

    sciCombFlags : dict
        Dictionary for the keyword flags of gemtools.gemcombine() function
        Based on these flags a set of stdCombFlags dictionary will be created for the standard advanced processing.

    transFlags : dict
        Dictionary for the keyword flags of gmos.gstransform() function.
        xxx

    skyFlags : dict
        Dictionary for the keyword flags of gmos.gsskysub() function

    skip_wavecal : bool
        Whether to skip interactive wavelength calibration.
        Useful when this is already done.


    Returns
    -------

    """

    print ("### Begin Processing GMOS/Longslit Images ###")
    print ("###")
    print ("=== Creating MasterCals ===")

    # From the work_directory:
    # Create the query dictionary of essential parameter=value pairs for Full and CenSp RoIs
    qd = {'Full': qd_full, 'CenSp': qd_censp}

    print (" --Creating Bias MasterCal--")

    # Set the task parameters.
    gemtools.gemextn.unlearn()  # Disarm a bug in gbias
    gmos.gbias.unlearn()

    regions = ['Full', 'CenSp']
    for r in regions:
        # The following SQL generates the list of full-frame files to process.
        SQL = fs.createQuery('bias', qd[r])
        biasFiles = fs.fileListQuery(dbFile, SQL, qd[r])

        # The str.join() funciton is needed to transform a python list into a
        # comma-separated string of file names that IRAF can understand.
        if len(biasFiles) > 1:
            # NT comment: sometimes if there are too many files, gmos.gbias() raises an error.
            # import pdb; pdb.set_trace()
            gmos.gbias(','.join(str(x) for x in biasFiles), 'MCbias' + r,
                       **biasFlags)

    # Clean up
    year_obs = qd_full['DateObs'].split('-')[0]
    if clean_files:
        iraf.imdel("gS{}*.fits".format(year_obs))

    ask_user("MC Bias done. Would you like to continue to proceed with GCAL Spectral Master Flats? (y/n): ",['y','yes'])

    print (" -- Creating GCAL Spectral Flat-Field MasterCals --")
    # Set the task parameters.
    qd['Full'].update({'DateObs': '*'})
    qd['CenSp'].update({'DateObs': '*'})
    gmos.gireduce.unlearn()
    gmos.gsflat.unlearn()
    # Normalize the spectral flats per CCD.
    # The response fitting should be done interactively.
    if flatFlags['fl_inter'] != 'yes':
        print("The response fitting should be done interactively. Please set flatFlags['fl_inter'] = 'yes'.")
        ask_user("Do you still want to proceed despite this important warning? (y/n): ", ['yes','y'])

    for r in regions:
        qr = qd[r]
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qr), qr)
        if len(flatFiles) > 0:
            gmos.gsflat(','.join(str(x) for x in flatFiles), 'MCflat' + r,
                        bias='MCbias' + r, **flatFlags)

    if clean_files:
        iraf.imdel('gS{}*.fits,gsS{}*.fits'.format(year_obs, year_obs))

    ask_user("GCAL Spectral Flat-Field MasterCals done. Would you like to continue to proceed with Basic Processing? (y/n): ",['y','yes'])

    print ("=== Processing Science Files ===")
    print (" -- Performing Basic Processing --")
    # Set task parameters.
    gmos.gsreduce.unlearn()
    sciFlags = sciFlags  # redundant but put here because NT likes it
    arcFlags = copy.deepcopy(sciFlags)
    arcFlags.update({'fl_flat': 'no', 'fl_vardq': 'no', 'fl_fulldq': 'no'})
    stdFlags = copy.deepcopy(sciFlags)
    stdFlags.update({'fl_fixpix': 'yes', 'fl_vardq': 'no', 'fl_fulldq': 'no'})

    # Perform basic reductions on all exposures for science targets.
    print ("  - Arc exposures -")
    for r in regions:
        qr = qd[r]
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qr), qr)
        if len(arcFiles) > 0:
            gmos.gsreduce(','.join(str(x) for x in arcFiles),
                          bias='MCbias' + r, **arcFlags)

    print ("  - Std star exposures -")
    r = 'CenSp'
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qd[r]), qd[r])
    if len(stdFiles) > 0:
        gmos.gsreduce(','.join(str(x) for x in stdFiles), bias='MCbias' + r,
                      flatim='MCflat' + r, **stdFlags)

    print ("  - Science exposures -")
    r = 'Full'
    sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qd[r]), qd[r])
    if len(sciFiles) > 0:
        gmos.gsreduce(','.join(str(x) for x in sciFiles), bias='MCbias' + r,
                      flatim='MCflat' + r, **sciFlags)

    # Clean up
    if clean_files:
        iraf.imdel('gS{}*.fits'.format(year_obs))

    ask_user("Basic processing done. Would you like to continue to determine wavelength calibration? (y/n): ",['y','yes'])

    print (" -- Determine wavelength calibration --")
    # Set task parameters
    gmos.gswavelength.unlearn()

    # The fit to the dispersion relation should be performed interactively.
    # Here we will use a previously determined result.
    if waveFlags['fl_inter'] != 'yes':
        print("The fit to the dispersion relation should be performed interactively. Please set waveFlags['fl_inter'] = 'yes'.")
        ask_user("Do you still want to proceed despite this important warning? (y/n): ", ['yes','y'])

    # Need to select specific wavecals to match science exposures.
    # NT: we do this now from the sciTargets + stdTarget input dictionaries
    # e.g.
        '''
        sciTargets = {
        'AM2306-721_a': {'arc': 'gsS20070623S0071', 'sky': '520:720'},
        'AM2306-72_b': {'arc': 'gsS20070623S0081', 'sky': '670:760,920:1020'},
        'AM2306-721_c': {'arc': 'gsS20070623S0091', 'sky': '170:380,920:1080'}
        }
        '''
    #prefix = 'gsS20070623S0'
    #for arc in ['071', '081', '091', '109']:
    #    gmos.gswavelength(prefix + arc, **waveFlags)
    prefix = 'gs'
    arc_files = []
    for key in sciTargets.keys():
        arc_files += [sciTargets[key]['arc']]
    for key in stdTarget.keys():
        arc_files += [stdTarget[key]['arc']]
    # import pdb; pdb.set_trace()
    if skip_wavecal is not True:
        for arc in arc_files:
            gmos.gswavelength(prefix + arc, **waveFlags)


    ### End of basic processing. Continue with advanced processing.
    ask_user("Wavelength solution done. Would you like to continue with advanced processing? (y/n): ",['y','yes'])

    print (" -- Performing Advanced Processing --")
    print (" -- Combine exposures, apply dispersion, subtract sky --")
    # Set task parameters.
    gemtools.gemcombine.unlearn()
    sciCombFlags = sciCombFlags
    stdCombFlags = copy.deepcopy(sciCombFlags)
    stdCombFlags.update({'fl_vardq': 'no', 'fl_dqprop': 'no'})
    gmos.gstransform.unlearn()

    # apply gtransform to standard
    # Process the Standard Star
    prefix = "gs"
    qs = qd['CenSp']
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qs), qs)
    std_name = stdTarget.keys()[0]
    if len(stdFiles) == 0:
        ValueError("No standard star associated. Please check parameters of search (e.g. RoI=CentSp)")
    import pdb; pdb.set_trace()
    if len(stdFiles) > 1:
        # import pdb; pdb.set_trace()
        gemtools.gemcombine(','.join(prefix + str(x) for x in stdFiles),
                                std_name, **stdCombFlags)
    else:
        os.system("cp {}.fits {}.fits".format(prefix + stdFiles[0], std_name))

    gmos.gstransform(std_name, wavtraname=prefix + stdTarget[std_name]['arc'], **transFlags)

    # The sky regions should be selected with care, using e.g. prows/pcols:
    #   pcols ("tAM2306b.fits[SCI]", 1100, 2040, wy1=40, wy2=320)
    print("The sky regions should be selected with care, using e.g. with prows/pcols (see tutorial).")
    answer = raw_input("Please provide the long_sample string to apply to gmos.gsskysub() for the standard star."
                       "e.g. '20:70,190:230'. Say 'no' for using the example as the default values.")
    if answer in ['n', 'no']:
        print("Using default long_sample set by stdTarget values {}.".format(stdTarget[std_name]['sky']))
        long_sample_std = stdTarget[std_name]['sky']
    else:
        long_sample_std = answer

    ask_user("Before proceeding it is important that you have set a good sky region for the standard.\n"
             "Thus far you have selected: {}\n Would you like to proceed with the current one? (y/n): ".format(long_sample_std), ['yes','y'])


    # apply sky substraction
    skyFlags = skyFlags
    gmos.gsskysub.unlearn()
    gmos.gsskysub('t{}'.format(std_name), long_sample=long_sample_std)

    # NT: make sure the process works ok until here before proceeding further. i.e. setting the sky region manually and correctly.
    # NT: seems to be working.

    print (" -- Extract Std spectrum --")
    # Extract the std spectrum using a large aperture.
    # It's important to trace the spectra interactively.
    gmos.gsextract.unlearn()
    extrFlags = {
        'apwidth': 3., 'fl_inter': 'yes', 'find': 'yes',
        'trace': 'yes', 'tfunction': 'chebyshev', 'torder': '6', 'tnsum': 20,
        'background': 'fit', 'bfunction': 'chebyshev', 'border': 2,
        'fl_vardq': 'no', 'logfile': 'gsextrLog.txt'
    }
    gmos.gsextract("st" + std_name, **extrFlags)

    stop
    print (" -- Derive the Flux calibration --")
    gmos.gsstandard.unlearn()
    sensFlags = {
        'fl_inter': 'yes', 'starname': 'l9239', 'caldir': 'onedstds$ctionewcal/',
        'observatory': 'Gemini-South', 'extinction': 'onedstds$ctioextinct.dat',
        'function': 'chebyshev', 'order': 9, 'verbose': 'no', 'logfile': 'gsstdLog.txt'
    }
    gmos.gsstandard('estLTT9239', sfile='std.txt', sfunction='sens', **sensFlags)

    # Process the science targets.
    # Use a dictionary to associate science targets with Arcs and sky regions.
    sciTargets = {
        'AM2306-721_a': {'arc': 'gsS20070623S0071', 'sky': '520:720'},
        'AM2306-72_b': {'arc': 'gsS20070623S0081', 'sky': '670:760,920:1020'},
        'AM2306-721_c': {'arc': 'gsS20070623S0091', 'sky': '170:380,920:1080'}
    }
    for targ, p in sciTargets.iteritems():
        qs = qd['Full']
        qs['Object'] = targ
        # Fix up the target name for the output file
        sciOut = targ.split('-')[0] + targ[-1]
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qs), qs)
        gemtools.gemcombine(','.join(prefix + str(x) for x in sciFiles),
                            sciOut, **sciCombFlags)
        gmos.gstransform(sciOut, wavtraname=p['arc'], **transFlags)
        gmos.gsskysub('t' + sciOut, long_sample=p['sky'], **skyFlags)

    # Clean up
    if clean_files:
        iraf.imdel("gsS{}*.fits".format(year_obs))

    ## Apply the sensitivity function.
    gmos.gscalibrate.unlearn()
    calibFlags = {
        'extinction': 'onedstds$ctioextinct.dat', 'fl_ext': 'yes', 'fl_scale': 'no',
        'sfunction': 'sens', 'fl_vardq': 'yes', 'logfile': 'gscalibrateLog.txt'
    }
    gmos.gscalibrate('stAM2306*', **calibFlags)
    calibFlags.update({'fl_vardq': 'no'})
    gmos.gscalibrate('estLTT9239', **calibFlags)

    print (" -- Extract Target Spectra --")
    onedspec.nsum = 4
    onedspec.sarith('cstAM2306b.fits[SCI]', 'copy', '', 'ecstAM2306b.ms',
                    apertures='222-346x4')

    print ("=== Finished Calibration Processing ===")

if __name__ == "__main__":
    gmos_ls_proc2()
