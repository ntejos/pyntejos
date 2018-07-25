#! /usr/bin/env python
#    2016-Jun-08  shaw@noao.edu
# Modified by ntejos to make it a script

import sys
from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos
import fileSelect as fs

bpm_gmos="/home/ntejos/python/pyntejos/pyntejos/gmos/bpm_gmos-s_EEV_v1_2x2_img_MEF.fits"

def ask_user(question, good_answers):
    answer = raw_input(question)
    if answer in good_answers
        pass
    else:
        raise ValueError("Stopping process because user said so.")

# def gmos_img_proc():
if 0:
    '''
    GMOS Data Reduction Cookbook companion script to the chapter:
      "Reduction of Images with PyRAF"

    PyRAF script to:
    Process GMOS images for Messier 8, in program GS-2006B-Q-18.

    The names for the relevant header keywords and their expected values are
    described in the Cookbook chapter entitled "Supplementary Material"

    Perform the following starting in the parent work directory:
      cd /path/to/work_directory

    Fetch the Static BPM file from the tutorial, place it in your work directory,
    and uncompress:
       bpm_gmos-s_EEV_v1_2x2_img_MEF.fits
       
    Place the fileSelect.py module in your work directory. 
    You may cut-n-paste lines from this script to your pyraf session, or 
    from the unix prompt: 
       python gmos_img_proc.py
    '''

    print ("### Begin Processing GMOS/MOS Images ###")
    print ("###")
    print ("=== Creating MasterCals ===")

    # This whole example depends upon first having built an sqlite3 database of metadata:
    #    cd ./raw
    #    python obslog.py obsLog.sqlite3
    dbFile='./raw/obsLog.sqlite3'

    # From the work_directory:
    # Create the query dictionary of essential parameter=value pairs.
    # Select bias exposures within ~2 months of the target observations:
    qd = {'use_me':1,
          'Instrument':'GMOS-S','CcdBin':'2 2','RoI':'Full','Object':'M8-%',
          'DateObs':'2006-09-01:2006-10-30'
          }
    print (" --Creating Bias MasterCal--")

    # Set the task parameters.
    gmos.gbias.unlearn()
    biasFlags = {
        'logfile':'biasLog.txt','rawpath':'./raw/','fl_vardq':'yes',
        'verbose':'no'
    }
    # The following SQL generates the list of files to process.
    SQL = fs.createQuery('bias', qd)
    biasFiles = fs.fileListQuery(dbFile, SQL, qd)
    
    # The str.join() function is needed to transform a python list into a string
    # filelist that IRAF can understand.
    if len(biasFiles) > 1:
        gmos.gbias(','.join(str(x) for x in biasFiles), 'MCbias.fits',
            **biasFlags)

    # Clean up
    iraf.imdel('gS2006*.fits')

    print (" --Creating Twilight Imaging Flat-Field MasterCal--")
    # Select flats obtained contemporaneously with the observations.
    qd.update({'DateObs':'2006-09-10:2006-10-10'})

    # Set the task parameters.
    gmos.giflat.unlearn()
    flatFlags = {
        'fl_scale':'yes','sctype':'mean','fl_vardq':'yes',
        'rawpath':'./raw/','logfile':'giflatLog.txt','verbose':'yes'
        }
    filters = ['Ha', 'HaC', 'SII', 'r', 'i']
    for f in filters:
        print "  Building twilight flat MasterCal for: %s" % (f)

        # Select filter name using a substring of the official designation.
        qd['Filter2'] = f + '_G%'
        mcName = 'MCflat_%s.fits' % (f)
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('twiFlat', qd), qd)
        if len(flatFiles) > 0:
            gmos.giflat(','.join(str(x) for x in flatFiles), mcName, 
                         bias='MCbias', **flatFlags)

    iraf.imdel('gS2006*.fits,rgS2006*.fits')
    stop
    print ("=== Processing Science Images ===")
    # Remove restriction on date range
    qd['DateObs'] = '*'
    prefix = 'rg'

    # Set task parameters.
    # Employ the imaging Static BPM for this set of detectors.
    gmos.gireduce.unlearn()
    sciFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_dark':'no',
        'fl_flat':'yes','logfile':'gireduceLog.txt','rawpath':'./raw/',
        'fl_vardq':'yes','bpm':'bpm_gmos-s_EEV_v1_2x2_img_MEF.fits','verbose':'no'
        }
    gemtools.gemextn.unlearn()    # disarms a bug in gmosaic
    gmos.gmosaic.unlearn()
    mosaicFlags = {
        'fl_paste':'no','fl_fixpix':'no','fl_clean':'yes','geointer':'nearest',
        'logfile':'gmosaicLog.txt','fl_vardq':'yes','fl_fulldq':'yes','verbose':'no'
        }
    # Reduce the science images, then mosaic the extensions in a loop
    for f in filters:
        print "    Processing science images for: %s" % (f)
        qd['Filter2'] = f + '_G%'
        flatFile = 'MCflat_' + f + '.fits'
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciImg', qd), qd)
        if len(sciFiles) > 0:
            gmos.gireduce (','.join(str(x) for x in sciFiles), bias='MCbias', 
                           flat1=flatFile, **sciFlags)
            for file in sciFiles:
                gmos.gmosaic (prefix+file, **mosaicFlags)

    iraf.imdelete('gS2006*.fits,rgS2006*.fits')

    ## Co-add the images, per position and filter.
    print (" -- Begin image co-addition --")

    # Use primarily the default task parameters.
    gemtools.imcoadd.unlearn()
    coaddFlags = {
        'fwhm':3,'datamax':6.e4,'geointer':'nearest','logfile':'imcoaddLog.txt'
        }
    targets = ['M8-1', 'M8-2', 'M8-3']
    prefix = 'mrg'
    for f in filters:
        print "  - Co-addding science images in filter: %s" % (f)
        qd['Filter2'] = f + '_G%'
        for t in targets:
            qd['Object'] = t + '%'
            print "  - Co-addding science images for position: %s" % (t)
            outImage = t + '_' + f + '.fits'
            coAddFiles = fs.fileListQuery(dbFile, fs.createQuery('sciImg', qd), qd)
            gemtools.imcoadd(','.join(prefix+str(x) for x in coAddFiles),
                    outimage=outImage, **coaddFlags)

    iraf.delete ("*_trn*,*_pos,*_cen")
    iraf.imdelete ("*badpix.pl,*_med.fits,*_mag.fits")
    #iraf.imdelete ("mrgS*.fits")

    print ("=== Finished Calibration Processing ===")

def gmos_img_proc2(dbFile="./raw/obsLog.sqlite3", qd={'use_me': 1,'Instrument': 'GMOS-S', 'CcdBin': '2 2', 'RoI': 'Full', 'Object': 'M8-%', 'DateObs': '2006-09-01:2006-10-30'},
                   bias_dateobs="2006-09-01:2006-10-30",
                   biasFlags={'logfile': 'biasLog.txt', 'rawpath': './raw/', 'fl_vardq': 'yes', 'verbose': 'yes'},
                   flat_dateobs='2006-09-10:2006-10-10',
                   flatFlags = {'fl_scale': 'yes', 'sctype': 'mean', 'fl_vardq': 'yes','rawpath': './raw/', 'logfile': 'giflatLog.txt', 'verbose': 'yes'},
                   filters = ['Ha', 'HaC', 'SII', 'r', 'i'],
                   sciFlags={'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias':'yes', 'fl_dark': 'no','fl_flat': 'yes', 'logfile':'gireduceLog.txt', 'rawpath': './raw/','fl_vardq': 'yes','bpm':bpm_gmos, 'verbose': 'yes'},
                   mosaicFlags = {'fl_paste': 'no', 'fl_fixpix': 'no', 'fl_clean': 'yes', 'geointer': 'nearest', 'logfile': 'gmosaicLog.txt', 'fl_vardq': 'yes', 'fl_fulldq': 'yes', 'verbose': 'yes'},
                   coaddFlags = {'fwhm': 3, 'datamax': 6.e4, 'geointer': 'nearest', 'logfile': 'imcoaddLog.txt'},
                   targets = ['M8-1', 'M8-2', 'M8-3'],
                   clean_files = False
                   ):
    """
    Parameters
    ----------
    dbFile : str
        Filename containing the SQL sqlite3 database created by obslog.py
        It must be placed in the ./raw/ directory
        Default is `./raw/obsLog.sqlite3`
    qd : dictionary
        Query Dictionary of essential parameter=value pairs.
        Select bias exposures within ~2 months of the target observations
        e.g.
         qd= {'use_me': 1,
          'Instrument': 'GMOS-S', 'CcdBin': '2 2', 'RoI': 'Full', 'Object': 'M8-%',
          'DateObs': '2006-09-01:2006-10-30'
          }
    bias_dateobs : str
        String representing the bias search Obsdate
        e.g. bias_dateobs = `2006-09-01:2006-10-30`

    biasFlags : dict
        Dictionary for the keyword flags of gmos.gbias() function

    flat_dateobs : str
        String representing the flat search Obsdate
        e.g. flat_dateobs = `2006-09-10:2006-10-10`

    flatFlags : dict
        Dictionary for the keyword flags of gmos.giflat() function
        e.g. flatFlags = {'fl_scale': 'yes', 'sctype': 'mean', 'fl_vardq': 'yes','rawpath': './raw/',
                         'logfile': 'giflatLog.txt', 'verbose': 'yes'}
    filters : list
        List of filter names to perform reduction
        e.g. filters=['Ha', 'HaC', 'SII', 'r', 'i']

    sciFlags : dict
        Dictionary for the keyword flags of gmos.gireduce() function

    mosaicFlags : dict
        Dictionary for the keyword flags of gmos.gimosaic() function

    coaddFlags : dict
        Dictionary for the keyword flags of gemtools.imcoadd() function

    targets : list
        List of names of target observations for the co-addition
        e.g.  targets = ['M8-1', 'M8-2', 'M8-3']

    clean_files : bool
        Whether to clean intermediate files from reduction process

    Returns
    -------
    Reduce GMOS imaging based on tutorial example.


    """
    print ("### Begin Processing GMOS/MOS Images ###")
    print ("###")
    print ("=== Creating MasterCals ===")

    # From the work_directory:
    # Create the query dictionary of essential parameter=value pairs.
    # Select bias exposures within ~2 months of the target observations:

    print (" --Creating Bias MasterCal--")
    qd.update({'DateObs': bias_dateobs})

    # Set the task parameters.
    gmos.gbias.unlearn()
    # The following SQL generates the list of files to process.
    SQL = fs.createQuery('bias', qd)
    biasFiles = fs.fileListQuery(dbFile, SQL, qd)

    # The str.join() function is needed to transform a python list into a string
    # filelist that IRAF can understand.
    if len(biasFiles) > 1:
        files_all = ','.join(str(x) for x in biasFiles)
        # import pdb; pdb.set_trace()
        gmos.gbias(files_all, 'MCbias.fits', **biasFlags)

    # Clean up
    year_obs = qd['DateObs'].split('-')[0]
    if clean_files:
        iraf.imdel('gS{}*.fits'.format(year_obs))

    ask_user("MC Bias done. Would you like to continue to proceed with Master Flats? (y/n): ",['y','yes'])

    print (" --Creating Twilight Imaging Flat-Field MasterCal--")
    # Select flats obtained contemporaneously with the observations.
    qd.update({'DateObs': flat_dateobs})

    # Set the task parameters.
    gmos.giflat.unlearn()

    #filters = ['Ha', 'HaC', 'SII', 'r', 'i']
    for f in filters:
        print "  Building twilight flat MasterCal for filter: %s" % (f)

        # Select filter name using a substring of the official designation.
        qd['Filter2'] = f + '_G%'
        mcName = 'MCflat_%s.fits' % (f)
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('twiFlat', qd), qd)
        if len(flatFiles) > 0:
            files_all = ','.join(str(x) for x in flatFiles)
            # import pdb; pdb.set_trace()
            gmos.giflat(files_all, mcName, bias='MCbias', **flatFlags)

    if clean_files:
        #
        #iraf.imdel('gS{}*.fits,rgS{}*.fits'.format(year_obs, year_obs))
        pass

    ask_user("MC Flats done. Would you like to continue to proceed with processing Science Images? (y/n): ", ['yes','y'])

    print ("=== Processing Science Images ===")
    # Remove restriction on date range
    qd['DateObs'] = '*'
    prefix = 'rg'

    gmos.gireduce.unlearn()
    gemtools.gemextn.unlearn()  # disarms a bug in gmosaic
    gmos.gmosaic.unlearn()
    # Reduce the science images, then mosaic the extensions in a loop
    for f in filters:
        print "    Processing science images for filter: %s" % (f)
        qd['Filter2'] = f + '_G%'
        flatFile = 'MCflat_' + f + '.fits'
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciImg', qd), qd)
        if len(sciFiles) > 0:
            # Make sure BPM table is in sciFlags for employing the imaging Static BPM for this set of detectors.
            # import pdb; pdb.set_trace()
            all_files = ','.join(str(x) for x in sciFiles)
            gmos.gireduce(all_files, bias='MCbias', flat1=flatFile, **sciFlags)
            for file in sciFiles:
                gmos.gmosaic(prefix + file, **mosaicFlags)
        else:
            print("No Science images found for filter {}. Check database.".format(f))
            import pdb; pdb.set_trace()

    if clean_files:
        iraf.imdelete('gS{}*.fits,rgS{}*.fits'.format(year_obs,year_obs))

    ask_user("Science Images done. Would you like to continue to proceed with image co-addition? (y/n): ", ['y','yes'])

    ## Co-add the images, per position and filter.
    print (" -- Begin image co-addition --")

    # Use primarily the default task parameters.
    gemtools.imcoadd.unlearn()
    prefix = 'mrg'
    for f in filters:
        print "  - Co-addding science images in filter: %s" % (f)
        qd['Filter2'] = f + '_G%'
        for t in targets:
            qd['Object'] = t + '%'
            print "  - Co-addding science images for position: %s" % (t)
            outImage = t + '_' + f + '.fits'
            coAddFiles = fs.fileListQuery(dbFile, fs.createQuery('sciImg', qd), qd)
            all_files = ','.join(prefix + str(x) for x in coAddFiles)
            if all_files == '':
                print('No files available for co-addition...')
                import pdb; pdb.set_trace()
            gemtools.imcoadd(all_files, outimage=outImage, **coaddFlags)

    ask_user("Co-addition done. Would you like to clean intermediate reduction files? (y/n): ", ['y','yes'])

    if clean_files:
        iraf.delete("*_trn*,*_pos,*_cen")
        iraf.imdelete("*badpix.pl,*_med.fits,*_mag.fits")
        # iraf.imdelete ("mrgS*.fits")

    print ("=== Finished Calibration Processing ===")


if __name__ == "__main__":
    gmos_img_proc2()
