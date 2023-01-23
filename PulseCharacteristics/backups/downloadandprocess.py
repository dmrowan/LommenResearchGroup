#!/usr/bin/env python

import datetime
import numpy as np
import pandas as pd
import os
from pipeline.pulsarpipe import runcmd
from pipeline.pulsarpipe import processdata
from astropy import log
import pipeline.pipeline_utils2 as pipeline_utils2
import shutil
import fnmatch

desc="""
Downloads and processes Crab data files
"""

obsid = pd.read_csv('par_info.csv')
goodobsid = obsid.dropna()
obsids = list(goodobsid['obsID'])
parfiles = list(goodobsid['par'])

#Download data
for n in range(len(obsids)):
    try: 
        os.chdir('%s'%obsids[n])
        os.chdir('..')
        log.info("Data already downloaded for ObsID %s, skipping"%obsids[n])
        continue
    except FileNotFoundError:
        log.info("Downloading data for ObsID %s"%obsids[n])
    cmd = 'pulsar_pipe2.py --download PSR_B0531+21 --obsID %s --user nicer_team --passwd sextant --k_file keyfile.txt'%obsids[n]
    print('CMD: ', cmd)
    try:
        os.system(cmd)
    except NameError:
        log.info('Unzipping failed the first time, trying again')
    os.system(cmd)

obsids = []
dirs =[name for name in os.listdir(".") if os.path.isdir(name)]
for name in dirs:
    if (fnmatch.fnmatch(name, '*_pipe')) or (name == 'tmp'):
        pass
    else:
        obsids.append(name)

#Process data
for n in range(len(obsids)):
    try: 
        os.chdir('%s'%obsids[n])
        os.chdir('..')
    except FileNotFoundError:
        continue
    try:
        os.chdir('%s_pipe'%obsids[n])
        os.chdir('..')
        log.info('Data for ObsID %s already processed, skipping'%obsids[n])
        continue
    except FileNotFoundError:
        log.info("Processing data for ObsID %s"%obsids[n])
    try: 
        processdata(str(obsids[n]), parfiles[n])
    except ValueError:
        log.info("This ObsID has no good data, skipping")
        with open('badObsIDs.txt', 'a') as f:
            print(obsids[n], file=f)
        #shutil.rmtree('%s_pipe'%obsids[n])
        continue

log.info('DONE')
