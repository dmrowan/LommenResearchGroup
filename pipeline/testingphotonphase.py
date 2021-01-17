#!/usr/bin/env python
from __future__ import print_function, division
from astropy.table import Table
from astropy import log
import datetime
import os
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import shutil
import subprocess
import time

#LommenResearchGroup imports
import pipeline_utils

#Dom Rowan 2020

desc = """
Pipeline for pulsar analysis
"""

def allprocedures(obsID, par, 
				  emin=0.25, emax=12,
				  trumpet=True, keith=True, 
				  clobber=False, crab=False):
	
    #Various clobber and path checks
    print('hello6')
    if ((os.path.isdir("{}_pipe".format(obsID))) 
        and (not clobber)):
        print('1')
        return 0
    elif ((os.path.isdir("{}_pipe".format(obsID)))
            and clobber):
        print('2')
        log.info("Removing {}_pipe".format(obsID))
        shutil.rmtree("{}_pipe".format(obsID))

    if not os.path.isdir(obsID):
        print('3')
        log.error("obsID not found")
        return -1

    #Some observations are missing event cl directory
    #Exit here if that's the case
    print('hello7')
    if not os.path.isdir(obsID+'/xti/event_cl'):
        log.error("no event cl dir for {}".format(obsID))
        return -1

    #Run nicerl2
    print('hello8')
    if (not pipeline_utils.check_nicerl2(obsID)) or clobber:
        pipeline_utils.run_nicerl2(obsID, trumpet=trumpet, clobber=clobber)

    #Add KP info
    print('hello9')
    pipeline_utils.run_add_kp(obsID)

    #Run psrpipe
    print('hello10')
    run_psrpipe(obsID, par,
                emin=emin, emax=emax,
                keith=keith, crab=crab)
    return 1

f= '1013010135'
par = '/students/pipeline/parfiles/crab/april2018.par'
allprocedures(f, par, clobber =  True, crab = False)
