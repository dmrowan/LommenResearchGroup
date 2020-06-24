#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import glob
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import subprocess

#LommenResearchGroup imports
from pipeline import pipeline_utils

#Dom Rowan 2020

desc="""
Functions for working with the crab data
"""

#Read par files to find start and finish dates
def crab_par_table(par_dir='/students/pipeline/parfiles/crab/'):
    
    #Crab all the par files 
    par_files = os.listdir(par_dir)
    par_files = [ os.path.join(par_dir, p) for p in par_files ]
    
    start = []
    finish = []
    #Read each
    for par in par_files:
        with open(par, 'r') as f:
            lines = f.readlines()

        #Find start and end, do date conversion
        for l in lines:
            if l.split()[0] == 'START':
                startval = float(l.split()[1].strip('\n'))
                start_dt = Time(startval, format='mjd').utc.datetime
                start.append(start_dt)
            elif l.split()[0] == 'FINISH':
                finishval = float(l.split()[1].strip('\n'))
                finish_dt = Time(finishval, format='mjd').utc.datetime
                finish.append(finish_dt)
            else:
                continue
		
    #Return in pandas df
    df = pd.DataFrame({'par':par_files, 'start':start, 'finish':finish})
    return df

#Match par files to obsIDs
def crab_par_match(par_dir='/students/pipeline/parfiles/crab', outdir='./', 
                  par_date_clearance=None, par_save='par_info'):
    
    #Get all par information
    par_df = pipeline_utils.crab_par_table(par_dir=par_dir)

    #Handling directory of observation
    if outdir == './':
        pre_split = os.getcwd()
    else:
        pre_split = outdir
    source_dir = pre_split.split('/')[-1]
    log.info("Checking sourcename and directory match")
    if 'PSR_B0531+21' != source_dir:
        check = pipeline_utils.outdircheck('PSR_B0531+21')
        if not check:
            return
    
    #Next get all crab obsIDs
    obsids = []
    for f in os.listdir(outdir):
        if (os.path.isdir(f)) and (not f.endswith('_pipe')) and (f!='tmp'):
            obsids.append(f)        

    #Need to read the mkfs to get date
    mkfs = [ os.path.join(outdir, ID, 'auxil/ni{0}.mkf'.format(ID))
             for ID in obsids ]

    obs_dates = []
    obs_par = []
    obs_clearance = [] #date leway clearance
    format_specifier = '%Y-%m-%dT%H:%M:%S'
    log.info("Matching observation dates with par files")
    #Go through each mkf
    for m in mkfs:
        if not os.path.isfile(m):
            log.error("No MKF found", m)
            obs_dates.append("-")
            obs_par.append("-")
            obs_clearance.append("-")
            continue
        tab = Table.read(m, hdu=1)
        date = tab.meta['DATE-OBS']
        date_formatted = datetime.datetime.strptime(date, format_specifier)
        obs_dates.append(date_formatted)
        par = None
        clearance = 0
        #Find the matching par file
        for i in range(len(par_df)):
            if par_df['start'][i] <= date_formatted <= par_df['finish'][i]:
                par = par_df['par'][i]

        #Use fuzzy bounds
        if (par is None) and (par_date_clearance is not None):
            diff = []
            for i in range(len(par_df)):
                if date_formatted < par_df['start'][i]:
                    diff.append((date_formatted - par_df['start'][i]).days)
                else:
                    assert(date_formatted > par_df['finish'][i])
                    diff.append((date_formatted - par_df['finish'][i]).days)

            diff_abs = list(map(abs, diff))
            idx_min = np.where(np.array(diff_abs) == min(diff_abs))[0][0]
            if abs(diff[idx_min]) <= par_date_clearance:
                par = par_df['par'][idx_min]
                clearance = diff[idx_min]

        obs_par.append(par)
        obs_clearance.append(clearance)
        
    #Save output in pandas df
    df_obs = pd.DataFrame({'obsID':obsids,
                           'date':obs_dates,
                           'par':obs_par,
                           'clearance':obs_clearance})

    #Write to file
    log.info("Writing par date match to file")
    if par_save is not None:
        df_obs.to_csv(par_save+'.csv', index=False)
        df_obs.to_latex(par_save+'.tex', index=False)

    return df_obs
