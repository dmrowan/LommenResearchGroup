#!/usr/bin/env python

from astropy.table import Table
from astropy import log
from astropy.time import Time
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
            if l.isspace():
                continue
            elif l.split()[0] == 'START':
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
def crab_par_match(username, passwd, 
                   par_dir='/students/pipeline/parfiles/crab', 
                   par_date_clearance=None, par_save='par_info'):

    log.info("Matching observation dates with par files")

    #Get available par information
    par_df = crab_par_table(par_dir=par_dir)

    #Query nasa site
    segment_df = pipeline_utils.print_nicer_segment(username=username, password=passwd)

    #Select crab data
    df_crab = segment_df[segment_df['Target Name']=='PSR_B0531+21']

    #reset pandas index
    df_crab = df_crab.reset_index(drop=True)

    obs_dates = []
    obs_par = []
    obs_clearance = []
    format_specifier='%Y-%m-%dT%H:%M:%S'
    for i in range(len(df_crab)):
        date_formatted = datetime.datetime.strptime(df_crab['Start TimeUTC'][i], format_specifier)
        obs_dates.append(date_formatted)
        par = None
        clearance=0

        for j in range(len(par_df)):
            if par_df['start'][j] <= date_formatted <= par_df['finish'][j]:
                par = par_df['par'][j]

        #Use fuzzy bounds
        if (par is None) and (par_date_clearance is not None):
            diff = []
            for k in range(len(par_df)):
                if date_formatted < par_df['start'][k]:
                    diff.append((date_formatted - par_df['start'][k]).days)
                else:
                    assert(date_formatted > par_df['finish'][k])
                    diff.append((date_formatted - par_df['finish'][k]).days)

            diff_abs = list(map(abs, diff))
            idx_min = np.where(np.array(diff_abs) == min(diff_abs))[0][0]
            if abs(diff[idx_min]) <= par_date_clearance:
                par = par_df['par'][idx_min]
                clearance = diff[idx_min]
        
        obs_par.append(par)
        obs_clearance.append(clearance)

    #Save output to pandas df
    df_obs = pd.DataFrame({'obsID':df_crab['Observation ID'],
                           'date':obs_dates,
                           'par':obs_par,
                           'clearance':obs_clearance})

    #Write to file
    log.info("Writing par date match to file")
    if par_save is not None:
        df_obs.to_csv(par_save+'.csv', index=False)
        df_obs.to_latex(par_save+'.tex', index=False)

    return df_obs
