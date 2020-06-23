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

#nicersoft psrpipe wrapper
def run_psrpipe(obsID, par, 
				emin=0.25, emax=12,
				min_sun=60,
				keith=True):

	log.info("Running pspipe with par " + par)

	if not os.path.isdir(obsID):
		print("obsID not found")
		return -1
	
	cmd = ['psrpipe.py', 
		   '--emin', str(emin), 
		   '--emax', str(emax),
		   '--minsun', str(min_sun)]

    #Add par file if it is given
	if par is not None:
		cmd.extend(['--par', par])
	else:
		log.info("No par file input: pulse phase will not be defined")

	if keith:
		cmd.append('--keith')

	cmd.append(obsID)

	subprocess.call(cmd)

#runs full pipe on single observation
def allprocedures(obsID, par, 
				  emin=0.25, emax=12,
				  trumpet=True, keith=True, 
				  clobber=False):
	
    #Various clobber and path checks
    if ((os.path.isdir("{}_pipe".format(obsID))) 
        and (not clobber)):
        return 0
    elif ((os.path.isdir("{}_pipe".format(obsID)))
            and clobber):
        log.info("Removing {}_pipe".format(obsID))
        shutil.rmtree("{}_pipe".format(obsID))

    if not os.path.isdir(obsID):
        log.error("obsID not found")
        return -1

    #Some observations are missing event cl directory
    #Exit here if that's the case
    if not os.path.isdir(obsID+'/xti/event_cl'):
        log.error("no event cl dir for {}".format(obsID))
        return -1

    #Run nicerl2
    if (not pipeline_utils.check_nicerl2(obsID)) or clobber:
        pipeline_utils.run_nicerl2(obsID, trumpet=trumpet, clobber=clobber)

    #Add KP info
    pipeline_utils.run_add_kp(obsID)

    #Run psrpipe
    run_psrpipe(obsID, par,
                emin=emin, emax=emax,
                keith=keith)
    return 1

#Multiprocessing wrapper
def wrapper(par, emin=0.25, emax=12,
			trumpet=True, keith=True, clobber=False,
            crab=False):

    #Create pool
    pool = mp.Pool(processes=mp.cpu_count()+2)
    jobs = []

    #We have to find par files for crab
    if crab:
        log.info("Processing Crab data: parsing par data")
        df = pd.read_csv(par, keep_default_na=False)
        par_dic = {}
        for i in range(len(df)):
            if df['par'][i] == '':
                continue
            else:
                par_dic[str(df['obsID'][i])] = df['par'][i]

    #Find observations to pipe
    for f in os.listdir("./"):
        if (os.path.isdir(f)) and (f.isnumeric()):

            #par check for crab
            if crab:
                if f not in list(par_dic.keys()):
                    log.info("No par for {}".format(f))
                    continue
                else:
                    par = par_dic[f]

            #Add job to pool
            job = pool.apply_async(allprocedures,
                                   (f, par, ),
                                   dict(emin=emin,
                                        emax=emax,
                                        trumpet=trumpet,
                                        keith=keith,
                                        clobber=clobber))
            jobs.append(job)

    #Run all jobs in parallel
    for job in jobs:
        job.get()

#niextract merge event file
#count rate cut applied in this step
def run_niextract(sourcename, max_date=None, output=None, 
                  cut=2.0, filterbinsize=16):

    log.info('Running niextract-events')

    source_dir = os.getcwd().split('/')[-1]
    if sourcename != source_dir:
        check = pipeline_utils.outdircheck(sourcename)
        if check == -1:
            return

    #Build list of files to merge
    filenames = []
    invalid_columns = []
    for f in os.listdir("./"):
        if os.path.isfile(f+"/cleanfilt.evt"):
            tab = Table.read(f+"/cleanfilt.evt", hdu=1)
            date = datetime.datetime.strptime(
                    tab.meta['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
            #Date check for merge
            if max_date is not None:
                if date > max_date:
                    print("Excluding "+f+" from "+str(date))
                    continue
            #Some evts are missing columns...remember to look into this later
            if len(tab.columns) == 15:
                filenames.append(f+"/cleanfilt.evt")
            else:
                invalid_columns.append(f)

    #Write event file list
    with open(sourcename+"_evtlist", 'w') as f1:
        for fn in filenames:
            f1.write(fn+"\n")

    #Default output
    if output is None:
        output = sourcename+"_combined.evt"

    #Call niextract-events to merge files
    cmd = ['niextract-events', 'filename=@{}_evtlist'.format(sourcename),
           'eventsout={}'.format(output), 'clobber=yes']
    subprocess.call(cmd)

    #Save invalid column info
    with open('invalid_columns.txt', 'w') as f2:
        for fn in invalid_columns:
            f2.write(fn+"\n")

    #Run count rate cut
    if cut is not None:
        run_cr_cut(output, cut, filterbinsize)

#Perform count rate cut
def run_cr_cut(merged_evt, cut, filterbinsize):
	log.info("Running cr_cut cut={0} filterbinsize={1}".format(
                    str(cut), str(filterbinsize)))

	cmd = ['cr_cut.py', 
		   '--cut', str(cut), 
		   '--filterbinsize', str(filterbinsize),
		   merged_evt]
	
	subprocess.call(cmd)

#Update directory to download/pipe new data
def update(sourcename, heasarc_user, heasarc_pwd, decryptkey,
		   par, emin=0.25, emax=12, cut=2, filterbinsize=16,
		   trumpet=True, keith=True, crab=False, silent_curl=False):

    #Check the directory is correct
    source_dir = os.getcwd().split('/')[-1]
    log.info("Checking sourcename")
    if sourcename != source_dir:
        check = pipeline_utils.outdircheck(sourcename)
        if not check:
            return

    #Download new data
    pipeline_utils.run_datadownload(sourcename, heasarc_user, 
            heasarc_pwd, './', 
            decryptkey, clobber=False, 
            silent_curl=silent_curl)

    #Launch multiprocessing wrapper for pip
    wrapper(par, emin=emin, emax=emax, trumpet=trumpet, keith=keith,
            clobber=False, crab=crab)

    #Merge and perform count rate cut
    run_niextract(sourcename, cut=cut, filterbinsize=filterbinsize)

    #Create merged evt backup
    merged_evt = sourcename+"_combined_cut.evt"
    message='Created with pulsar_pipe.update'
    pipeline_utils.product_backup(merged_evt, backupdir='evt_backups',
                                  message=message)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)

    #Arguments for accessing data from NASA site
    parser.add_argument("--user", help="heasarc username", 
                        default=None, type=str)
    parser.add_argument("--passwd", help="heasarc password", 
                        default=None, type=str)
    parser.add_argument("--outdir", help="output dir for download", 
                        default="./", type=str)
    parser.add_argument("-k", help="decryption key", dest='key',
                        action=pipeline_utils.PasswordPromptAction, 
                        type=str, required=False)
    parser.add_argument("--k_file", help="decryption key from file",
                        type=str, default=None)

    #Arguments for observation filtering
    parser.add_argument("--emin", help="Minimum energy to include",
                        type=float, default=0.25)
    parser.add_argument("--emax", help="Maximum energy to include", 
                        type=float, default=12.0)
    parser.add_argument("--par", help="Par file to use for phases", 
                        type=str, default="")
    parser.add_argument("--cut", help="Count rate for cut in cts/sec", 
                        default=2.0, type=float)
    parser.add_argument("--filterbinsize", help="Bin size in seconds", 
                        default=16.0, type=float)
    parser.add_argument("--trumpet", help="Run nicerl2 with trumpet cut",
                        default=True, type=bool)
    parser.add_argument("--keith", 
                        help="Use standard filters from Keith Gendreau",
                        default=True, type=bool)
    parser.add_argument("--clobber", help="overwrite download files", 
                        default=False, action='store_true')

    #Arguments for different function calls
    parser.add_argument("--update", default=None, type=str,
                        help="Update single source with full procedure")
    parser.add_argument("--obsID", help="ObsID for single psrpipe call", 
                        default=None, type=str)
    parser.add_argument("--download", help="Run data download for source", 
                        default=None, type=str)
    parser.add_argument("--backup", help="Create EVT backup for file",
                        default=None, type=str)
    parser.add_argument("--message", help="Log message for evt backup",
                        default="", type=str)
    #Two additional arguments for merge call
    parser.add_argument("--merge", help="combine evt with niextract-events", 
                        default=None, type=str)
    parser.add_argument("--max_date", help="Max date to merge",
                        default=None, type=str)
    parser.add_argument("--output", help="Output name for merge",
                        default=None, type=str)

    #Option to silence curl progress bar
    parser.add_argument("--silent_curl", help="no curl progress bar",
                        default=False, action='store_true')

    args = parser.parse_args()

    #Key handling
    if args.k_file is not None:
        with open(args.k_file, 'r') as f:
            args.key = f.readlines()[0].strip("\n")

    #Update directory and pipe all new data
    if args.update is not None:
        update(args.update, args.user, args.passwd, 
               args.key, args.par, 
               emin=args.emin, emax=args.emax, 
               cut=args.cut, 
               filterbinsize=args.filterbinsize, 
               trumpet=args.trumpet,
               keith=args.keith,
               silent_curl=args.silent_curl)

    #Option to call pipeline on a single observation
    elif args.obsID is not None:
        allprocedures(args.obsID, args.par,
                      args.emin, args.emax, 
                      trumpet=args.trumpet, 
                      keith=args.keith,
                      clobber=args.clobber)
        if os.path.isfile(os.path.join(args.obsID+'_pipe', 'cleanfilt.evt')):
            check = input("Would you like to run count rate cut [y/n] -- ")
            if check in ['y', 'Y']:
                run_cr_cut(os.path.join(args.obsID+'_pipe', 'cleanfilt.evt'),
                                        args.cut, args.filterbinsize)

    #Option to only download the data with no filtering
    elif args.download is not None:
        assert(all([arg is not None for arg in [
                    args.user, args.passwd, args.key]]))
        pipeline_utils.run_datadownload(args.download,
                                        args.user, args.passwd, 
                                        args.outdir, 
                                        args.key, clobber=args.clobber,
                                        silent_curl=args.silent_curl)
    #Option to just merge
    #Optional arguments to merge up to specific date and give output name
    elif args.merge is not None:
        if args.max_date is None:
            run_niextract(args.merge, 
                          cut=args.cut, filterbinsize=args.filterbinsize,
                          output=args.output)
        else:
            dt = datetime.datetime.strptime(args.max_date, '%Y-%m-%d')
            run_niextract(args.merge, max_date=dt, 
                          cut=args.cut, filterbinsize=args.filterbinsize,
                          output=args.output)

    #Back up event file
    elif args.backup is not None:
        pipeline_utils.product_backup(args.backup, backupdir='evt_backups',
                                      message=args.message)

    else:
        log.warning("No method selected! Try 'pulsar_pipe.py --help'")

