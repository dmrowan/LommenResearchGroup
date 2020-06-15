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
import subprocess
import getpass
import time

#LommenResearchGroup imports
import pipeline_utils

#Dom Rowan 2020

desc = """
Pipeline for pulsar analysis
"""

def crab_par_match(par_dir='/students/pipeline/parfiles/crab', outdir='./', 
                  par_date_clearance=None, par_save='par_info'):
    
    par_df = pipeline_utils.crab_par_table(par_dir=par_dir)

    if outdir == './':
        pre_split = os.getcwd()
    else:
        pre_split = outdir
    source_dir = pre_split.split('/')[-1]
    log.info("Checking sourcename and directory match")
    if 'PSR_B0531+21' != source_dir:
        check = outdircheck('PSR_B0531+21')
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
        for i in range(len(par_df)):
            if par_df['start'][i] <= date_formatted <= par_df['finish'][i]:
                par = par_df['par'][i]

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
        
    df_obs = pd.DataFrame({'obsID':obsids,
                           'date':obs_dates,
                           'par':obs_par,
                           'clearance':obs_clearance})

    log.info("Writing par date match to file")
    if par_save is not None:
        df_obs.to_csv(par_save+'.csv', index=False)
        df_obs.to_latex(par_save+'.tex', index=False)

    return df_obs


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

	if par is not None:
		cmd.extend(['--par', par])
	else:
		log.info("No par file input: pulse phase will not be defined")

	if keith:
		cmd.append('--keith')

	cmd.append(obsID)

	subprocess.call(cmd)

def allprocedures(obsID, par, 
				  emin=0.25, emax=12,
				  trumpet=True, keith=True, 
				  clobber=False):
	
    if ((os.path.isdir("{}_pipe".format(obsID))) 
        and (not clobber)):
        return 0
    elif ((os.path.isdir("{}_pipe".format(obsID)))
            and clobber):
        log.info("Removing {}_pipe".format(obsID))
        os.rmdir("{}_pipe".format(obsID))

    if not os.path.isdir(obsID):
        log.error("obsID not found")
        return -1

    if not os.path.isdir(obsID+'/xti/event_cl'):
        log.error("no event cl dir for {}".format(obsID))
        return -1

    if (not pipeline_utils.check_nicerl2(obsID)) or clobber:
        pipeline_utils.run_nicerl2(obsID, trumpet=trumpet, clobber=clobber)

    pipeline_utils.run_add_kp(obsID)

    run_psrpipe(obsID, par,
                emin=emin, emax=emax,
                keith=keith)
    return 1

def wrapper(par, emin=0.25, emax=12,
			trumpet=True, keith=True, clobber=False,
            crab=False):

    pool = mp.Pool(processes=mp.cpu_count()+2)
    jobs = []
    times = []

    if crab:
        log.info("Processing Crab data: parsing par data")
        df = pd.read_csv(par, keep_default_na=False)
        par_dic = {}
        for i in range(len(df)):
            if df['par'][i] == '':
                continue
            else:
                par_dic[str(df['obsID'][i])] = df['par'][i]

    for f in os.listdir("./"):
        if (os.path.isdir(f)) and (not f.endswith('_pipe')) and (f!='tmp'):

            if crab:
                if f not in list(par_dic.keys()):
                    log.info("No par for {}".format(f))
                    continue
                else:
                    par = par_dic[f]

            job = pool.apply_async(allprocedures,
                                   (f, par, ),
                                   dict(emin=emin,
                                        emax=emax,
                                        trumpet=trumpet,
                                        keith=keith,
                                        clobber=clobber))
            jobs.append(job)

    for job in jobs:
        job.get()

def run_niextract(sourcename, max_date=None):
	source_dir = os.getcwd().split('/')[-1]
	if sourcename != source_dir:
		check = outdircheck(sourcename)
		if check == -1:
			return
	filenames = []
	invalid_columns = []
	for f in os.listdir("./"):
		if os.path.isfile(f+"/cleanfilt.evt"):
			tab = Table.read(f+"/cleanfilt.evt", hdu=1)
			date = datetime.datetime.strptime(
                    tab.meta['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
			if max_date is not None:
				if date > max_date:
					print("Excluding "+f+" at date "+str(date))
					continue

			if len(tab.columns) == 15:
				filenames.append(f+"/cleanfilt.evt")
			else:
				invalid_columns.append(f)

	with open(sourcename+"_evtlist", 'w') as f1:
		for fn in filenames:
			f1.write(fn+"\n")
	output = sourcename+"_combined.evt"
	cmd = ['niextract-events', 'filename=@{}_evtlist'.format(sourcename),
		   'eventsout={}'.format(output), 'clobber=yes']
	subprocess.call(cmd)

	with open('invalid_columns.txt', 'w') as f2:
		for fn in invalid_columns:
			f2.write(fn+"\n")

def run_cr_cut(merged_evt, cut, filterbinsize):
	log.info("Running cr_cut")

	cmd = ['cr_cut.py', 
		   '--cut', str(cut), 
		   '--filterbinsize', str(filterbinsize),
		   merged_evt]
	
	subprocess.call(cmd)

def update(sourcename, heasarc_user, heasarc_pwd, decryptkey,
		   par, emin=0.25, emax=12, cut=2, filterbinsize=16,
		   trumpet=True, keith=True, crab=False):

    source_dir = os.getcwd().split('/')[-1]
    log.info("Checking sourcename")
    if sourcename != source_dir:
        check = outdircheck(sourcename)
        if not check:
            return

    pipeline_utils.run_datadownload(sourcename, heasarc_user, 
            heasarc_pwd, './', 
            decryptkey, clobber=False)

    wrapper(par, emin=emin, emax=emax, trumpet=trumpet, keith=keith,
            clobber=False, crab=crab)

    run_niextract(sourcename)

    run_cr_cut(sourcename+"_combined.evt", cut, filterbinsize)

    #Create merged evt backup
    merged_evt = sourcename+"_combined_cut.evt"
    message='Created with pulsar_pipe.update'
    pipeline_utils.product_backup(merged_evt, backupdir='evt_backups',
                                  message=message)

def cronjob(heasarc_user, heasarc_pwd, decryptkey):

    subprocess.call(['/bin/bash', '-i', '-c', 'heainit'])
    fname = "/homes/pipeline/logs/"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    with open(fname, 'w') as f:
        f.write("Completed at " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

    ### Running on PSR_B1937+21 ###
    os.chdir("/students/pipeline/heasoft6.27/PSR_B1937+21")
    par_1937 = "/students/pipeline/parfiles/PSR_B1937+21.par.nancaytzr"
    update("PSR_B1937+21", heasarc_user, heasarc_pwd, './', decryptkey, 
    par_1937)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("--obsID", help="ObsID for single psrpipe call", 
						default=None, type=str)
	parser.add_argument("--download", help="Run ni_datadownload", 
						default=False, action='store_true')
	parser.add_argument("--sourcename", help="source name for downloading",
						default=None, type=str)
	parser.add_argument("--user", help="heasarc username", 
						default='nicer_team', type=str)
	parser.add_argument("--passwd", help="heasarc password", 
						default='sextant', type=str)
	parser.add_argument("--outdir", help="output dir for download", 
						default="./", type=str)
	parser.add_argument("--clobber", help="overwrite download files", 
						default=False, action='store_true')
	parser.add_argument("-k", help="decryption key", dest='key',
						action=pipeline_utils.PasswordPromptAction, 
						type=str, required=False)
	parser.add_argument("--k_file", help="decryption key from file",
						type=str, default=None)
	parser.add_argument("--mp", default=False, action='store_true',
						help="Run all procedures with multiprocessing")
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
	
	parser.add_argument("--combine", help="combine evt/pha", 
						default=False, action='store_true')

	parser.add_argument("--update", default=None, type=str)
						help="Update single source with full procedure")
	parser.add_argument("--cron", help="Run cron job", 
						default=False, action='store_true')
	parser.add_argument("--trumpet", help="Run nicerl2 with trumpet cut",
						default=True, type=bool)
	parser.add_argument("--keith", help="Use standard filters from Keith Gendreau",
						default=True, type=bool)

	parser.add_argument("--max_date", help="Max date to merge",
						default=None, type=str)
	args = parser.parse_args()

    #Key handling
	if args.k_file is not None:
		with open(args.k_file, 'r') as f:
			args.key = f.readlines()[0].strip("\n")

	if args.download:
		assert(all([arg is not None for arg in [
					args.sourcename, args.user, args.passwd, args.key]]))
		pipeline_utils.run_datadownload(args.sourcename, 
										args.user, args.passwd, 
										args.outdir, 
										args.key, clobber=args.clobber)
	elif args.combine:
		if args.max_date is None:
			run_niextract(args.sourcename)
		else:
			dt = datetime.datetime.strptime(args.max_date, '%Y-%m-%d')
			run_niextract(args.sourcename, max_date=dt)
	elif args.update is not None:
		update(args.update, args.user, args.passwd, 
			   args.key, args.par, 
			   emin=args.emin, emax=args.emax, 
			   cut=args.cut, 
			   filterbinsize=args.filterbinsize, 
			   trumpet=args.trumpet,
			   keith=args.keith)
	elif args.cron:
		cronjob(args.user, args.passwd, args.key)

