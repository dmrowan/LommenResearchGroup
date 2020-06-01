#!/usr/bin/env python
from __future__ import print_function, division
from astropy.table import Table
from astropy import log
import datetime
import os
import argparse
import multiprocessing as mp
import subprocess
import getpass
import time
from datetime import datetime

#LommenResearchGroup imports
import pipeline_utils

#Dom Rowan 2020

desc = """
Calls psrpipe with standard options with mulitprocessing
"""

"""
Note about python versions:

Depending on package versions for various PINT dependencies, some of the psrpipe related code may throw syntax errors*. One way i solve this in the past was to jump to python 2.7.15. Therefore, this code does not use f strings so it is compatible with both python 2 and 3.

* the error is in some pint TOA thing
"""

#This is temporary until i get the par file for crab
def nicerl2_wrapper(clobber=False, trumpet=True):
	pool = mp.Pool(processes=mp.cpu_count()+2)
	jobs = []
	times = []
	for f in os.listdir("./"):
		if (os.path.isdir(f)) and (not f.endswith('_pipe')):
			if f != 'tmp':
				job = pool.apply_async(run_nicerl2, (f,), 
									   dict(trumpet=trumpet,clobber=clobber))
				jobs.append(job)
	for job in jobs:
		output = job.get()
		times.append(output)

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
		cmd.append('--par', par)
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
	
	if ((os.path.isfile("{}_pipe/cleanfilt.evt".format(obsID))) 
	    and (not clobber)):
		return 0
	elif ((os.path.isfile("{}_pipe/cleanfilt.evt".format(obsID)))
		  and clobber):
		log.info("Removing {}_pipe".format(obsID))
		os.rmdir("{}_pipe".format(obsID))

	if not os.path.isdir(obsID):
		log.error("obsID not found")
		return -1

	pipeline_utils.run_nicerl2(obsID, trumpet=trumpet)

	pipeline_utils.run_add_kp(obsID)

	run_psrpipe(obsID, par,
				emin=emin, emax=emax,
				keith=keith)
	return 1

def wrapper(par, emin=0.25, emax=12,
			trumpet=True, keith=True, clobber=False):

	pool = mp.Pool(processes=mp.cpu_count()+2)
	jobs = []
	times = []
	for f in os.listdir("./"):
		if (os.path.isdir(f)) and (not f.endswith('_pipe')):
			if f != 'tmp':
				job = pool.apply_async(allprocedures,
									   (f, par, ),
									   dict(emin=emin,
										    emax=emax,
											trumpet=trumpet,
											keith=keith,
											clobber=clobber))
				jobs.append(job)

	for job in jobs:
		output = job.get()
		times.append(output)


def update(sourcename, heasarc_user, heasarc_pwd, outdir, decryptkey,
		   par, emin=0.25, emax=12, cut=2, filterbinsize=16,
		   trumpet=True, keith=True):

	source_dir = os.getcwd().split('/')[-1]
	log.info("Checking sourcename")
	if sourcename != source_dir:
		check = outdircheck(sourcename)
		if check == -1:
			return

	pipeline_utils.run_datadownload(sourcename, heasarc_user, 
									heasarc_pwd, outdir, 
									decryptkey, clobber=False)

	for f in os.listdir(outdir):
		if (os.path.isdir(f)) and (not f.endswith('_pipe')):
			if not os.path.isdir(f+"_pipe"):
				allprocedures(f,  par, emin=emin, emax=emax,
							  trumpet=trumpet, keith=keith, 
							  clobber=False)

	run_niextract(sourcename)

	run_cr_cut(sourcename+"_combined.evt", cut, filterbinsize)

	
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
			date = datetime.strptime(tab.meta['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
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

#I want to change this to only perform cr_cut on merged event file
def run_cr_cut(merged_evt, cut, filterbinsize):
	log.info("Running cr_cut")

	cmd = ['cr_cut.py', 
		   '--cut', str(cut), 
		   '--filterbinsize', str(filterbinsize),
		   merged_evt]
	
	subprocess.call(cmd)

#Need to change this to work with decrypt key
def cronjob(heasarc_user, heasarc_pwd, decryptkey, emin, emax, mask,
			cormin, cut, filterbinsize, 
			trumpet, keith):

	fname = "/homes/pipeline/logs/"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
	with open(fname, 'w') as f:
		f.write("Completed at " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
		
	### Running on PSR_B1821-24 ###
	os.chdir("/students/pipeline/PSR_B1821-24")
	par_1821 = "/students/pipeline/parfiles/PSR_B1821-24.par"
	update("PSR_B1821-24", heasarc_user, heasarc_pwd, './', decryptkey, 
		   emin, emax, mask, par_1821, cormin, cut, 
		   filterbinsize, trumpet)

	### Running on PSR_B1937+21 ###
	os.chdir("/students/pipeline/PSR_B1937+21")
	par_1937 = "/students/pipeline/parfiles/PSR_B1937+21.par"
	update("PSR_B1937+21", heasarc_user, heasarc_pwd, './', decryptkey,
		   emin, emax, mask, par_1937, cormin, cut, 
		   filterbinsize, trumpet)

	### Running on PSR_J0218+4232 ###
	os.chdir("/students/pipeline/PSR_J0218+4232")
	par_0218 = "/students/pipeline/parfiles/PSR_J0218+4232.par"
	update("PSR_J0218+4232", heasarc_user, heasarc_pwd, './', decryptkey,
			emin, emax, mask, par_0218, cormin, cut, 
			filterbinsize, _trumpet)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("--obsID", help="ObsID for single psrpipe call", 
						default=None, type=str)
	parser.add_argument("--download", help="Run ni_datadownload", 
						default=False, action='store_true')
	parser.add_argument("--sourcename", help="source name for downloading",
						default=None, type=str)
	parser.add_argument("--user", help="heasarc username", 
						default=None, type=str)
	parser.add_argument("--passwd", help="heasarc password", 
						default=None, type=str)
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

	parser.add_argument("--update", default=False, action='store_true',
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
	elif args.mp:
		log.info("Launching wrapper")
		wrapper(args.par, emin=args.emin, emax=args.emax,  
				trumpet=args.trumpet, 
				keith=args.keith)
	elif args.combine:
		if args.max_date is None:
			run_niextract(args.sourcename)
		else:
			dt = datetime.strptime(args.max_date, '%Y-%m-%d')
			run_niextract(args.sourcename, max_date=dt)
	elif args.update:
		update(args.sourcename, args.user, args.passwd, 
			   args.outdir, args.key, args.par, 
			   emin=args.emin, emax=args.emax, 
			   cut=args.cut, 
			   filterbinsize=args.filterbinsize, 
			   trumpet=args.trumpet,
			   keith=args.keith)
	elif args.cron:
		cronjob(args.user, args.passwd, args.key, 
				args.emin, args.emax, args.mask,
				args.cormin, args.cut, args.filterbinsize, 
				args.no_trumpet, args.keith)


