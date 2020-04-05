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

desc = """
Calls psrpipe with standard options with mulitprocessing
"""

class PasswordPromptAction(argparse.Action):
	def __init__(self, option_strings, dest=None, nargs=0, 
				 default=None, required=False, type=None, 
				 metavar=None, help=None):
		super(PasswordPromptAction, self).__init__(
				option_strings=option_strings, 
				dest=dest,
				nargs=nargs,
				default=default,
				required=required,
				metavar=metavar,
				type=type,
				help=help)
	
	def __call__(self, parser, args, values, option_string=None):
		password = getpass.getpass()
		setattr(args, self.dest, password)


def outdircheck(sourcename):
	print("WARNING: source name {}\
		   inconsistent with output dir".format(sourcename))
	check = raw_input("Continue -- [y/n] ")
	if any([check == ['y', 'Y']]):
		return -1
	else:
		return 0

def run_datadownload(sourcename, heasarc_user, heasarc_pwd, outdir,
					 decryptkey, clobber=False, unzip=True):
	if outdir == "./":
		if not os.getcwd().endswith(sourcename):
			if outdircheck(sourcename) == -1:
				return -1
	else:
		if not (outdir.endswith(sourcename or sourcename+"/")):
			if outdircheck(soucename) == -1:
				return -1
			
	cmd = ['ni_data_download.py', 
		   sourcename, heasarc_user, heasarc_pwd,
		   '--outdir', outdir,
		   '--decryptkey', decryptkey]
	if clobber:
		cmd.append('--clobber')
	if unzip:
		cmd.append('--unzip')

	subprocess.call(cmd)

def run_nicerl2(obsID, clobber=True, no_trumpet=False):
	log.info("Running nicerl2")
	assert(os.path.isdir(obsID))
	
	if not os.path.isdir("tmp/"+obsID+"_pfiles"):
		os.mkdir("tmp/"+obsID+"_pfiles")
	
	abspath = os.path.abspath("tmp/"+obsID+"_pfiles")
	os.environ['PFILES'] = ( abspath + 
							';/packages/heasoft-6.26.1/x86_64-pc-linux-gnu-libc2.23/syspfiles')

	log.info("Set pfiles to" + os.environ['PFILES'])

	cmd = ['nicerl2', obsID]
	if clobber:
		cmd.append("clobber=YES")
	if no_trumpet:
		log.info("Not using trumpet filter")
		cmd.append("trumpetfilt=NO")
	subprocess.call(cmd)

#no need to run niprefilter2 explicitly because its included in nicerl2
def run_niprefilter2(obsID):
	log.info("Running niprefilter2")
	infile = "infile=./{f}/auxil/ni{f}.mkf".format(f=obsID)
	outfile = "outfile=./{f}/auxil/ni{f}.mkf".format(f=obsID)
	cmd = ['niprefilter2', 'indir=./{}'.format(obsID),
		   infile, outfile,
		   'clobber=YES']

	subprocess.call(cmd)

def run_add_kp(obsID):
	log.info("Running add_kp")
	
	cmd = ['add_kp.py', "{f}/auxil/ni{f}.mkf".format(f=obsID)]
	print(cmd)
	subprocess.call(cmd)

def run_psrpipe(obsID, emin, emax, mask,
				par, cormin, 
				keith=False):
	log.info("Running pspipe with par " + par)
	start = time.time()
	if not os.path.isdir(obsID):
		print("obsID not found")
		return -1
	
	cmd = ['psrpipe.py', 
		   '--par', par,
		   #'--mask', str(mask),
		   '--emin', str(emin), 
		   '--emax', str(emax),
		   '--minsun', str(60)]
	if keith:
		cmd.append('--keith')

	cmd.append(obsID)

	print(cmd)
	subprocess.call(cmd)

	end = time.time()
	elapsed = "{:10.4f}".format(end-start)
	return "{}: {} seconds".format(obsID, elapsed)

def run_cr_cut(obsID, cut, filterbinsize):
	log.info("Running cr_cut")
	evtfile = "{}_pipe/cleanfilt.evt".format(obsID)
	lcfile = "{}_pipe/cleanfilt.lc".format(obsID)
	if (not os.path.isfile(evtfile)) or (not os.path.isfile(lcfile)):
		log.info("cr_cut error: evt or lcfile not found")
	cmd = ['cr_cut.py', 
		   #'--lcfile', lcfile,
		   '--cut', str(cut), 
		   '--filterbinsize', str(filterbinsize),
		   evtfile]
	
	subprocess.call(cmd)

def allprocedures(obsID, emin, emax, 
				  mask, par, cormin, 
				  cut, filterbinsize,
				  no_trumpet=False, keith=False):
	
	if os.path.isfile("{}_pipe/cleanfilt_cut.evt".format(obsID)):
		return 0

	start = time.time()
	if not os.path.isdir(obsID):
		print("obsID not found")
		return -1
	

	#run_nicerl2(obsID, no_trumpet=no_trumpet)
	#run_niprefilter2(obsID)
	run_add_kp(obsID)
	run_psrpipe(obsID, emin, emax, mask, 
				par, cormin, 
				keith=keith)
	run_cr_cut(obsID, cut, filterbinsize)
	end = time.time()
	elapsed = "{:10.4f}".format(end-start)
	return "{}: {} seconds".format(obsID, elapsed)

def wrapper(emin, emax, mask, par, cormin, cut, filterbinsize,
			output="OutputTimes.txt",
			no_trumpet=False, keith=False):

	#pool = mp.Pool(processes=1)
	pool = mp.Pool(processes=mp.cpu_count()+2)
	jobs = []
	times = []
	for f in os.listdir("./"):
		if (os.path.isdir(f)) and (not f.endswith('_pipe')):
			if f != 'tmp':
				job = pool.apply_async(allprocedures,
									   (f, emin, emax, mask, par, cormin,
											 cut, filterbinsize,),
									   dict(no_trumpet=no_trumpet,
											keith=keith))
				jobs.append(job)
	for job in jobs:
		output = job.get()
		times.append(output)


def update(sourcename, heasarc_user, heasarc_pwd, outdir, decryptkey,
		   emin, emax, mask, par, cormin, cut, filterbinsize,
		   no_trumpet, keith):
	source_dir = os.getcwd().split('/')[-1]
	log.info("Checking sourcename")
	if sourcename != source_dir:
		check = outdircheck(sourcename)
		if check == -1:
			return
	run_datadownload(sourcename, heasarc_user, heasarc_pwd, outdir, 
					 decryptkey, clobber=False, unzip=True)
	for f in os.listdir(outdir):
		if (os.path.isdir(f)) and (not f.endswith('_pipe')):
			if not os.path.isdir(f+"_pipe"):
				allprocedures(f, emin, emax, mask, par, cormin, 
							  cut, filterbinsize,
							  no_trumpet=no_trumpet, keith=keith)
	run_niextract(sourcename)
	#merge_spectra(sourcename)

	
def run_niextract(sourcename, max_date=None):
	source_dir = os.getcwd().split('/')[-1]
	if sourcename != source_dir:
		check = outdircheck(sourcename)
		if check == -1:
			return
	filenames = []
	invalid_columns = []
	for f in os.listdir("./"):
		if os.path.isfile(f+"/cleanfilt_cut.evt"):
			tab = Table.read(f+"/cleanfilt_cut.evt", hdu=1)
			date = datetime.strptime(tab.meta['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
			if max_date is not None:
				if date > max_date:
					print("Excluding "+f+" at date "+str(date))
					continue

			if len(tab.columns) == 15:
				filenames.append(f+"/cleanfilt_cut.evt")
			else:
				invalid_columns.append(f)
	#print(invalid_columns)
	print(len(invalid_columns))
	with open(sourcename+"_evtlist", 'w') as f1:
		for fn in filenames:
			f1.write(fn+"\n")
	output = sourcename+"_combined.evt"
	cmd = ['niextract-events', 'filename=@{}_evtlist'.format(sourcename),
		   'eventsout={}'.format(output), 'clobber=yes']
	subprocess.call(cmd)

	for f in invalid_columns:
		print(f, "not removing")
		subprocess.call(['rm', '-rf', f])

def merge_spectra(sourcename):
	filenames = []
	for f in os.listdir("./"):
		if os.path.isfile(f+"_pipe/cleanfilt.pha"):
			filenames.append(f+"_pipe/cleanfilt.pha")
	with open(sourcename+"_phalist", 'w') as f:
		for fn in filenames:
			f.write(fn+"\n")
	cmd = ['addspec', 'infil={}_phalist'.format(sourcename), 
		   'outfil=combined.pha']
	subprocess.call(cmd)

#Need to change this to work with decrypt key
def cronjob(heasarc_user, heasarc_pwd, decryptkey, emin, emax, mask,
			cormin, cut, filterbinsize, 
			no_trumpet, keith):

	fname = "/homes/pipeline/logs/"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
	with open(fname, 'w') as f:
		f.write("Completed at " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
		
	### Running on PSR_B1821-24 ###
	os.chdir("/students/pipeline/PSR_B1821-24")
	par_1821 = "/students/pipeline/parfiles/PSR_B1821-24.par"
	update("PSR_B1821-24", heasarc_user, heasarc_pwd, './', decryptkey, 
		   emin, emax, mask, par_1821, cormin, cut, 
		   filterbinsize, no_trumpet)

	### Running on PSR_B1937+21 ###
	os.chdir("/students/pipeline/PSR_B1937+21")
	par_1937 = "/students/pipeline/parfiles/PSR_B1937+21.par"
	update("PSR_B1937+21", heasarc_user, heasarc_pwd, './', decryptkey,
		   emin, emax, mask, par_1937, cormin, cut, 
		   filterbinsize, no_trumpet)

	### Running on PSR_J0218+4232 ###
	os.chdir("/students/pipeline/PSR_J0218+4232")
	par_0218 = "/students/pipeline/parfiles/PSR_J0218+4232.par"
	update("PSR_J0218+4232", heasarc_user, heasarc_pwd, './', decryptkey,
			emin, emax, mask, par_0218, cormin, cut, 
			filterbinsize, no_trumpet)

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
						action=PasswordPromptAction, 
						type=str, required=False)
	parser.add_argument("--k_file", help="decryption key from file",
						type=str, default=None)
	parser.add_argument("--mp", default=False, action='store_true',
						help="Run all procedures with multiprocessing")
	parser.add_argument("--emin", help="Minimum energy to include",
						type=float, default=0.25)
	parser.add_argument("--emax", help="Maximum energy to include", 
						type=float, default=12.0)
	parser.add_argument("--mask", help="Mask these IDS", type=int, default=-1)
	parser.add_argument("--cormin", 
						help="Set minimum cutoff rigidity for nimaketime", 
						type=float, default=3.0)
	parser.add_argument("--par", help="Par file to use for phases", 
						type=str, default="")
	parser.add_argument("--cut", help="Count rate for cut in cts/sec", 
						default=2.0, type=float)
	parser.add_argument("--filterbinsize", help="Bin size in seconds", 
						default=8.0, type=float)
	
	parser.add_argument("--combine", help="combine evt/pha", 
						default=False, action='store_true')

	parser.add_argument("--update", default=False, action='store_true',
						help="Update single source with full procedure")
	parser.add_argument("--cron", help="Run cron job", 
						default=False, action='store_true')
	parser.add_argument("--no_trumpet", help="Run nicerl2 without trumpet cut",
						default=False, action='store_true')
	parser.add_argument("--keith", help="Use standard filters from Keith Gendreau",
						default=False, action='store_true')

	parser.add_argument("--max_date", help="Max date to merge",
						default=None, type=str)
	args = parser.parse_args()

	if args.k_file is not None:
		with open(args.k_file, 'r') as f:
			args.key = f.readlines()[0].strip("\n")

	if args.download:
		assert(all([arg is not None for arg in [
					args.sourcename, args.user, args.passwd, args.key]]))
		run_datadownload(args.sourcename, args.user, args.passwd, args.outdir, 
						 args.key, clobber=args.clobber, unzip=True)
	elif args.mp:
		log.info("Launching wrapper")
		wrapper(args.emin, args.emax, args.mask, args.par, args.cormin, 
				args.cut, args.filterbinsize, 
				output="OutputTimes.txt", 
				no_trumpet=args.no_trumpet, 
				keith=args.keith)
	elif args.combine:
		if args.max_date is None:
			run_niextract(args.sourcename)
		else:
			dt = datetime.strptime(args.max_date, '%Y-%m-%d')
			run_niextract(args.sourcename, max_date=dt)
	elif args.update:
		update(args.sourcename, args.user, args.passwd, 
			   args.outdir, args.key, args.emin, args.emax, 
			   args.mask, args.par, args.cormin, args.cut, 
			   args.filterbinsize, 
			   args.no_trumpet,
			   args.keith)
	elif args.cron:
		cronjob(args.user, args.passwd, args.key, 
				args.emin, args.emax, args.mask,
				args.cormin, args.cut, args.filterbinsize, 
				args.no_trumpet, args.keith)


	"""

	for d in os.listdir("./"):
		if os.path.isdir(d):
			if d.endswith("_pipe"):
				if os.path.isfile(d+"/cleanfilt.evt") and not os.path.isfile(d+"/cleanfilt_cut.evt"):
					cmd = ['rm', '-rf', d]
					obsID = d[:-5]
					print(obsID)
					subprocess.call(cmd)
					allprocedures(obsID, args.emin, args.emax, args.mask, 
								  args.par, args.cormin, 
								  args.cut, args.filterbinsize)

	allprocedures("1070030241", args.emin, args.emax, args.mask, 
				  args.par, args.cormin, args.cut, args.filterbinsize, 
				  keith=True)
	"""
	

