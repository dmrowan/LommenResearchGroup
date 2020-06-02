#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
import argparse
from glob import glob
import subprocess 
import multiprocessing as mp
import os
import shutil

#Dom Rowan 2019

desc="""
Functions for processing and extracting background data
"""

#Define a function to check if nicerl2 has been run on an obsID
def check_nicerl2(obsID):
	event_cl_files = os.listdir(obsID+"/xti/event_cl")
	for mpu in [ "mpu"+str(i) for i in range(7) ]:
		if not any( [ mpu in f for f in event_cl_files ]):
			return 0
	return 1

def run_nicerl2(obsID, clobber=True,no_br_filter=True):
	log.info("Running nicerl2 for "+obsID)
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
	if no_br_filter:
		cmd.append("br_earth=0")
	
	subprocess.call(cmd)


def nicerl2_wrapper():

	#pool = mp.Pool(processes=mp.cpu_count()+2)
	pool = mp.Pool(processes=1)
	jobs = []

	for f in os.listdir("./"):
		if not f.endswith("_pipe"):
			if os.path.isdir(f) and (f != 'tmp') and (not check_nicerl2(f)):
				#job = pool.apply_async(run_nicerl2, (f,))
				#jobs.append(job)
				run_nicerl2(f)
	
	for job in jobs:
		job.get()

def bkg_pipe_wrapper():

	pool = mp.Pool(processes=mp.cpu_count()+2)
	jobs = []

	for f in os.listdir("./"):
		if not f.endswith("_pipe"):
			if os.path.isdir(f) and (f != 'tmp') and check_nicerl2(f):
				if not os.path.isfile(f+"_pipe/cleanfilt.evt"):
					bkg_pipe(f)

def run_add_kp(obsID):
	log.info("Running add_kp")

	cmd = ['add_kp.py', "{f}/auxil/ni{f}.mkf".format(f=obsID), '--kp', '/homes/pipeline/kp.fits']
	subprocess.call(cmd)


#Heavily based on Paul's bkgpipe.py
def bkg_pipe(obsID):

	assert(check_nicerl2(obsID))

	run_add_kp(obsID)

	#Define some constants
	minfpm = 38

	pipedir = "{0}_pipe".format(obsID)
	if not os.path.exists(pipedir):
		os.makedirs(pipedir)
	
	mkf = glob(os.path.join(obsID, 'auxil/ni*.mkf'))[0]
	log.info('MKF File: {0}'.format(mkf))
	shutil.copy(mkf,pipedir)

	att = glob(os.path.join(obsID, 'auxil/ni*.att'))[0]
	log.info('ATT File: {0}'.format(att))
	shutil.copy(att,pipedir)

	orb = glob(os.path.join(obsID, 'auxil/ni*.orb'))[0]
	log.info('ORB File: {0}'.format(orb))
	shutil.copy(orb,pipedir)
	
	bkf = os.path.join(pipedir, obsID)+'_prefilt.bkf'
	log.info('BKF File: {0}'.format(bkf))



	gti = os.path.join(pipedir, 'tot.gti')

	overonly_string = ['FPM_OVERONLY_COUNT<1',
					   'FPM_OVERONLY_COUNT<(1.52*COR_SAX**(-0.633))']
	cor_sax_string = ['(COR_SAX.gt.(1.914*KP**0.684+0.25))']
	kp_string = ['KP.lt.5']
	sunshine_string = ['(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)'.format(60)]
	extra_expr = overonly_string+cor_sax_string+kp_string+sunshine_string
	extra_expr = "("+" && ".join("%s" %expr for expr in extra_expr) + ")"

	maxunder=200.0

	cmd = ['nimaketime', 'infile={0}'.format(mkf),
		   'outfile={0}'.format(gti), 'nicersaafilt=YES',
		   'saafilt=NO', 'trackfilt=YES', 'ang_dist=0.015', 'elv=20',
		   'br_earth=0', 'min_fpm={0}'.format(minfpm), 
		   'underonly_range=0-{0}'.format(maxunder),
		   'expr={0}'.format(extra_expr), 
		   'outexprfile={0}'.format(
				   os.path.join(pipedir, "bkgpipe_expr.txt")), 
		   'clobber=YES']
	
	subprocess.call(cmd)

	evfiles = glob(os.path.join(obsID,'xti/event_cl/ni*mpu7_cl.evt'))
	evlistname=os.path.join(pipedir, 'evfiles.txt')
	with open(evlistname, 'w') as f:
		for ev in evfiles:
			f.write(ev+'\n')
	evfiles.sort()
	log.info("Cleaned Event Files: {0}".format(
				"\n		\n".join(evfiles)))

	
	evtfilename = os.path.join(pipedir, 'cleanfilt.evt')
	cmd = ["niextract-events", 
		   'filename=@{0}[EVENT_FLAGS=bx1x000]'.format(evlistname),
		   "eventsout={0}".format(evtfilename), 
		   "timefile={0}".format(gti), "gti=GTI", "clobber=yes"]

	subprocess.call(cmd)
	
	##Final step filters the mkf files using the GTI
	gti_expr = 'gtifilter("{0}")'.format(gti)

	filtered_mkf = os.path.join(pipedir, "mkf_filtered.mkf")
	cmd = ["ftcopy", '{0}[{1}]'.format(mkf, gti_expr), filtered_mkf, 
		   "clobber=yes", "history=yes"]
	subprocess.call(cmd)

def niextract_all_events(verbose=True):
	region_nums = [1,2,3,4,5,6,8]
	region_names = [ "BKGD_RXTE_{0}".format(i) for i in region_nums ]
	for d in region_names:
		assert(os.path.isdir(d))

	filenames = []
	mkfnames = []
	for d in region_names:
		for f in os.listdir(d):
			if os.path.isfile(d+"/"+f+"/cleanfilt.evt"):
				filenames.append(d+"/"+f+"/cleanfilt.evt")
				mkfnames.append(d+"/"+f+"/mkf_filtered.mkf")
	
	with open("bkgd_all_evt_list", 'w') as f:
		for fn in filenames:
			f.write(fn+"\n")
	with open("bkgd_all_mkf_list", 'w') as f1:
		for fn in mkfnames:
			f1.write(fn+"\n")

	output='bkgd_combined.evt'
	cmd = ['niextract-events', 'filename=@bkgd_all_evt_list',
		   'eventsout={}'.format(output), 'clobber=yes']
	if verbose: cmd.append('chatter=5')
	subprocess.call(cmd)

	cmd = ['ftmerge', 'infile=@bkgd_all_mkf_list',
		   'outfile=bkgd_merged.mkf', 'clobber=YES']
	if verbose: cmd.append('chatter=5')
	subprocess.call(cmd)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("--download", help="nidata download for input source", default=None)
	parser.add_argument("--nicerl2", help="run nicerl2 on current directory", 
						default=False, action='store_true')
	parser.add_argument("--pipe", help="run pipe on current directory", 
						default=False, action='store_true')
	parser.add_argument("--extract", help="extract all events", default=False, action='store_true')
	args = parser.parse_args()

	if args.download is not None:
		print("I haven't added this yet")
	
	elif args.nicerl2:
		nicerl2_wrapper()

	elif args.pipe:
		bkg_pipe_wrapper()

	elif args.extract:
		niextract_all_events()


