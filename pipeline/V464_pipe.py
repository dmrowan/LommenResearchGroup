#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
import argparse
from glob import glob
import subprocess 
import multiprocessing as mp
import os
import shutil

import pipeline_utils

#Dom Rowan 2020

desc="""
Pipeline procedure for V4641_Sgr
This is the source used for horizon crossing project
"""

def V464_pipe(obsID, clobber=False):

	if (clobber) or (not pipeline_utils.check_nicerl2(obsID)):
		pipeline_utils.run_nicerl2(obsID, clobber=True, horizon_filter=True)

	#Define some constants
	minfpm = 7
	maxunder=200.0

	pipedir = "{0}_pipe".format(obsID)

	if not os.path.exists(pipedir):
		os.makedirs(pipedir)
	elif clobber:
		shutil.rmtree(pipedir)
		os.makedirs(pipedir)
	else:
		log.info("ObsID already cleaned, exiting")
		return 1
	
	pipeline_utils.run_add_kp(obsID)
	mkf = glob(os.path.join(obsID, 'auxil/ni*.mkf'))[0]
	log.info('MKF File: {0}'.format(mkf))
	shutil.copy(mkf,pipedir)

	att = glob(os.path.join(obsID, 'auxil/ni*.att'))[0]
	log.info('ATT File: {0}'.format(att))
	shutil.copy(att,pipedir)

	orb = glob(os.path.join(obsID, 'auxil/ni*.orb'))[0]
	log.info('ORB File: {0}'.format(orb))
	shutil.copy(orb,pipedir)
	
	gti = os.path.join(pipedir, 'tot.gti')

	#Adding expressions not defined in nicerl2
	cor_sax_string = ['(COR_SAX.gt.(1.914*KP**0.684+0.25))']
	kp_string = ['KP.lt.5']
	extra_expr = cor_sax_string+kp_string
	extra_expr = "("+" && ".join("%s" %expr for expr in extra_expr) + ")"


	cmd = ['nimaketime', 'infile={0}'.format(mkf),
		   'outfile={0}'.format(gti), 'nicersaafilt=NO',
		   'saafilt=NO', 'trackfilt=NO', 'ang_dist=180', 'elv=0',
		   'st_valid=NO',
		   'br_earth=0', 'min_fpm={0}'.format(minfpm), 
		   'underonly_range=0-{0}'.format(maxunder),
		   'expr={0}'.format(extra_expr), 
		   'outexprfile={0}'.format(
				   os.path.join(pipedir, "V464pipe_expr.txt")), 
		   'clobber=YES', 'chatter=5']
	
	subprocess.call(cmd)

	#Grab event files
	evfiles = glob(os.path.join(obsID,'xti/event_cl/ni*mpu7_cl.evt'))
	evlistname=os.path.join(pipedir, 'evfiles.txt')
	with open(evlistname, 'w') as f:
		for ev in evfiles:
			f.write(ev+'\n')
	evfiles.sort()
	log.info("Cleaned Event Files: {0}".format(
				"\n		\n".join(evfiles)))

	
	#Run niextract events to make cleanfilt.evt
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

def horizon_crossing_wrapper(user, passwd, outdir, key, clobber=False):

	if clobber:
		for obs in ['2200300101', '2200300102']:
			if os.path.isdir(obs):
				shutil.rmtree(obs)
		
			if os.path.isdir(obs+"_pipe"):
				shutil.rmtree(obs+"_pipe")

	pipeline_utils.run_datadownload(
			'V4641_Sgr',user, passwd, outdir,
			key, clobber=clobber, 
			obsIDs=['2200300101', '2200300102'])

	V464_pipe('2200300101')
	V464_pipe('2200300102')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("--download", help="Download data for source",
						default=False, action='store_true')
	parser.add_argument("--user", help="username for data download",
						default=None, type=str)
	parser.add_argument("--passwd", help="passwd for data download",
						default=None, type=str)
	parser.add_argument("--outdir", help="outdir for data download",
						default='./', type=str)
	parser.add_argument("-k", help='decryption key', dest='key',
						action=pipeline_utils.PasswordPromptAction,
						type=str, required=False)
	parser.add_argument("--k_file", help='decryption key from file',
						type=str, default=None)
	parser.add_argument("--clobber", help="overwrite existing",
						default=False, action='store_true')
	parser.add_argument("--horizon", 
						help="re-run for horizon crossing project",
						default=False, action='store_true')
	parser.add_argument("--trackfilt",
						help="when re-running for horizon crossing project",
						default=False, action='store_true')
						

	args = parser.parse_args()

	if args.k_file is not None:
		with open(args.k_file, 'r') as f:
			args.key = f.readlines()[0].strip('\n')

	if args.download:
		pipeline_utils.run_datadownload(
				'V4641_Sgr',args.user, args.passwd, args.outdir,
				args.key, clobber=args.clobber, 
				obsIDs=['2200300101', '2200300102'])

	
	if args.horizon:
		horizon_crossing_wrapper(args.user, args.passwd, args.outdir, 
								 args.key, clobber=args.clobber)

