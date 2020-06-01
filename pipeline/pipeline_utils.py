#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
import argparse
import subprocess
import os

#Dom Rowan 2020

desc="""
Pipeline utilities for new summer 2020 system
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


#Define a function to check if nicerl2 has been run on an obsID
def check_nicerl2(obsID):
	event_cl_files = os.listdir(obsID+"/xti/event_cl")
	for mpu in [ "mpu"+str(i) for i in range(7) ]:
		if not any( [ mpu in f for f in event_cl_files ]):
			return 0
	return 1

#Command line call of nicerl2
def run_nicerl2(obsID, clobber=False, br_filter=True, horizon_filter=False,
				trumpet=True):
	'''
	if br_filter is True, the bright earth filter is used with default value
	else br_earth=0

	if horizon_filter=True multiple conditions are changed for horizon
		crossing project
	else default values are used
	'''

	assert(os.path.isdir(obsID))

	log.info("Running nicerl2 for "+obsID)

	if not os.path.isdir("tmp/"):
		log.error("No tmp directory for pfiles exists, exiting")
		return 0

	if not os.path.isdir("tmp/"+obsID+"_pfiles"):
		log.warning("Creating pfile dir for "+obsID)
		os.mkdir("tmp/"+obsID+"_pfiles")

	abspath = os.path.abspath("tmp/"+obsID+"_pfiles")
	os.environ['PFILES'] = ( abspath +
					        ';/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23/syspfiles')

	log.info("Set pfiles to" + os.environ['PFILES'])

	cmd = ['nicerl2', obsID]
	if clobber:
		cmd.append("clobber=YES")
	

	if not br_filter:
		cmd.append("br_earth=0")

	if horizon_filter:
		cmd.append('elv=0')
		cmd.append('br_earth=0')
		cmd.append('nicersaafilt=NO')
		cmd.append('trackfilt=NO')
		cmd.append('st_valid=NO')
		cmd.append('ang_dist=0.015')

	if not trumpet:
		log.info("Not using trumpet filter")
		cmd.append("trumpetfilt=NO")

	subprocess.call(cmd)


def run_add_kp(obsID):
	log.info("Running add_kp with potsdam values")
	cmd = ['add_kp.py', "{f}/auxil/ni{f}.mkf".format(f=obsID),
		   '--potsdam']


def outdircheck(sourcename):
	log.warning("WARNING: source name {}\
	inconsistent with output dir".format(sourcename))
	check = raw_input("Continue -- [y/n] ")
	if any([check == ['y', 'Y']]):
		return 1
	else:
		return 0

def run_datadownload(sourcename, heasarc_user, heasarc_pwd, outdir,
					 decryptkey, clobber=False):

	if outdir == "./":
		if not os.getcwd().endswith(sourcename):
			if not outdircheck(sourcename):
				return 0
	else:
		if not (outdir.endswith(sourcename or sourcename+"/")):
			if not outdircheck(soucename):
				return 0

	cmd = ['ni_data_download.py', sourcename,
		   heasarc_user, heasarc_pwd, '--outdir', outdir,
		   '--decryptkey', decryptkey, '--unzip']

	if clobber:
		cmd.append('--clobber')

	subprocess.call(cmd)

