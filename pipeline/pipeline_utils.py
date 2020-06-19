#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
from astropy.time import Time
import argparse
from bs4 import BeautifulSoup
import datetime
import getpass
import numpy as np
import pandas as pd
import pint
import requests
import subprocess
import os

#LommenResearchGroup imports
import niutils

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
        log.error("No tmp directory for pfiles exists")
        response = input("Would you like to create a tmp directory? [y/n] ")
        if response in ['Y', 'y']:
            os.mkdir('tmp')
            log.info("tmp directory created")
        else:
            log.error("Exiting")
            return 0

    if not os.path.isdir("tmp/"+obsID+"_pfiles"):
        log.warning("Creating pfile dir for "+obsID)
        os.mkdir("tmp/"+obsID+"_pfiles")

    #Set environmental variables for pfiles
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
        cmd.append('ang_dist=180')

    if not trumpet:
        log.info("Not using trumpet filter")
        cmd.append("trumpetfilt=NO")

    subprocess.call(cmd)

def run_add_kp(obsID):
	log.info("Running add_kp with potsdam values")
	cmd = ['add_kp.py', "{f}/auxil/ni{f}.mkf".format(f=obsID),
		   '--potsdam']

	subprocess.call(cmd)


def outdircheck(sourcename):
	log.warning("source name {}\
	inconsistent with output dir".format(sourcename))
	check = input("Continue -- [y/n] ")
	if check in ['y', 'Y']:
		return True
	else:
		return False

def run_datadownload(sourcename, heasarc_user, heasarc_pwd, outdir,
					 decryptkey, clobber=False, obsIDs=None):

	if outdir == "./":
		if not os.getcwd().endswith(sourcename):
			if not outdircheck(sourcename):
				return 0
	else:
		if not (outdir.endswith(sourcename or sourcename+"/")):
			if not outdircheck(soucename):
				return 0

	cmd = ['custom_data_download.py', sourcename,
		   heasarc_user, heasarc_pwd, '--outdir', outdir,
		   '--decryptkey', decryptkey, '--unzip']

	if clobber:
		cmd.append('--clobber')

	if obsIDs is not None:
		cmd.append('--obsIDs')
		cmd.extend(obsIDs)

	subprocess.call(cmd)

#For large files the curl fail3
#This breaks the download into chunks
def download_parts(obsid, url, cleanup=False):
    
    log.info("Attempting to curl {0} in parts".format(obsid))

    #Initial byte range
    byte_range = [0, 999999999]
    fname_i = 0
    fname = '{0}.tar.part{1}'.format(obsid, fname_i)

    cmd = ['curl', '--retry', '10', url, '--ftp-ssl', '-k', 
           '--create-dirs', '-o', fname, '--range',
           '{0}-{1}'.format(str(byte_range[0]), str(byte_range[1]))]

    subprocess.call(cmd)

    while os.path.getsize(fname) == 1e9:
        byte_range = [ b+1e9 for b in byte_range ]
        fname_i += 1
        fname = '{0}.tar.part{1}'.format(obsid, fname_i)
        cmd = ['curl', '--retry', '10', url, '--ftp-ssl', '-k', 
               '--create-dirs', '-o', fname, '--range',
               '{0}-{1}'.format(str(int(byte_range[0])), str(int(byte_range[1])))]

        subprocess.call(cmd)

    log.info("Data downloaded in {} parts".format(fname_i+1))

    #Merge into 1 tar file
    #subprocess.call isn't working here so I'm using os.system
    cmd = 'cat {0}.tar.part* > {0}.tar'.format(obsid)
    os.system(cmd)

    assert(os.path.isfile('{0}.tar'.format(obsid)))
    log.info('tar merged for {0}'.format(obsid))

    if cleanup:
        log.info('cleaning partial tars')
        for i in range(fname_i+1):
            os.remove('{0}.tar.part{1}'.format(obsid, i))

    return 0 

def version_check():
	print("Heasoft version: ", 
		  subprocess.check_output('echo $HEADAS', shell=True))

	print("Last NICERsift git pull",
		  subprocess.check_output(
			  'stat -c %y /homes/pipeline/nicersoft/.git/FETCH_HEAD', 
			  shell=True))

	print("pint version: ", pint.__version__)


#Code from ni_data_download
#Can't do an import because code is written as executable for some reason...
def print_nicer_segment(url = 'https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html',
                        username = None, password=None):
    """
    This prints out the segment detail table in text format
    usage: % print_nicer_segment(username = "nicer_user_name" password = "nicer password")
    outputs: prints the nicer segment table to the terminal
    :param url: location of the segment detail page
    :param username: nicer team username
    :param password: nicer team password
    :return:
    """
    if (not username) or (not password):
        raise ValueError("must supply username and password to access the NICER obs page")
    req = requests.get(url, auth=(username, password))
    if req.status_code != 200:
        raise ValueError('Problem accessing {0} with ({1}, {2}) \nReturn code: {3}'.format(
            url, username, password, req.status_code))
        
    soup = BeautifulSoup(req.text, 'lxml')
    tabs = soup.find_all('table')[1]    
    df = pd.read_html(str(tabs))
    return df[0]


def get_exposure(sources):
	df = print_nicer_segment(username='nicer_team', password='sextant')

	if not niutils.check_iter(sources):
		sources = [sources]
	else:
		pass

	exp_list = []
	n_obsIDs = []
	for s in sources:
		df_selection = df[df['Target Name']==s]
		n_obsIDs.append(len(df_selection))
		exp_list.append(np.sum(df_selection['Good Expo[s]']))

	return exp_list, n_obsIDs


def crab_par_table(par_dir='/students/pipeline/parfiles/crab/'):
    
    par_files = os.listdir(par_dir)
    par_files = [ os.path.join(par_dir, p) for p in par_files ]
    
    start = []
    finish = []
    for par in par_files:
        with open(par, 'r') as f:
            lines = f.readlines()

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
		
    df = pd.DataFrame({'par':par_files, 'start':start, 'finish':finish})
    return df

#Manages product backups for evt and mkfs
def product_backup(fname, backupdir=None, date='now', message=""):

    #Are we working with evt or mkf:
    extension = os.path.splitext(fname)[1].lstrip('.')
    if extension == 'evt':
        mode='event'
    elif extension == 'mkf':
        mode='filter'
    else:
        raise TypeError("input file must be evt or mkf")

    #Use extension to determine backup directory
    if backupdir is None:
        backupdir = '{}_backups'.format(extension)

    
    #First do file and dir checks
    if not os.path.isfile(fname):
        raise FileNotFoundError("{} file not found".format(mode))

    #Option to create backupdir if it doesn't exist
    if not os.path.isdir(backupdir):
        log.warning("Backup directory doesn't exist")
        check = input("Would you like to create backup directory [y/n] -- ")
        if check in ['Y', 'y']:
            os.mkdir(backupdir)
        else:
            log.error("Exiting")
            return -1

    #Collect date information
    if date=='now':
        backup_date = datetime.datetime.now()
    elif date=='modified':
        backup_date = datetime.datetime.fromtimestamp(os.path.getmtime(fname))
    else:
        raise ValueError("Invalid date mode {}".format(date))
    #Format date
    format_specifier = '%Y-%m-%d_%H:%M:%S'
    backup_date = backup_date.strftime(format_specifier)
   
    #If the backup log doesn't exist, create and write first few lines
    backup_log = os.path.join(backupdir, 'log.txt')
    if not os.path.isfile(backup_log):
        with open(backup_log, 'w') as f:
            f.write("{} file backups\n".format(mode))
            f.write("Filename\tDate Created\tMessage\n")

    #Find the most recent backup in log
    with open(backup_log, 'r') as f:
        lines = f.readlines()
        last = lines[-1].split()[0]

    #Set backup name
    if last == 'Filename':
        backup_name = os.path.join(backupdir, 'backup_0.{}'.format(extension))
    else:
        new_num = str(int(last.strip('.{}'.format(extension))[-1])+1)
        backup_name = os.path.join(backupdir, 'backup_{0}.{1}'.format(
                new_num, extension))
    
    if os.path.isfile(backup_name):
        raise FileExistsError("Backup name already exists. This means that the log file and backup directory are out of sync")

    #Copy evt file to backupdir
    cmd = ['cp', fname, backup_name]
    subprocess.call(cmd)


    #Add info to log
    with open(backup_log, 'a') as f:
        f.write('\n')
        f.write('{0}\t\t{1}\t\t{2}\t\t'.format(
            os.path.basename(os.path.normpath(backup_name)),
            backup_date,
            message))

    log.info("Backup {} file created".format(mode))

