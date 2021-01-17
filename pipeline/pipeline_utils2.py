#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
from astropy.table import Table
import argparse
from bs4 import BeautifulSoup
import datetime
import getpass
import multiprocessing as mp
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
		password = 'sextant'
                #getpass.getpass()
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

    print('hello17')
    if not os.path.isdir("tmp/"):
        log.error("No tmp directory for pfiles exists")
        response = input("Would you like to create a tmp directory? [y/n] ")
        if response in ['Y', 'y']:
            os.mkdir('tmp')
            log.info("tmp directory created")
        else:
            log.error("Exiting")
            return 0

    print('hello18')
    if not os.path.isdir("tmp/"+obsID+"_pfiles"):
        log.warning("Creating pfile dir for "+obsID)
        os.mkdir("tmp/"+obsID+"_pfiles")

    #Set environmental variables for pfiles
    abspath = os.path.abspath("tmp/"+obsID+"_pfiles")
    os.environ['PFILES'] = ( abspath +
        ';/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23/syspfiles')

    log.info("Set pfiles to" + os.environ['PFILES'])

    print('hello19')
    cmd = ['nicerl2', obsID]
    if clobber:
        cmd.append("clobber=YES")

    if not br_filter:
        cmd.append("br_earth=0")

    print('hello20')
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
	check = 'y'
	if check in ['y', 'Y']:
		return True
	else:
		return False

def run_datadownload(sourcename, heasarc_user, heasarc_pwd, outdir,
					 decryptkey, clobber=False, obsIDs=None,
                     silent_curl=False):

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
        
        if not niutils.check_iter(obsIDs):
            obsIDs = [obsIDs]

        cmd.append('--obsIDs')
        cmd.extend(obsIDs)
    if silent_curl:
        cmd.append("--silent_curl")

    subprocess.call(cmd)

#For large files the curl fails
#This breaks the download into chunks
def download_parts(obsid, url, cleanup=False, silent=False):
    
    log.info("Attempting to curl {0} in parts".format(obsid))

    #Initial byte range
    byte_range = [0, 999999999]
    fname_i = 0
    fname = '{0}.tar.part{1}'.format(obsid, fname_i)

    cmd = ['curl', '--retry', '10', url, '--ftp-ssl', '-k', 
           '--create-dirs', '-o', fname, '--range',
           '{0}-{1}'.format(str(byte_range[0]), str(byte_range[1]))]

    if silent:
        cmd.append("-s")

    subprocess.call(cmd)

    if not os.path.isfile(fname):
        log.error("Download in parts with for {0} from {1} failed".format(
                obsid, url))
        return -1

    while os.path.getsize(fname) == 1e9:
        byte_range = [ b+1e9 for b in byte_range ]
        fname_i += 1
        fname = '{0}.tar.part{1}'.format(obsid, fname_i)
        cmd = ['curl', '--retry', '10', url, '--ftp-ssl', '-k', 
               '--create-dirs', '-o', fname, '--range',
               '{0}-{1}'.format(str(int(byte_range[0])), 
                                str(int(byte_range[1])))]

        if silent:
            cmd.append("-s")

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

#Some evts don't have 15 columns
#This checks which columns are missing and says if it's cause for concern
def check_invalid_columns(evt_path):
    required_columns = ['TIME', 'RAWX', 'RAWY', 'PHA', 'PHA_FAST', 
                        'DET_ID', 'DEADTIME', 'EVENT_FLAGS', 'TICK',
                        'MPU_A_TEMP', 'MPU_UNDER_COUNT', 'PI_FAST',
                        'PI', 'PI_RATIO', 'PULSE_PHASE']

    tab = Table.read(evt_path, hdu=1)

    existing_columns = tab.colnames

    missing_columns = [ c for c in required_columns 
                        if c not in existing_columns ]

    if len(missing_columns) != 0:
        if len(tab) == 0:
            log.info(
                    "Missing columns because no events left after filtering")
        else:
            log.warning("{0} is missing {1} column(s): {2}".format(
                    evt_path, len(missing_columns), missing_columns))

#paul's code from ni_data_download
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


def get_exposure(sources, user, passwd):
	df = print_nicer_segment(username=user, password=passwd)

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

def check_recent_obs(source, user, passwd, expected_dir='./', ncheck=5):

    #Get full target summary
    df_full = print_nicer_segment(username=user, password=passwd)

    #Grab the source of interest
    df_source = df_full[df_full['Target Name']==source]

    #reset df index
    df_source = df_source.reset_index(drop=True)

    dt_objects = [ datetime.datetime.strptime(d, '%Y-%m-%dT%H:%M:%S')
                   for d in df_source['Start TimeUTC'] ]
    now = datetime.datetime.now()
    days_ago = np.array([ (now-d).days for d in dt_objects ])

    #This is a nice way to pull the minimum ncheck values
    idx = np.argpartition(days_ago, ncheck)

    obs_list = []
    results_list = []
    for i in idx[:ncheck]:
        obsID = str(df_source['Observation ID'][i])
        date = df_source['Start TimeUTC'][i]

        obsdir = os.path.join(expected_dir, obsID)
        pipedir = os.path.join(expected_dir, obsID+'_pipe')
        evt_path = os.path.join(pipedir, 'cleanfilt.evt')

        #[observation dir exists, pipe dir exists, evt exists, evt length]
        results = ['X', 'X', 'X', 0]

        if os.path.isdir(obsdir):
            results[0] = u'\u2713'
            if os.path.isdir(pipedir):
                results[1] = u'\u2713'
                if os.path.isfile(evt_path):
                    results[2] = u'\u2713'
                    tab = Table.read(evt_path, hdu=1)
                    results[3] = len(tab)

        print("Observation {0} from {1}:\n\t\tObservation directory: {2}\n\t\tPipe directory: {3}\n\t\tEvent file exists: {4} w/ length {5}".format(
                obsID, date, *results))

        obs_list.append(obsID)
        results_list.append(results)

    df_out = pd.DataFrame(results_list, index=obs_list, 
                          columns=['obsdir', 'pipedir', 'evt', 'len(evt)'])
    return df_out

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


def run_photonphase(evt, orbfile, par, ephem='DE421'):

    if not os.path.isfile(evt):
        raise FileNotFoundError("event file not found")

    if not os.path.isfile(orbfile):
        raise FileNotFoundError("orbit file not found")

    if not os.path.isfile(par):
        raise FileNotFoundError("par file not found")

    cmd = ['photonphase', '--ephem', ephem, '--orb', orbfile,
           '--addphase', evt, par]

    subprocess.call(cmd)

    return 0


def split_photonphase(evt, orbfile, par, split_len=100000, use_mp=False):

    log.info("Running split photonphase routine")

    if not os.path.isfile(evt):
        raise FileNotFoundError("event file not found")

    if not os.path.isfile(orbfile):
        raise FileNotFoundError("orbit file not found")

    if not os.path.isfile(par):
        raise FileNotFoundError("par file not found")

    #Read in full table and get length
    full_table = Table.read(evt, hdu=1)
    total_length = len(full_table)

    #Split files and get names
    if total_length <= split_len:
        run_photonphase(evt, orbfile, par)
        return 0
    else:
        #make a copy of the evt
        log.info("Copying {0} to {1}".format(
                evt, evt.replace('.evt', '_nophase.evt')))
        cmd = ['cp', evt, evt.replace('.evt', '_nophase.evt')]
        subprocess.call(cmd)
        fnames = split_evt(evt, split_len)

    #Check that the split completed successfully
    if not check_split(evt, fnames):
        log.error("fselect split for {0}, exiting".format(evt))
        return -1

    #Option for multiprocessing
    # we dont want to do this in the pipe wrapper 
    if use_mp:
        #Set up multiprocessing pool
        log.info("Using multiprocessing for photonphase")
        pool = mp.Pool(processes=mp.cpu_count()+2)
        jobs = []

        #Add each split evt to pool
        for f in fnames:
            job = pool.apply_async(
                    run_photonphase,
                    (f, orbfile, par,))
            jobs.append(job)

        #run jobs
        for job in jobs:
            job.get()

    else: #no multiprocessing
        for f in fnames:
            run_photonphase(f, orbfile, par)

    #Check that all the files have phase columns
    for f in fnames:
        tab = Table.read(f, hdu=1)
        if 'PULSE_PHASE' not in list(tab.columns):
            log.error("no pulse phase for split evt {0}, exiting".format(f))
            return -1

    #Write file list for merge
    split_file_list = os.path.join(os.path.split(evt)[0], 'split_evt_list')
    with open(split_file_list, 'w') as f:
        for x in fnames:
            f.write(x+'\n')
    
    log.info("Merging with niextract-events")
    #Call niextract to merge
    cmd = ['niextract-events', 'filename=@{0}'.format(split_file_list),
           'eventsout={0}'.format(evt.replace('.evt', '_phase.evt')),
           'clobber=YES']

    subprocess.call(cmd)

    log.info("Copying {0} to {1}".format(evt.replace('.evt', '_phase.evt'),
                                         evt))
    cmd = ['cp', evt.replace('.evt', '_phase.evt'), evt]
    subprocess.call(cmd)

def row_expression(row_range):
    return "#row >= {0} && #row <= {1}".format(*row_range)

def split_evt(evt, split_len, clobber=True, cleanup=False):
    
    log.info("splitting {} into parts".format(evt))

    #Check path
    if not os.path.isfile(evt):
        raise FileNotFoundError("Event file not found")

    #Make dir for split evt files
    pipe_dir_path = os.path.split(evt)[0]
    split_evts_dir_path = os.path.join(pipe_dir_path, 'split_evts')

    if os.path.isdir(split_evts_dir_path):
        log.warning("split_evts dir already exists, removing")
        shutil.rmtree(split_evts_dir_path)

    log.info("Creating dir for split evts")
    os.mkdir(split_evts_dir_path)
            
    row_range = [1, split_len]

    #Get split evt file name
    evt_number = 0
    evt_base = os.path.split(evt)[1]
    evt_split = os.path.join(
            split_evts_dir_path, 
            evt_base.replace('.evt', '_split{0}.evt'.format(evt_number)))

    cmd = ['fselect', evt, evt_split, row_expression(row_range)]
    subprocess.call(cmd)

    tab = Table.read(evt_split, hdu=1)

    #Save all the filenames
    fnames = [evt_split]

    while len(tab) == split_len:
        row_range = [v+split_len for v in row_range]
        evt_number += 1
        evt_split = os.path.join(
                split_evts_dir_path, 
                evt_base.replace('.evt', '_split{0}.evt'.format(evt_number)))

        fnames.append(evt_split)
        cmd = ['fselect', evt, evt_split, row_expression(row_range)]
        subprocess.call(cmd)
        tab = Table.read(evt_split, hdu=1)

    log.info("{0} split into {1} parts".format(evt, len(fnames)))

    return fnames

#Function to check if split completed
def check_split(evt, fnames):
    
    if not os.path.isfile(evt):
        raise FileNotFoundError("Event file not found")

    if not niutils.check_iter(fnames):
        raise TypeError("fnames must be iterable")

    for f in fnames:
        if not os.path.isfile(f):
            raise FileNotFoundError("{0} not found".format(f))

    full_tab = Table.read(evt, hdu=1)
    total_length = len(full_tab)

    summed = 0
    for split_evt in fnames:
        tab = Table.read(split_evt, hdu=1)
        summed += len(tab)

    return summed == total_length
