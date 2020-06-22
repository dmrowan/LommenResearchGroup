#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
import argparse
import datetime
from glob import glob
import subprocess 
import multiprocessing as mp
import os
import shutil

#LommenResearchGroup imports
import pipeline_utils

#Dom Rowan 2020

desc="""
Functions for processing and extracting background data
"""


#Multiprocessing wrapper for background pipe
def wrapper(br_earth=0, elv=15, sun_angle=20, clobber=False):

    #Initialize pool
    pool = mp.Pool(processes=mp.cpu_count()+2)
    jobs = []

    #Iterate through directory to search for obsIDs
    for f in os.listdir("./"):
        if os.path.isdir(f) and (f.isnumeric()):
            #Apply background pipe
            job = pool.apply_async(bkg_pipe, (f,),
                                   dict(br_earth=br_earth,
                                        elv=elv,
                                        sun_angle=sun_angle,
                                        clobber=clobber))
            jobs.append(job)

    #run each job in the pool
    for job in jobs:
        job.get()


#Background pipe for single observation
def bkg_pipe(obsID, br_earth=0, elv=15, sun_angle=20, clobber=False):

    if not os.path.isdir(obsID):
        raise FileNotFoundError("observation not found")

    #If the cleaned event file already exists and clobber=False, return
    if os.path.isfile(obsID+'_pipe/cleanfilt.evt') and (not clobber):
        return 0

    #if nicerl2 has been run go to next step
    if pipeline_utils.check_nicerl2(obsID) and (not clobber):
        pass
    #otherwise run nicerl2
    else:
        pipeline_utils.run_nicerl2(obsID, br_filter=False, clobber=clobber)

    #Add geomagnetic index information
    pipeline_utils.run_add_kp(obsID)

    #Create pipe directory
    pipedir = "{0}_pipe".format(obsID)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)
	
    #copy relevant files to the pipe directory
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

    #space weather filtering expressions
    cor_sax_string = ['(COR_SAX.gt.(1.914*KP**0.684+0.25))']
    kp_string = ['KP.lt.5']

    #sun angle criteria
    sunshine_string = ['(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)'.format(sun_angle)]
    extra_expr = cor_sax_string+kp_string+sunshine_string
    extra_expr = "("+" && ".join("%s" %expr for expr in extra_expr) + ")"

    maxunder=200.0

    #nimaketime command
    cmd = ['nimaketime', 'infile={0}'.format(mkf),
           'outfile={0}'.format(gti), 'nicersaafilt=YES',
           'saafilt=NO', 'trackfilt=YES', 'ang_dist=0.015', 
           'elv={0}'.format(elv),
           'br_earth={0}'.format(br_earth), 'min_fpm=7',
           'underonly_range=0-{0}'.format(maxunder),
           'expr={0}'.format(extra_expr), 
           'outexprfile={0}'.format(
                    os.path.join(pipedir, "bkgpipe_expr.txt")), 
           'clobber=YES']
	
    subprocess.call(cmd)

    #Gather event files
    evfiles = glob(os.path.join(obsID,'xti/event_cl/ni*mpu7_cl.evt'))
    evlistname=os.path.join(pipedir, 'evfiles.txt')
    with open(evlistname, 'w') as f:
        for ev in evfiles:
            f.write(ev+'\n')
            evfiles.sort()
            log.info("Cleaned Event Files: {0}".format(
                     "\n		\n".join(evfiles)))

	
    #Extract events using the gti
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

#merge event file for all observations
def niextract_all_events(verbose=False, max_date=None, 
                         output=None):

    log.info("Merging event and filter files")

    #check existence of each target directory
	region_nums = [1,2,3,4,5,6,8]
	region_names = [ "BKGD_RXTE_{0}".format(i) for i in region_nums ]
	for d in region_names:
		assert(os.path.isdir(d))

	filenames = []
	mkfnames = []
    #Go through each region
	for d in region_names:
		for f in os.listdir(d):
			if os.path.isfile(d+"/"+f+"/cleanfilt.evt"):
                #Date filtering
                if max_date is not None:
                    #Get table date from meta
                    tab = Table.read(f+"/cleanfilt.evt", hdu=1)
                    date = datetime.datetime.strptime(
                            tab.meta['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
                    if date > max_date:
                        print("Excluding {0} from {1}".format(f, str(date)))
                        continue
                    
				filenames.append(d+"/"+f+"/cleanfilt.evt")
				mkfnames.append(d+"/"+f+"/mkf_filtered.mkf")
	
    #Write files with evt and mkf lists
	with open("bkgd_all_evt_list", 'w') as f:
		for fn in filenames:
			f.write(fn+"\n")
	with open("bkgd_all_mkf_list", 'w') as f1:
		for fn in mkfnames:
			f1.write(fn+"\n")

    #we handle the output name like this to make it easier for
    # command line calling
    if output is None:
        output='bkgd_combined.evt'

    #run niextract to merge events
	cmd = ['niextract-events', 'filename=@bkgd_all_evt_list',
		   'eventsout={}'.format(output), 'clobber=yes']
	if verbose: cmd.append('chatter=5')
	subprocess.call(cmd)

    #run ftmerge to merge filter files
    mkf_output = output.replace('.evt', '.mkf')
	cmd = ['ftmerge', 'infile=@bkgd_all_mkf_list',
		   'outfile={}'.format(mkf_output), 'clobber=YES']
	if verbose: cmd.append('chatter=5')
	subprocess.call(cmd)

#update dataset for all background regions
def update(heasarc_user, heasarc_pwd, decryptkey, br_earth=0,
           elv=15, sun_angle=20, silent_curl=False):

    region_names = [ "BKGD_RXTE_{0}".format(i) for i in [1,2,3,4,5,6,8]

    #Check that the region directories exist
    for d in region_names:
        if not os.path.isdir(d):
            raise FileNotFoundError("Missing source directory {}".format(d))

    #Go through each, download and pipe
    for source in region_names:
        os.chdir(source) #Move into region directory
        pipeline_utils.run_datadownload(source, heasarc_user, 
                                        heasarc_pwd, './',
                                        decryptkey, clobber=False,
                                        silent_curl=silent_curl)

        wrapper(br_earth=br_earth, elv=elv, sun_angle=sun_angle, 
                clobber=False)
        os.chdir('../')

    #merge evt and mkf
    niextract_all_events()

    #Create backup files for evt and mkf
    merged_evt = 'bkgd_combined.evt'
    merged_mkf = 'bkgd_merged.mkf'
    message='Created with background_pipe.update'

    pipeline_utils.product_backup(merged_evt, backupdir='evt_backups',
                                  message=message)
    pipeline_utils.product_backup(merged_mkf, backupdir='mkf_backups',
                                  message=message)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)

    #Arguments for accessing data from NASA site
    parser.add_argument("--user", help="heasarc username",
                        default='nicer_team', type=str)
    parser.add_argument("--passwd", help="heasarc password",
                        default='sextant', type=str)
    parser.add_argument("--outdir", help="outdir for download", 
                        default='./', type=str)
    parser.add_argument("-k", help="decryption key", dest='key',
                        action=pipeline_utils.PasswordPromptAction,
                        type=str, required=False)
    parser.add_argument("--k_file", help="decryption key from file",
                        type=str, default=None)

    parser.add_argument("--clobber", help="overwrite existing data", 
                        default=False, action='store_true')

    #Arguments for filtering parameters
    parser.add_argument("--br_earth", help="minimum bright earth angle",
                        default=0, type=int)
    parser.add_argument("--elv", help="minimum elevation angle",
                        default=15, type=int)
    parser.add_argument("--sun_angle", help="minimum sun angle",
                        default=20, type=int)

    #Arguments for different function calls
    parser.add_argument("--update", help="Run pipe on all regions",
                        default=False, action='store_true')
    parser.add_argument("--obsID", help="Argument for single pipeline call",
                        default=None, type=str)
    parser.add_argument("--download", help="Run data download for source",
                        default=None, type=str)
    parser.add_argument("--backup", help="create backup for file",
                        default=None, type=str)
    parser.add_argument("--message", help="Log message for evt/mkf backup",
                        default=None, type=str)
    #Arguments for merge call
    parser.add_argument("--merge", help="Combine evt/mkf",
                        default=False, action='store_true')
    parser.add_argument("--output", help="output name for merge",
                        default=None, type=str)
    parser.add_argument("--max_date", help="max date to merge",
                        default=None, type=str)
    args = parser.parse_args()

    #Option to silence curl progess bar
    parser.add_argument("--silent_curl", help="no curl progress bar",
                        default=False, action='store_true')

    if args.k_file is not None:
        with open(args.k_file, 'r') as f:
            args.key = f.readlines()[0].strip("\n")

    if args.update:
        update(args.user, args.passwd, args.key, 
               br_earth=args.br_earth, elv=args.elv,
               sun_angle=args.sun_angle, 
               silent_curl=args.silent_curl)

    elif args.obsID is not None:
        bkg_pipe(obsID, br_earth=args.br_earth, elv=args.elv, 
                 sun_angle=args.sun_angle, clobber=args.clobber):

    elif args.download:
        sources = ['BKGD_RXTE_{}'.format(i) for i in [1,2,3,4,5,6,8]]
        for s in sources:
            pipeline_utils.run_datadownload(s, args.user, args.passwd,
                                            s, args.key, 
                                            clobber=args.clobber,
                                            silent_curl=args.silent_curl)

    #Option to just merge files
    #optional arguments to merge up to specific date and give output name
    elif args.merge is not None:
        if args.max_date is None:
            run_niextract(output=args.output)
        else:
            dt = datetime.datetime.strptime(args.max_date, '%Y-%m-%d')
            run_niextract(output=args.output, max_date=dt)

    elif args.backup is not None:
        if args.backup.endswith('.evt'):
            evt = args.backup
            mkf = args.backup.replace('.evt', '.mkf')
        elif args.backup.endswith('.mkf'):
            evt = args.backup.replace('.mkf', '.evt')
            mkf = args.backup
        else:
            evt = args.backup+'.evt'
            mkf = args.backup+'.mkf')

        pipeline_utils.product_backup(evt, backupdir='evt_backups',
                                      message=args.message)
        pipeline_utils.product_backup(mkf, backupdir='mkf_backups',
                                      message=args.message)

    else:
        log.warning("No method selected! Try 'background_pipe.py --help'")
