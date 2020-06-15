#!/usr/bin/env python

from __future__ import print_function, division
from astropy import log
import argparse
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


def wrapper(br_earth=0, elv=15, sun_angle=20, clobber=False):

    pool = mp.Pool(processes=mp.cpu_count()+2)
    jobs = []

    for f in os.listdir("./"):
        if not f.endswith("_pipe"):
            if os.path.isdir(f) and (f != 'tmp'):
                if not os.path.isfile(f+"_pipe/cleanfilt.evt"):
                    job = pool.apply_async(bkg_pipe, (f,),
                                           dict(br_earth=br_earth,
                                                elv=elv,
                                                sun_angle=sun_angle,
                                                clobber=clobber))
                    jobs.append(job)

    for job in jobs:
        job.get()


def bkg_pipe(obsID, br_earth=0, elv=15, sun_angle=20, clobber=False):

    if pipeline_utils.check_nicerl2(obsID) and (not clobber):
        pass
    elif pipeline_utils.check_nicerl2(obsID) and clobber:
        pipeline_utils.run_nicerl2(obsID, br_filter=False, clobber=True)
    else:
        pipeline_utils.run_nicerl2(obsID, br_filter=False)

    pipeline_utils.run_add_kp(obsID)

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

    cor_sax_string = ['(COR_SAX.gt.(1.914*KP**0.684+0.25))']
    kp_string = ['KP.lt.5']

    sunshine_string = ['(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)'.format(sun_angle)]
    extra_expr = cor_sax_string+kp_string+sunshine_string
    extra_expr = "("+" && ".join("%s" %expr for expr in extra_expr) + ")"

    maxunder=200.0

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

def update(heasarc_user, heasarc_pwd, decryptkey):

    region_nums = [1,2,3,4,5,6,8]
    region_names = [ "BKGD_RXTE_{0}".format(i) for i in region_nums ]

    #Check that the region directories exist
    for d in region_names:
        if not os.path.isdir(d):
            raise FileNotFoundError("Missing source directory {}".format(d))

    #Go through each, download and pipe
    for source in region_names:
        os.chdir(source)
        pipeline_utils.run_datadownload(source, heasarc_user, 
                                        heasarc_pwd, './',
                                        decryptkey, clobber=False)

        wrapper(br_earth=0, elv=15, sun_angle=20, clobber=False)
        os.chdir('../')

    #merge evt and mkf
    niextract_all_events()

    merged_evt = 'bkgd_combined.evt'
    merged_mkf = 'bkgd_merged.mkf'
    message='Created with background_pipe.update'

    pipeline_utils.product_backup(merged_evt, backupdir='evt_backups',
                                  message=message)
    pipeline_utils.product_backup(merged_mkf, backupdir='mkf_backups',
                                  message=message)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
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
    parser.add_argument("--update", help="Run pipe on all regions",
                        default=False, action='store_true')
    args = parser.parse_args()

    if args.k_file is not None:
        with open(args.k_file, 'r') as f:
            args.key = f.readlines()[0].strip("\n")

    if args.update:
        update(args.user, args.passwd, args.key)
