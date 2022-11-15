#!/usr/bin/env python

from __future__ import print_function, division
from astropy.table import Table
from astropy import log
import datetime
import os, sys
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import shutil
import subprocess
import time
import pipeline.pipeline_utils2 as pipeline_utils2
from os import path
from glob import glob
from subprocess import check_call
from astropy.io import fits
import tempfile
from nicer.values import *
from nicer.plotutils import find_hot_detectors
import fitsio


# CMD should be a list of strings since it is not processed by a shell
def runcmd(cmd):
    log.info('CMD: '+" ".join(cmd))
    log.info(cmd)
    check_call(cmd,env=os.environ)

obsID=str(input("ObsID: "))
par = str(input("parfile: "))

emin=0.25
emax=12
trumpet=True
keith=True
clobber=False
crab=True
min_sun=60
mask = [14,34,54]
cormin = 1.5
angdist=0.005

log.info('Make sure there is actually data')
if not os.path.isdir(obsID+'/xti/event_cl'):
    log.error("no event cl dir for {}".format(obsID))

pipeline_utils2.run_nicerl2(obsID, trumpet=trumpet, clobber=clobber)
pipeline_utils2.run_add_kp(obsID)

pipedir = "{0}_{1}".format(obsID,'pipe')
if not os.path.exists(pipedir):
    os.makedirs(pipedir)
"""
cmd = ['nicerql.py', '--save', '--filtall', '--lcbinsize', '4.0', '--allspec', '--alllc', '--lclog', '--useftools', '--extraphkshootrate', '--emin', '0.25', '--emax', '12.0', '--bkg', '--map', '--obsdir', obsID, '--basename', '%s/%s_prefilt'%(pipedir, obsID), '--mask', '-1', '--keith', '--par', par]
runcmd(cmd)

"""
cmd = ['nimaketime', 'infile=%s/auxil/ni%s.mkf'%(obsID, obsID), 'outfile=%s/tot.gti'%pipedir, 'nicersaafilt=YES', 'saafilt=NO', 'trackfilt=YES', 'ang_dist = 0.005', 'elv=20.0', 'br_earth=30.0','cor_range=1.5-', 'min_fpm=38', 'underonly_range=0-200.0', 'ingtis=NONE', 'clobber=yes', 'expr=(ST_VALID.eq.1 && FPM_OVERONLY_COUNT<1 && FPM_OVERONLY_COUNT<(1.52*COR_SAX**(-0.633)) && (COR_SAX.gt.(1.914*KP**0.684+0.25)).and.KP.lt.5 && (SUN_ANGLE.gt.60.or.SUNSHINE.eq.0))', 'outexprfile=%s/psrpipe_expr.txt'%pipedir]
runcmd(cmd)

evfiles = glob(path.join(obsID,'xti/event_cl/ni*mpu7_cl.evt'))
mkfile = glob(path.join(obsID,'auxil/ni*.mkf'))[0]
shutil.copy(mkfile,pipedir)
hkfiles = glob(path.join(obsID,'xti/hk/ni*.hk'))
hkfiles.sort()

evlistname=path.join(pipedir,'evfiles.txt')
fout = open(evlistname,'w')
for en in evfiles:
    print('{0}'.format(en),file=fout)
fout.close()
evfilt_expr = 'PI={0}:{1},EVENT_FLAGS=bx1x000'.format(int(emin*KEV_TO_PI), int(emax*KEV_TO_PI))

cmd = ['niextract-events', "filename=@{0}[{1}]".format(evlistname,evfilt_expr), 'eventsout=%s/intermediate.evt'%pipedir, 'timefile=%s/tot.gti'%pipedir, 'gti=GTI', 'clobber=yes']
runcmd(cmd)

"""
tb = Table.read("%s/intermediate.evt"%pipedir, hdu=1)
bad_dets = find_hot_detectors(tb)
if bad_dets is not None:
    log.info('Found hot detectors {0}!!'.format(bad_dets))

evfilt_expr = '(EVENT_FLAGS==bx1x000)'
if bad_dets is not None:
    fout = open(path.join(pipedir,"bad_dets.txt"),"w")
    print("{0}".format(bad_dets),file=fout)
    fout.close()
    for detid in bad_dets:
        evfilt_expr += ".and.(DET_ID!={0})".format(detid)
"""

cmd = ['ftcopy', '%s/intermediate.evt[(EVENT_FLAGS==bx1x000).and.(DET_ID!=14).and.(DET_ID!=34).and.(DET_ID!=54)]'%pipedir, '%s/cleanfilt.evt'%pipedir, 'clobber=yes', 'history=yes']
runcmd(cmd)

cmd = ['fltime', 'infile=%s/auxil/ni%s.mkf[1]'%(obsID, obsID), 'gtifile=%s/tot.gti'%pipedir, 'outfile=%s/cleanfilt.mkf'%pipedir, 'clobber=yes']
runcmd(cmd)

orbfile = glob(path.join(obsID,'auxil/ni*.orb'))[0]
log.info('Orbit File: {0}'.format(orbfile))
shutil.copy(orbfile,pipedir)

if os.path.isfile("%s/cleanfilt2.evt"%pipedir):
    os.remove("%s/cleanfilt2.evt"%pipedir)

cmd = ["barycorr", "%s/cleanfilt.evt"%pipedir, "%s/cleanfilt2.evt"%pipedir, "%s/ni%s.orb"%(pipedir, obsID)]
runcmd(cmd)

h = fitsio.read_header('%s/intermediate.evt'%pipedir, ext=1) 
n_events = h["NAXIS2"]
log.info('There are %s events for this ObsID'%n_events)
if n_events == 0:
    raise ValueError("Zero valid events for this ObsID")

if n_events < 1000000:
    cmd = ["photonphase3.py", '--polycos', '--counts', '100000', '--overwrite', '%s/cleanfilt2.evt'%pipedir, par]
if n_events >= 1000000:
    cmd = ["photonphase3.py", '--polycos', '--overwrite', '%s/cleanfilt2.evt'%pipedir, par]
runcmd(cmd)

cmd = ['extractor', '%s/cleanfilt2.evt'%pipedir, 'eventsout=none', 'imgfile=none', 'phafile=%s/cleanfilt.pha'%pipedir, 'fitsbinlc=%s/cleanfilt.lc'%pipedir, 'binlc=4.0', 'regionfile=none', 'timefile=none', 'xcolf=RAWX', 'ycolf=RAWY', 'tcol=TIME', 'ecol=PI', 'xcolh=RAWX', 'ycolh=RAWY', 'gti=GTI']
runcmd(cmd)

"""
cmd = ["nicerql.py", "--save",
           "--orb", "%s/ni%s.orb"%(pipedir, obsID),
           "--map",
           "--sci", "--eng", '%s/cleanfilt2.evt'%pipedir,"--allspec","--alllc",
           "--lcbinsize", "4.0",
           "--filterbinsize", "16.0",
           "--mkf", "%s/cleanfilt.mkf"%(pipedir), "--bkg",
           "--basename", '%s/%s_cleanfilt'%(pipedir, obsID), 
           "--par", par, "--keith"]
runcmd(cmd)
"""
print('DONE')


