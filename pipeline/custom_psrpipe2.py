#!/usr/bin/env python
from __future__ import print_function, division
import os, sys
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
from astropy.table import Table
from astropy.io import fits
import tempfile

from nicer.values import *
from nicer.plotutils import find_hot_detectors

#LommenResearchGroup imports
from pipeline import pipeline_utils2

desc = """
Pipeline process NICER data.

Output will be written in current working directory in directories that end in '_pipe'.

Several diagnostic plots are produced, and the following data processing steps are run:
* Select good times according to the following:
  * (ANG_DIST.lt.0.015).and.(ELV>30.0)
  * (MODE.eq.1).and.(SUBMODE_AZ.eq.2).and.(SUBMODE_EL.eq.2)
  * SAT_LAT, SAT_LONG not in SAA or polar horn regions specified by region files
  * If --dark is specified then also filter on SUNSHINE.eq.0
Optionally, you can filter on overshoot rate or rate of bad ratio events

* Select events according to:
  * EVENT_FLAGS=bx1x000 (SLOW-only or SLOW+FAST events)
  * PI in the specified range (default is 0.25-12.0 keV)
  * Remove events from any DET_ID specified by the --mask parameter

 The final output is 'cleanfilt.evt' and extracted PHA and lightcurve FITS
 files produced by extrator. The event file will have a PULSE_PHASE column
 computed using PINT if the --par file is specified on the command line.
"""
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=12.0)", type=float, default=12.0)
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--filtpolar",help="Turn on  filtering polar horn regions from data",default=False,action='store_true')
parser.add_argument("--cormin",help="Set minimum cutoff rigidity (COR_SAX) for nimaketime filtering (default=no COR filtering, typical value = 4)",default=None)
parser.add_argument("--kpmax",help="Set maximum KP value for nimaketime filtering (default=no KP filtering, typical value = 5)",default=None)
parser.add_argument("--minfpm",help="Set minimum of FPMs active for nimaketime filtering (default=38)",default=38)
parser.add_argument("--maxovershoot",help="Select data where overshoot rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--badcut",help="Select data where bad ratio event rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--angdist",help="Set threshold for ANG_DIST in call to nimaketime (degrees, default=0.015)", type=float, default=0.005)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
parser.add_argument("--shrinkelvcut", help="Shrink ELV cut to 20 deg and BR_EARTH cut to 30.0 deg to get more data (this is now ignored since it is the default)", action='store_true')
parser.add_argument("--dark", help="Apply SUNSHINE=0 filter to get only data in Earth shadow", action='store_true')
parser.add_argument("--nounderfilt", help="Don't filter good times based on UNDERONLY rate", action='store_true')
parser.add_argument("--minsun",help="Set minimum sun angle (SUN_ANGLE) for nimaketime filtering (default=no SUN_ANGLE filtering, typical values = 60, 70, 80, 90 deg). Note: Allows dark time at any Sun angle!",default=None)
parser.add_argument("--day", help="Apply SUNSHINE=1 filter to get only data in ISS-day", action='store_true')
parser.add_argument("--par", help="Par file to use for phases")
parser.add_argument("--ephem", help="Ephem to use with photonphase", default="DE421")
parser.add_argument("--outdir", help="Add name to output directories (by default: directories end in '_pipe')", default='pipe')
parser.add_argument("--merge", help="Merge all ObsIDs provided into single event list, lightcurve and spectrum (outputdir called 'merged')", action='store_true')
parser.add_argument("--crcut", help="perform count rate cut on merged event file (only if --merge)", action='store_true')
parser.add_argument("--lcbinsize", help="Lightcurve bin size (sec, default=4.0)", type=float, default=4.0)
parser.add_argument("--filterbinsize", help="Bin size for Count rate and Overshoot rate filtering (sec, default=16.0)", type=float, default=16.0)
parser.add_argument("--keith", help="Standard filters used by Keith Gendreau for Space-Weather backgrounds (Masks detectors 14, 34, 54; cormin 1.5; custom cut on overshoot rate + COR_SAX)", action='store_true')
# parser.add_argument("--crabnorm", help="normalize the spectrum with the crab (only if --merge)", action='store_true')

#DOM ADDED THIS
parser.add_argument("--crab", help='Use multiprocessing for photonphase step',
                    default=False, action='store_true')
args = parser.parse_args()

os.environ['HEADASNOQUERY'] = ' '
os.environ['HEADASPROMPT'] = '/dev/null'

# Checking the presence of HEASOFT
try:
    check_call('nicerversion',env=os.environ)
except:
    print("You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!", file=sys.stderr)
    exit()

# Create a temporary dir for pfiles and set the environment
tempdir = tempfile.mkdtemp()
os.mkdir(tempdir+'/pfiles')
headas = os.environ['HEADAS']
os.environ["PFILES"] = tempdir+'/pfiles;'+headas+'/syspfiles'
#print(os.environ["PFILES"])
# Called before each exit() and at end of program to clean up:
#shutil.rmtree(tempdir)

# For some reason if the script is called via #!/usr/bin/env python
# it does not inherit LD_LIBRARY_PATH so ftools don't run.
#print(os.environ['LD_LIBRARY_PATH'])
#print(os.environ['HEADAS'])
#os.environ['LD_LIBRARY_PATH'] = path.join(os.environ['HEADAS'],'lib')
#print(os.environ['LD_LIBRARY_PATH'])

all_evfiles = []
all_orbfiles = []

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    log.info(cmd)
    check_call(cmd,env=os.environ)

# Check if outdir contains 'None', 'NONE', or 'none' (causes bug in ni-extractevents)
if args.outdir:
    names = ['none', 'None', 'NONE']
    if any(st in args.outdir for st in names):
        log.error("Due to a current bug in ni-extractevents, --outdir cannot contain 'none', 'None', or 'NONE'.  Exiting...")
        shutil.rmtree(tempdir)
        exit()

# Check if ObsIDs are listed or to be read from file
if len(args.indirs)==1:
    if args.indirs[0].startswith('@'):
        inputdirs = args.indirs[0].split('@')[1]
        log.info('Reading input ObsID list: {}'.format(inputdirs))
        all_obsids = np.loadtxt(inputdirs,dtype=str)
    else:
        all_obsids = args.indirs
else:
    all_obsids = args.indirs

# Start processing all ObsIDs
for obsdir in all_obsids:

    # Set up a basename and make a work directory
    if args.obsid is not None:
        basename = args.obsid
    else:
        # Trim trailing / if needed
        if obsdir.endswith("/"):
            obsdir = obsdir[:-1]
        basename = path.basename(obsdir)

    # Make directory for working files and output
    pipedir = "{0}_{1}".format(basename,args.outdir)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info('Making initial QL plots')

    # Get event filenames (could be just one)
    evfiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_cl.evt'))
    evfiles.sort()
    log.info('Cleaned Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    # Get filter file
    mkfile = glob(path.join(obsdir,'auxil/ni*.mkf'))[0]
    log.info('MKF File: {0}'.format(mkfile))
    # Check if MKF file contains the new columns (try opening one of the new columns)
    try:
        mkftest = fits.open(mkfile)[1].data['FPM_OVERONLY_COUNT']
    except:
        log.error("New *mkf files needed in {}/auxil/. Please use niprefilter2.".format(obsdir))
        exit()
    try:
        mkftest = fits.open(mkfile)[1].data['KP']
        has_KP = True
    except:
        log.warning('No KP column in MKF file. Will not use any KP based filters!')
        has_KP = False

    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(mkfile,pipedir)

    if args.keith:
        args.mask   = [14,34,54]
        args.cormin = 1.5

    if args.dark and args.day:
        log.warning("Both --dark and --day are requested")
        args.dark = False
        args.day = False

    # # Get ufa file (unfiltered events)
    # ufafiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_ufa.evt*'))
    # ufafiles.sort()
    # log.info('Unfiltered Event Files: {0}'.format("\n" + "    \n".join(ufafiles)))

    # ufaevents = 0
    # for ufafile in ufafiles:
    #     hdulist = fits.open(ufafile, memmap=True)
    #     nevents = hdulist[1].header['NAXIS2']
    #     hdulist.close()
    #     ufaevents = ufaevents + nevents

    # if ufaevents < 10000000:
    #     cmd = ["nicerql.py", "--save", "--filtall", "--lcbinsize", "{}".format(args.lcbinsize), ##  "--allspec","--alllc",
    #            "--lclog", "--useftools", "--extraphkshootrate",
    #            "--eventshootrate",
    #            "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
    #            "--sci", "--eng", "--bkg", "--map", "--obsdir", obsdir,
    #            "--basename", path.join(pipedir,basename)+'_prefilt']
    # else:
    #     cmd = ["nicerql.py", "--save", "--filtall", "--lcbinsize", "{}".format(args.lcbinsize), ## "--allspec","--alllc",
    #            "--lclog", "--useftools", "--extraphkshootrate",
    #            "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
    #            "--sci", "--eng", "--bkg", "--map", "--obsdir", obsdir,
    #            "--basename", path.join(pipedir,basename)+'_prefilt']

    """
    cmd = ["nicerql.py", "--save", "--filtall",
           "--lcbinsize", "{}".format(args.lcbinsize), "--lclog", "--useftools",
           "--filterbinsize", "{}".format(args.filterbinsize),
           "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
           "--sci", "--eng", "--map", "--obsdir", obsdir,
           "--basename", path.join(pipedir,basename)+'_prefilt']
    #Dom removed bkg argument from nicerql

    if args.mask is not None:
        cmd.append("--mask")
        for detid in args.mask:
            cmd.append("{0}".format(detid))
    if args.keith:
        cmd.append("--keith")

#    if args.par is not None:
#        cmd.append("--par")
#        cmd.append("{0}".format(args.par))
    # if (args.maxovershoot>0) or (args.badcut>0):
    if (args.badcut>0):
        cmd.append("--writebkf")

    runcmd(cmd)
    """
    # Get orbit file
    orbfile = glob(path.join(obsdir,'auxil/ni*.orb'))[0]
    log.info('Orbit File: {0}'.format(orbfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(orbfile,pipedir)

    #  Get ATT hk files
    attfile = glob(path.join(obsdir,'auxil/ni*.att'))[0]
    log.info('ATT HK File: {0}'.format(attfile))

    #  Get BKF file for filtering based on background indicators
    # if (args.maxovershoot>0) or (args.badcut>0):
    if (args.badcut>0):
        bkffile = path.join(pipedir,basename)+'_prefilt.bkf'
        log.info('BKF File: {0}'.format(bkffile))
    else:
        bkffile=None

    #  Get MPU hk files
    hkfiles = glob(path.join(obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()
    log.info('MPU HK Files: {0}'.format("\n" + "    \n".join(hkfiles)))

    # Create any additional GTIs beyond what nimaketime does...
    extragtis="NONE"
    if args.filtpolar:
        saafile = path.join(datadir,'saa.reg')
        mkf_expr = 'regfilter("{0}",SAT_LON,SAT_LAT)'.format(saafile)
        phfile = path.join(datadir,'polarhorns.reg')
        mkf_expr += '.and.regfilter("{0}",SAT_LON,SAT_LAT)'.format(phfile)
        gtiname2 = path.join(pipedir,'extra.gti')
        cmd = ["pset", "maketime", "expr={0}".format(mkf_expr)]
        runcmd(cmd)
        cmd = ["maketime", "infile={0}".format(mkfile), "outfile={0}".format(gtiname2),
            "compact=no", "time=TIME",  "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname2,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue
        extragtis = gtiname2

    gtiname3 = None
    # # Create GTI from overshoot file using overshoot rate
    # if args.maxovershoot > 0:
    #     gtiname3 = path.join(pipedir,'bkf.gti')
    #     bkf_expr = 'EV_OVERSHOOT.lt.{0}'.format(args.maxovershoot)
    #     cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
    #         "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
    #     runcmd(cmd)
    #     if len(Table.read(gtiname3,hdu=1))==0:
    #         log.error('No good time left after filtering!')
    #         continue


    # Create GTI from overshoot file using overshoot rate
    if args.maxovershoot > 0:
        gticolumns = path.join(datadir,'gti_columns.txt')
        gtiheader = path.join(datadir,'gti_header.txt')

        if not os.path.isfile(gtiheader) or not os.path.isfile(gticolumns):
            log.error('The files gti_header.txt or gti_columns.txt are missing. Check the {} directory'.format(os.path.abspath(datadir)))
            log.error('No filtering on MaxOverShoot will be performed. Continuing...')
            gtiname3 = None
            continue
        else:
            ovbinfile = path.join(pipedir,basename)+'_prefilt_ovbin.fits'
            log.info("Filtering overshoots in {} with {} max overshoots/bin".format(ovbinfile,args.maxovershoot))
            gtiname3 = path.join(pipedir,'max_overshoots.gti')
            gtidata  = path.join(pipedir,'max_overshoots.txt')
            filtovbinfile = path.join(pipedir,basename)+'_cleanfilt_maxovbin.fits'

            ## STEP 2 -- making cut with ftcopy, and conditions, e.g. < 20 and =/= 0.0
            cmd = ['ftcopy', '{}[1][FPM_OVERONLY_COUNT<{} && FPM_OVERONLY_COUNT!=0]'.format(ovbinfile,args.maxovershoot/52.0),
                   '{}'.format(filtovbinfile),'clobber=yes']
            runcmd(cmd)

            ## STEP 3 - calculate start and end times of remaining bins
            cmd = ['ftcalc','{}'.format(filtovbinfile),'{}'.format(filtovbinfile),
                   'TSTART', '\"TIME-(0.5*{0})+#TIMEZERO\"'.format(args.filterbinsize),'clobber=yes']
            log.info('CMD: '+" ".join(cmd))
            os.system(" ".join(cmd))
            cmd = ['ftcalc','{}'.format(filtovbinfile),'{}'.format(filtovbinfile),
                   'TEND', '\"TIME+(0.5*{0})+#TIMEZERO\"'.format(args.filterbinsize),'clobber=yes']
            log.info('CMD: '+" ".join(cmd))
            os.system(" ".join(cmd))

            ## STEP 4 - dumping the TSTART and TEND into text file
            cmd = ['ftlist', '{}[1]'.format(filtovbinfile), 'columns=TSTART,TEND',
                   'rownum=no','colheader=no', 'opt=t', '>', '{}'.format(gtidata)]
            log.info('CMD: '+" ".join(cmd))
            os.system(" ".join(cmd))
            # runcmd(cmd)

            ##  STEP 6 - Making the GTI file from the text file
            cmd = ['ftcreate', '{}'.format(gticolumns),'{}'.format(gtidata), '{}'.format(gtiname3),
                   'headfile={}'.format(gtiheader), 'extname="GTI"', 'clobber=yes']
            runcmd(cmd)

            if len(Table.read(gtiname3,hdu=1))==0:
                log.error('No good time left after filtering!')
                continue
    else:
        gtiname3 = None


    # Create GTI from overshoot file using bad event lightcurve
    if args.badcut > 0:
        gtiname3 = path.join(pipedir,'bkf.gti')
        bkf_expr = 'BAD_RATIO.lt.{0}'.format(args.badcut)
        cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
            "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname3,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue

    # If either of the bkf filters were used, include that GTI
    # in the extragtis passed to nimaketime
    if gtiname3 is not None:
        if extragtis == "NONE":
            extragtis = gtiname3
        else:
            extragtis = extragtis + ',{0}'.format(gtiname3)

    # Make final merged GTI using nimaketime
    gtiname_merged = path.join(pipedir,"tot.gti")


    # Manage extra_expr for nimaketime (ST_VALID, DARK/DAY, and FPM_OVER_ONLY filters from --KEITH)
    list_extra_expr = ['ST_VALID.eq.1']

    if args.dark:
        list_extra_expr.append('SUNSHINE.eq.0')
    if args.day:
        list_extra_expr.append('SUNSHINE.eq.1')

    if args.keith:
        list_extra_expr.append('FPM_OVERONLY_COUNT<1')
        list_extra_expr.append('FPM_OVERONLY_COUNT<(1.52*COR_SAX**(-0.633))')
        #list_extra_expr.append('FPM_UNDERONLY_COUNT<200')
        if has_KP:
            list_extra_expr.append('(COR_SAX.gt.(1.914*KP**0.684+0.25)).and.KP.lt.5')

    if args.kpmax and has_KP:
        list_extra_expr.append('KP.lt.{0}'.format(args.kpmax))

    if args.minsun:
        # Exclude data that is at a Sun angle less that some value, unless it is in eclipse
        list_extra_expr.append('(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)'.format(args.minsun))

    extra_expr = "(" + " && ".join("%s" %expr for expr in list_extra_expr) + ")"

    cor_string="-"
    if args.cormin is not None:
        cor_string = "{0}-".format(args.cormin)
## Default is now 20 and 30 since that is standard processing as of 2019 May.  Option now ignored
    elvcut = 20.0
    brcut = 30.0
#    if args.shrinkelvcut:
#        # Keith suggests that these cuts can give more data without hurting data quality
#        elvcut = 20.0
#        brcut = 30.0
    maxunder = 200.0
    if args.nounderfilt:
        maxunder=650.0
    cmd = ["nimaketime",  "infile={0}".format(mkfile),
        'outfile={0}'.format(gtiname_merged), 'nicersaafilt=YES',
        'saafilt=NO', 'trackfilt=YES', 'ang_dist={0:.3f}'.format(args.angdist), 'elv={0}'.format(elvcut),
        'br_earth={0}'.format(brcut), 'cor_range={0}'.format(cor_string), 'min_fpm={0}'.format(args.minfpm),
        'underonly_range=0-{0}'.format(maxunder),
        'ingtis={0}'.format(extragtis), "clobber=yes", 
        'expr={0}'.format(extra_expr),
        'outexprfile={0}'.format(path.join(pipedir,"psrpipe_expr.txt"))]
    runcmd(cmd)

    # nimaketime gives an output GTI file with 1 row and START==STOP in the
    # case where no good time is selected.  This differs from the normal
    # maketime, which produces a GTI file with no rows in that case

    ###  Extract filtered events and apply merged GTI
    filteredname = path.join(pipedir,"cleanfilt.evt")
    intermediatename = path.join(pipedir,"intermediate.evt")

    # Build input file for niextract-events
    evlistname=path.join(pipedir,'evfiles.txt')
    fout = open(evlistname,'w')
    for en in evfiles:
        print('{0}'.format(en),file=fout)
    fout.close()

    # Build selection expression for niextract-events
    # Select events with PI in the selected range, require SLOW trigger (FAST optional)
    # and filter all undershoot, overshoot, and force trigger events
    evfilt_expr = 'PI={0}:{1},EVENT_FLAGS=bx1x000'.format(
        int(args.emin*KEV_TO_PI), int(args.emax*KEV_TO_PI))

    cmd = ["niextract-events", "filename=@{0}[{1}]".format(evlistname,evfilt_expr),
        "eventsout={0}".format(intermediatename), "timefile={0}".format(gtiname_merged),
        "gti=GTI", "clobber=yes"]
    runcmd(cmd)

    # Here we can automatically mask detectors by parsing the intermediate file
    bad_dets = None
    if args.mask is not None and args.mask[0] < 0:
        etable = Table.read(intermediatename,hdu=1)
        # Apply TIMEZERO if needed
        if 'TIMEZERO' in etable.meta:
            log.info('Applying TIMEZERO of {0} to etable'.format(etable.meta['TIMEZERO']))
            etable['TIME'] += etable.meta['TIMEZERO']
            etable.meta['TIMEZERO'] = 0.0
        log.info('Auto-masking detectors')
        bad_dets = find_hot_detectors(etable)
        if bad_dets is not None:
            log.info('Found hot detectors {0}!!'.format(bad_dets))
        # Make intermediate eng plot to show bad detectors
        cmd = ["nicerql.py", "--save",
               "--eng", intermediatename,
               "--lcbinsize", "{}".format(args.lcbinsize),
               "--filterbinsize", "{}".format(args.filterbinsize),
               "--basename", path.join(pipedir,basename)+"_intermediate"]

        if args.keith:
            cmd.append("--keith")

        runcmd(cmd)

    # Now filter any bad detectors
    evfilt_expr = '(EVENT_FLAGS==bx1x000)'
    # Filter any explicitly specified masked detectors
    if args.mask is not None:
        for detid in args.mask:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)
    # Filter any automatically identified hot detectors
    if bad_dets is not None:
        fout = open(path.join(pipedir,"bad_dets.txt"),"w")
        print("{0}".format(bad_dets),file=fout)
        fout.close()
        for detid in bad_dets:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)

    cmd = ["ftcopy", "{0}[{1}]".format(intermediatename,evfilt_expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
    # Remove intermediate file
    #os.remove(intermediatename)

    # Check that there are events left after hot detector cut
    if len(Table.read(filteredname,hdu=1))==0:
        log.error('No events left after hot detector filtering!')
        continue

    # Make cleanfilt.mkf file from ObsID .mkf file and merged_GTI
    cleanfilt_mkf = path.join(pipedir,"cleanfilt.mkf")
    log.info('Applying the GTI filtering to the *mkf file')
    cmd = ["fltime", "infile={}[1]".format(mkfile), "gtifile={}".format(gtiname_merged),"outfile={}".format(cleanfilt_mkf),"clobber=yes"]
    runcmd(cmd)

    # Make final clean plot
    #Dom is commenting this out 6/7 
    """
    cmd = ["nicerql.py", "--save",
           "--orb", path.join(pipedir,path.basename(orbfile)),
           "--map",
           "--sci", "--eng", filteredname,"--allspec","--alllc",
           "--lcbinsize", "{}".format(args.lcbinsize),
           "--filterbinsize", "{}".format(args.filterbinsize),
           "--mkf", cleanfilt_mkf, "--bkg",
           "--basename", path.join(pipedir,basename)+"_cleanfilt"]
    
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    if args.keith:
        cmd.append("--keith")
    runcmd(cmd)
    """
    # Add phases
    if args.par is not None:
        plotfile = path.join(pipedir,"phaseogram.png")
        log.info('Applying phases to {0}'.format(filteredname))
        log.info('Event file has TIMZERO < 0, so not applying --fix in photonphase!')
        #Dom removed the plot 6/8
        #cmd = ["photonphase", "--ephem", args.ephem, "--orb", orbfile, 
        #       "--addphase", filteredname, args.par]
        #runcmd(cmd)
        
        print('hi hi hi hiii running barycorr now')

        cmd = ["barycorr", filteredname, "%s/cleanfilt2.evt"%pipedir, orbfile]
        runcmd(cmd)
        
        print('hellooooo running photonphase3')

        etable = Table.read(intermediatename,hdu=1)
        if len(etable['TIME'] < 1000000):
            cmd = ["photonphase3.py", '--polycos', '--counts', '100000', 
                   '--overwrite', '%s/cleanfilt2.evt'%pipedir, args.par]
        else:
            cmd = ["photonphase3.py", '--polycos', '--overwrite', 
                   '%s/cleanfilt2.evt'%pipedir, args.par]
        runcmd(cmd)

        #Dom changing to split photonphase
        #if args.crab is True, use multiprocessing
        #pipeline_utils.split_photonphase(
        #        filteredname, orbfile, args.par, use_mp=args.crab)


    # Extract simple PHA file and light curve
    phafile = path.splitext(filteredname)[0] + ".pha"
    lcfile = path.splitext(filteredname)[0] + ".lc"
    cmd = ["extractor", filteredname, "eventsout=none", "imgfile=none",
        "phafile={0}".format(phafile), "fitsbinlc={0}".format(lcfile),
        "binlc={}".format(args.lcbinsize), "regionfile=none", "timefile=none",
        "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "xcolh=RAWX",
        "ycolh=RAWY", "gti=GTI"]
    runcmd(cmd)

    # if --merge option, Add clean evt file to list of files to merge
    if args.merge:
        all_evfiles.append(filteredname)
        all_orbfiles.append(orbfile)


# Merging all ObsIDs
if args.merge and (len(all_evfiles)>1) :
    ## save list of orbit files
    orbfiles_list = path.join(os.getcwd(),"list_orbfiles.txt")
    np.savetxt(orbfiles_list,all_orbfiles,fmt=['%s'])

    ## Call merge.py
    cmd = ["merge.py"]
    for evt in all_evfiles:
        cmd.append(evt)
    cmd.append("merged")
    cmd.append("merged_{0}".format(args.outdir))
    cmd.append("--clobber")
    cmd.append("--lcbinsize")
    cmd.append("{}".format(args.lcbinsize))

    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
        cmd.append("--orb")
        cmd.append("{0}".format(orbfiles_list))
    if args.crcut:
        cmd.append("--cut")
        cmd.append("--filterbinsize")
        cmd.append("{}".format(args.filterbinsize))

    # if args.crabnorm:
    #     cmd.append("--crabnorm")
    # if args.dark:
    #     cmd.append("--dark")
    runcmd(cmd)

else:
    if args.crcut:
        log.warning("Count rate cuts are only performed on merged event files (add the option --merge)")

    # # Perform Crab Normalization (in the case of single ObsID bring processed)
    # --- NOT IMPLEMENTED --- #
    # if args.crabnorm:
    #     log.info("Performing normalization of the cleanfilt spectrum with crab residuals")
    #     normed_spec = path.join(phafile.strip('.pha'),"_crabnorm.pha")
    #     if args.dark:
    #         cmd = ["mathpha","{}/{}".format(phafile,CRAB_RES_NIGHT),"R","{}".format(normed_spec),"{} % POISS-0 0".format(phafile)]
    #         runcmd(cmd)
    #     else:
    #         cmd = ["mathpha","{}/{}".format(phafile,CRAB_RES_TOT),"R","{}".format(normed_spec),"{} % POISS-0 0".format(phafile)]
    #         runcmd(cmd)

shutil.rmtree(tempdir)
