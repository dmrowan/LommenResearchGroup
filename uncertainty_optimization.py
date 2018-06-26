#!/usr/bin/env python

# Yet another copy of photon_toa.py for me to mess around with.

from __future__ import division, print_function, absolute_import
import os,sys
import argparse
import numpy as np
from astropy import log
import astropy.units as u
import astropy.io.fits as pyfits
import pint.residuals
from pint.event_toas import load_NICER_TOAs
from pint.event_toas import load_RXTE_TOAs
from pint.event_toas import load_XMM_TOAs
from pint.plot_utils import phaseogram_binned
from pint.observatory.nicer_obs import NICERObs
from pint.observatory.rxte_obs import RXTEObs
import pint.toa, pint.models
from pint.eventstats import hmw, hm, h2sig
from astropy.time import Time, TimeDelta
from pint.templates.lctemplate import LCTemplate,prim_io
from pint.templates import lcfitters
import cPickle
import cStringIO
from collections import deque
import time

import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import scipy
from astropy import log
from glob import glob
from nicer.values import *
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord






def local_load_NICER_TOAs(eventname):
    """ Local override to add MET field to each TOA object."""
    # TODO -- add this to PINT method ?
    tl  = load_NICER_TOAs(eventname)
    f = pyfits.open(eventname)
    mets = f['events'].data.field('time')
    f.close()
    for t,met in zip(tl,mets):
        t.met = met
    return tl

desc="""Generate TOAs from photon event data."""

parser=argparse.ArgumentParser(description=desc)
parser.add_argument("eventname",help="FITS file to read events from")
parser.add_argument("templatename",help="Name of file to read template from")
parser.add_argument("parname",help="Timing model file name")
parser.add_argument("--orbfile",help="Name of orbit file", default=None)
parser.add_argument("--ephem",help="Planetary ephemeris to use (default=DE421)", default="DE421")
parser.add_argument("--plot",help="Show phaseogram plot.", action='store_true', default=False)
parser.add_argument("--plotfile",help="Output figure file name (default=None)", default=None)
parser.add_argument("--fitbg",help="Fit an overall background level (e.g. for changing particle background level (default=False).",action='store_true',default=False)
parser.add_argument("--unbinned",help="Fit position with unbinned likelihood.  Don't use for large data sets. (default=False)",action='store_true',default=False)
#parser.add_argument("--fix",help="Adjust times to fix 1.0 second offset in NICER data (default=False)", action='store_true',default=False)
parser.add_argument("--tint",help="Integrate for tint seconds for each TOA, or until the total integration exceeds maxint.  The algorithm is based on GTI, so the integration will slightly exceed tint (default None; see maxint.)",default=None)
parser.add_argument("--maxint",help="Maximum time interval to accumulate exposure for a single TOA (default=2*86400s)",default=2*86400.)
parser.add_argument("--minexp",help="Minimum exposure (s) for which to include a TOA (default=None).",default=None)

## Parse arguments
args = parser.parse_args()

# Load PINT model objects
modelin = pint.models.get_model(args.parname)
log.info(str(modelin))

# check for consistency between ephemeris and options
if 'SolarSystemShapiro' in modelin.components.keys():
    planets=True
else:
    planets=False

# Load Template objects
try:
    template = cPickle.load(file(args.templatename))
except:
    primitives,norms = prim_io(args.templatename)
    template = LCTemplate(primitives,norms)
print(template)

# Load photons as PINT toas, and weights, if specified
# Here I might loop over the files specified
# Read event file header to figure out what instrument is is from
hdr = pyfits.getheader(args.eventname,ext=1)

log.info('Event file TELESCOPE = {0}, INSTRUMENT = {1}'.format(hdr['TELESCOP'],
    hdr['INSTRUME']))
if hdr['TELESCOP'] == 'NICER':
    # Instantiate NICERObs once so it gets added to the observatory registry
    if args.orbfile is not None:
        log.info('Setting up NICER observatory')
        NICERObs(name='NICER',FPorbname=args.orbfile,tt2tdb_mode='none')
    # Read event file and return list of TOA objects
    try:
        tl  = local_load_NICER_TOAs(args.eventname)
    except KeyError:
        log.error('Failed to load NICER TOAs. Make sure orbit file is specified on command line!')
        raise
elif hdr['TELESCOP'] == 'XTE':
    # Instantiate RXTEObs once so it gets added to the observatory registry
    if args.orbfile is not None:
        # Determine what observatory type is.
        log.info('Setting up RXTE observatory')
        RXTEObs(name='RXTE',FPorbname=args.orbfile,tt2tdb_mode='none')
    # Read event file and return list of TOA objects
    tl  = load_RXTE_TOAs(args.eventname)
elif hdr['TELESCOP'].startswith('XMM'):
    # Not loading orbit file here, since that is not yet supported.
    tl  = load_XMM_TOAs(args.eventname)
else:
    log.error("FITS file not recognized, TELESCOPE = {0}, INSTRUMENT = {1}".format(
        hdr['TELESCOP'], hdr['INSTRUME']))
    sys.exit(1)

# Now convert to TOAs object and compute TDBs and posvels
ts = pint.toa.get_TOAs_list(tl,ephem=args.ephem,planets=planets,include_bipm=False,include_gps=False)
del tl
ts.filename = args.eventname
#if args.fix:
#	if hdr['TIMEZERO'] < 0.0:
#        log.error('TIMEZERO<0 and --fix: You are trying to apply the 1-s offet twice!')
#    ts.adjust_TOAs(TimeDelta(np.ones(len(ts.table))*-1.0*u.s,scale='tt'))

print(ts.get_summary())
mjds = ts.get_mjds()
print(mjds.min(),mjds.max())

# Compute model phase for each TOA
phss = modelin.phase(ts.table)[1].value # discard units
# ensure all postive
phases = np.where(phss < 0.0, phss + 1.0, phss)
tdbs = ts.table['tdb']
h = float(hm(phases))
print("Htest : {0:.2f} ({1:.2f} sigma)".format(h,h2sig(h)))
if args.plot:
    phaseogram_binned(mjds,phases,bins=100,plotfile = args.plotfile)

# get exposure information
try:
    f = pyfits.open(args.eventname)
    exposure = f[1].header['exposure']
    f.close()
except:
    exposure = 0


def estimate_toa(mjds,phases,tdbs):
    """ Return a pint TOA object for the provided times and phases."""

    # Given some subset of the event times, phases, and weights, compute
    # the TOA based on a reference event near the middle of the span.
    # Build the TOA as a PINT TOA() object
    lcf = lcfitters.LCFitter(template,phases)
    if args.fitbg:
        for i in xrange(2):
            lcf.fit_position(unbinned=False)
            lcf.fit_background(unbinned=False)
    dphi,dphierr = lcf.fit_position(unbinned=args.unbinned)
    log.info('Measured phase shift dphi={0}, dphierr={1}'.format(dphi,dphierr))

    # find MJD closest to center of observation and turn it into a TOA
    argmid = np.searchsorted(mjds,0.5*(mjds.min()+mjds.max()))
    tmid = tdbs[argmid]
    tplus = tmid + TimeDelta(1*u.s,scale='tdb')
    toamid = pint.toa.TOA(tmid)
    toaplus = pint.toa.TOA(tplus)
    toas = pint.toa.TOAs(toalist=[toamid,toaplus])
    toas.compute_TDBs()
    toas.compute_posvels(ephem=args.ephem,planets=planets)
    phsi,phsf = modelin.phase(toas.table)
    fbary = (phsi[1]-phsi[0]) + (phsf[1]-phsf[0])
    fbary._unit = u.Hz
    # First delta is to get time of phase 0.0 of initial model
    # Second term corrects for the measured phase offset to align with template
    print("printing dphi: {}".format(dphi))
    print("printing fbary: {}".format(fbary))
    print("printing TimeDelta(dphi/fbary,scale='tdb'): {}".format(TimeDelta(dphi/fbary,scale='tdb')))
    tfinal = tmid + TimeDelta(-phsf[0].value/fbary,scale='tdb') + TimeDelta(dphi/fbary,scale='tdb')
    print("printing tfinal from trial5_uncertainty.py: {}".format(tfinal))
    # Use PINT's TOA writer to save the TOA
    nsrc = lcf.template.norm()*len(lcf.phases)
    nbkg = (1-lcf.template.norm())*len(lcf.phases)
    toafinal = pint.toa.TOA(tfinal,
            nsrc='%.2f'%nsrc,nbkg='%.2f'%nbkg,exposure='%.2f'%exposure,dphi='%.5f'%dphi)
    log.info("Src rate = {0} c/s, Bkg rate = {1} c/s".format(nsrc/exposure, nbkg/exposure))
    print("dphierr = {0}, fbary = {1}".format(dphierr, fbary))
    return toafinal, dphierr/fbary*1e6

def make_toa(mjds, phases, tdbs):
	if args.tint is None:
	    # do a single TOA for table
	    toafinal,toafinal_err = estimate_toa(mjds,phases,tdbs)
	    print("hopefully this is an astropy quantity: {0}".format(toafinal_err))
	    if 'OBS_ID' in hdr:
		# Add ObsID to the TOA flags
		toafinal.flags['-obsid'] = hdr['OBS_ID']
	    toafinal = [toafinal]
	    toafinal_err = [toafinal_err]
	for t in toafinal:
	    t.flags["-t"] = hdr['TELESCOP']
	toas = pint.toa.TOAs(toalist=toafinal)
	print("--------------------------------------------------")
	print("Print statements in trial5_uncertainty.py")
	toas.table['error'][:] = np.asarray(toafinal_err)
	print("Trying after toas.table['error'][:]")
	print(toas.table['error'])
	print("--------------------------------------------------")
	sio = cStringIO.StringIO()
	toas.write_TOA_file(sio,name='nicer',format='tempo2')
	output = sio.getvalue()
	output = output.replace('barycenter','@')
	print("This is right before printing output variable")
	print(output)

def get_TOA_accuracy(phases):
    # Given some subset of the event times, phases, and weights, compute
    # the TOA based on a reference event near the middle of the span.
    # Build the TOA as a PINT TOA() object
    lcf = lcfitters.LCFitter(template,phases)
    if args.fitbg:
        for i in xrange(2):
            lcf.fit_position(unbinned=False)
            lcf.fit_background(unbinned=False)
    dphi,dphierr = lcf.fit_position(unbinned=args.unbinned)
    #print("The TOA error is: {}".format(dphierr/fbary)) 
    return (dphierr/fbary)



etable = Table.read(args.eventname,hdu=1)
# Compute model phase for each TOA
phss = modelin.phase(ts.table)[1].value # discard units
# ensure all postive
phasesinitial = np.where(phss < 0.0, phss + 1.0, phss)

# TODO Test with original phases because my final TOA does not seem to be working
print("------------------------------------------")
print("INITIAL TOA")
#make_toa(mjds, phasesinitial, tdbs)
print("------------------------------------------")

TOA_best = estimate_toa(mjds,phasesinitial,tdbs)
eminbest = 0.0
emaxbest = 100.0
print("Initial TOA accuracy: {}".format(TOA_best))
emins = np.arange(0.25, 4.00, 0.02)
start_time = time.time()
#with log.log_to_file('myfile.log'):	
for emin in emins:
    emaxs = np.arange(emin+0.10, 12.01, 0.10)
    for emax in emaxs:    # currently etable doesn't exist, and PULSE_PHASE wouldn't be in it either way
	idx = np.where(np.logical_and(etable['PI']*PI_TO_KEV>emin, etable['PI']*PI_TO_KEV<emax))[0]
	phases = etable['PULSE_PHASE'][idx]
	if len(phases)==0:
		continue
	loop_start = time.time()
	TOA, TOA_error = estimate_toa(mjds, phases, tdbs)
	print("Energy Range: {0} -> {1} gives: {2}".format(emin, emax, TOA_error))
	loop_end = time.time()
	print("time elapsed per iteration: {}".format(loop_end - loop_start))
	#Converting from astropy quantity to python scalar for operand
	if float(TOA_error*(1.0*u.Hz)) <= TOA_best:
		TOA_best= float((TOA_error*(1.0*u.Hz)))
		eminbest=emin
		emaxbest=emax
		print("New Best Energy Range: {0} -> {1}, gives {2}".format(eminbest, emaxbest, TOA_best))
idx = np.where(np.logical_and(etable['PI']*PI_TO_KEV>eminbest, etable['PI']*PI_TO_KEV<emaxbest))[0]
phases = etable['PULSE_PHASE'][idx]
print("Double check which mask are we using: {} -> {}".format(eminbest,emaxbest))
#print(ts.get_summary())
#mjds = ts.get_mjds()
#print(mjds.min(),mjds.max())

#tdbs = ts.table['tdb']
#toa_final, toa_final_error = estimate_toa(mjds, phases, tdbs)
make_toa(mjds, phases, tdbs)
final_time = time.time()
print("time elapsed: {}".format(final_time - start_time))
print("Final TOA_Accuracy = {0}".format(TOA_best*u.us))
print("Best emin: {0}, emax: {1}".format(eminbest,emaxbest))
