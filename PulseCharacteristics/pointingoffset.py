#!/usr/bin/env python

#imports from nicerql.py
from __future__ import (print_function, division, unicode_literals, absolute_import)
import sys
# Hack to add this to pythonpath
#sys.path.append('/Users/paulr/src/NICERsoft')
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.stats import mad_std, sigma_clipped_stats
from astropy.time import Time
from os import path
from nicer.values import *
from nicer.cartographer import *
from nicer.plotutils import *
from nicer.sci_plots import sci_plots
from nicer.eng_plots import eng_plots
from nicer.eng_plots import plot_all_spectra
from nicer.eng_plots import plot_all_lc
from glob import glob

#imports from plotutils.py
import numpy as np
import matplotlib.pyplot as plot
import matplotlib as mpl
import copy
import scipy
from scipy import ndimage
from astropy import log
from glob import glob
from nicer.values import *
from os import path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std, sigma_clipped_stats
import astropy.io.fits as pyfits
import astropy.units as u

from argparse import Namespace
import pandas as pd
from nicer.NicerFileSet import *

from nicer.plotutils import *
from nicer.sci_plots import sci_plots
from nicer.eng_plots import eng_plots
from nicer.eng_plots import plot_all_spectra
from nicer.eng_plots import plot_all_lc
from glob import glob
import sys
from nicer.bkg_plots import *
from nicer.fitsutils import *
from InteractiveLC import *
from nicer.NicerFileSet import *

#get filenames
fnames = pd.read_csv('crabfilenames2.txt', header = None)
print(fnames[0][0])
fname = str(fnames[0][0])
print(fname)

#arguments
args = Namespace(alllc=False, allspec=False, applygti=None, basename=None, bkg=True, emax=-1.0, emin=-1.0, eng=False, eventshootrate=False, extraphkshootrate=False, filtall=False, filterbinsize=16.0, filtovershoot=False, filtratio=False, filtswtrig=False, filtundershoot=False, foldfreq=0.0, gtirows=None, infiles=[], interactive=False,keith=False, lcbinsize=1.0, lclog=False, map=False, mask=None, merge=False, mkf=None, nyquist=100.0, object=None, obsdir=fname, orb=None, par=None, pi=False, powspec=False, pslog=False, save=False, sci=False, sps=None, tskip=0.0, useftools=False, writebkf=False, writeps=False)

#from nicerql.py
data = NicerFileSet(args)
etable = data.etable
gtitable = data.gtitable
mktable = data.mktable
basename = data.basename
ovbintable = data.ovbintable
if args.pi or not ('PI' in etable.colnames):
    log.info('Adding PI')
    calfile = path.join(datadir,'gaincal_linear.txt')
    pi = calc_pi(etable,calfile)
    etable['PI'] = pi
if gtitable is not None:
    startmet = gtitable['START'][0]
    stopmet = gtitable['STOP'][0]
    duration = stopmet-startmet
    # Add 1 bin to make sure last bin covers last events
    lc_elapsed_bins = np.arange(0.0,duration+args.lcbinsize,args.lcbinsize)
    lc_met_bins = startmet+lc_elapsed_bins
    cumtimes = [ 0.0 ]
    cumtime = lc_elapsed_bins[-1]+args.lcbinsize
    for i in range(1,len(gtitable['START'])):
        startmet = gtitable['START'][i]
        stopmet = gtitable['STOP'][i]
        duration = stopmet-startmet
        myelapsedbins = np.arange(0.0,duration+args.lcbinsize,args.lcbinsize)
        lc_elapsed_bins = np.append(lc_elapsed_bins,cumtime+myelapsedbins)
        lc_met_bins = np.append(lc_met_bins,np.arange(startmet,stopmet+args.lcbinsize,args.lcbinsize))
        mylcduration = myelapsedbins[-1]+args.lcbinsize
        cumtimes.append(cumtime)
        cumtime += mylcduration
    gtitable['CUMTIME'] = np.array(cumtimes)
if args.bkg:
    #if hkmet is None:
    #    log.error("Can't make background plots without MPU HKP files")
    #else:
    #     if eventovershoots is not None:
    #         figure4 = bkg_plots(etable, data, gtitable, args, mktable, data.eventshoottable)
    #     else:
    #         figure4 = bkg_plots(etable, data, gtitable, args, mktable, data.hkshoottable)

   # figure4 = bkg_plots(etable, gtitable, args, mktable, ovbintable)
   # figure4.set_size_inches(16,12)
    if args.save:
        log.info('Writing bkg plot {0}'.format(basename))
        figure4.savefig('{0}_bkg.png'.format(basename), dpi = 100)

#from plotutils.py
def plot_total_count_hist(etable, ax_rate, ax_counts):
    'Plots event count per ID as a histogram with event count and countrate on y axes'

    num_events, colors = hist_use(etable)

    tc = ax_counts.bar(IDS, num_events, color = colors)

    ax_rate.set_ylabel('c/s')
    cntmin, cntmax = ax_counts.get_ylim()
    ax_rate.set_ylim((cntmin/etable.meta['EXPOSURE'],cntmax/etable.meta['EXPOSURE']))
    #countrate.set_ylim([np.min(rate)-20,np.max(rate)+20])

    ax_counts.set_xlabel('DET_ID')
    ax_counts.set_ylabel('# of Events')
    plot.locator_params(nticks = 20)
    plot.title('Total (Filtered) Event Count by Detector')
    #total_counts.set_ylim([np.min(num_events)-20, np.max(num_events)+20])

    return num_events

def gti_colormap():
    colornames = ['black','green','red','blue','magenta','orange','cyan','yellow','gray']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')
    return colornames, cmap, norm

#pointing, from plotutils.py
time, pointing, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['ANG_DIST'], gtitable)
colornames, cmap, norm = gti_colormap()
plot.scatter(time, pointing, c = np.fmod(cc,len(colornames)), cmap = cmap,norm=norm, marker = '+', label='Pointing Offset')
plot.legend(loc = 2)
plot.ylabel('Angle (deg)')
plot.yscale('log')
plot.axhline(1.0/60.0,c='r')
plot.ylim((0.0001,100.0))

rows = np.where((pointing >= (1/60 - 0.0005)) & (pointing <= (1/60 + 0.0005)))
times = time[rows]
t = []
for n in times:
    t.append(n + mktable['TIME'][0])
t2 = []
for n in times:
    t2.append(n + gtitable['START'][0])
print(t)
print(t2)
plot.show()
