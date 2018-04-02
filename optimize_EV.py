 #!/usr/bin/env python
from __future__ import (print_function, division, unicode_literals, absolute_import)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import scipy
import argparse
from astropy import log
from glob import glob
from nicer.values import *
from os import path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from pint.eventstats import hm

parser = argparse.ArgumentParser(description="Read event file with phase and optimize cuts")
parser.add_argument("evfiles", help="Name of event files to process", nargs='+')
args = parser.parse_args()

for ef in args.evfiles:
    tlist = []
    for fn in args.evfiles:
        log.info('Reading file {0}'.format(fn))
        tlist.append(Table.read(fn,hdu=1))
    log.info('Concatenating files')
    if len(tlist) == 1:
        etable = tlist[0]
    else:
        etable = vstack(tlist,metadata_conflicts='silent')
    del tlist

    phasesinitial = etable['PULSE_PHASE']
    #TOA_accuracy = get_TOA_accuracy(phasesinitial)
    eminbest = 0.0
    emaxbest = 100.0
    #print("Initial Htest = {0}".format(TOA_accuracy))


    emins = np.arange(0.25,4.0,0.025)
    for emin in emins:
        emaxs = np.arange(emin,12.0,0.1)
        for emax in emaxs:
            idx = np.where(np.logical_and(etable['PI']*PI_TO_KEV>emin,
                etable['PI']*PI_TO_KEV<emax))[0]
            phases = etable['PULSE_PHASE'][idx]
            print("Energy Range: {} -> {}".format(emin,emax))
	    #if len(phases)==0:
            #    continue
            #TOA_accuracy = get_TOA_accuracy(phases)
            #if TOA_accuracy >= TOA_best:
            #    TOA_best= TOA_accuracy
            #    eminbest=emin
            #    emaxbest=emax

