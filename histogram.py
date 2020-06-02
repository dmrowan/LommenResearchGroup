#!/usr/bin/env python

from __future__ import print_function, division
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
from astropy.table import Table
import sys
from nicer.values import *
import argparse

desc = "This creates both a phase, and an energy histogram for an evt file. (In order to make the phase histogram, it will need to have been run through photon_phase. I'll at some point implement a check to see if the column exists so that it can be run with or without the addition of the pulse profile)"

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--evt",  help = "The evt file that you'd like to graph")
parser.add_argument("--min",  help = "Minimum energy range to be considered", type = float, default = .25)
parser.add_argument("--max",  help = "Maximum energy range to be considered", type = float, default = 0)
parser.add_argument("--bin",  help = "Number of bins to use in the histogram", type = int, default =20)

args = parser.parse_args()


bkfstable = Table.read(args.evt,hdu=1)
#print (bkfstable.info)
energy = bkfstable['PI']/100.  #the units are wrong in the energy column, so corrected by 1/100
phase = bkfstable['PULSE_PHASE']

left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.8   # the amount of width reserved for space between subplots,
               # expressed as a fraction of the average axis width
hspace = 0.8   # the amount of height reserved for space between subplots,
               # expressed as a fraction of the average axis height
plt.subplots_adjust(left, bottom, right, top,
                wspace, hspace)

plt.figure(1)

plt.subplot(211)
energy_hist = plt.hist(energy, bins = args.bin, range= (args.min,args.max))
plt.title(str(args.evt) + " " + "Histogram")
plt.ylabel("Number of photons")
plt.xlabel("Photon energy (kev)")

plt.subplot(212)
plt.hist(phase, bins = args.bin, histtype = "step")
plt.title("pulse_profile")
plt.ylabel("Number of photons")
plt.xlabel("pulse phase")

#plt.savefig(str(sys.argv[1][:-4]) + ".png")
#plt.plot(pi, phase)
plt.show()

