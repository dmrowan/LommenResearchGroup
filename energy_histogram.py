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

#This script will take in an evt file and produce a histogram 
desc = """
This creates a nergy histogram for one event file.
"""

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--evt", help="Input event file .evt")
parser.add_argument("--max", help="Input max eV", type=float, default=15)
parser.add_argument("--min", help="Input min eV", type=float, default=0)
parser.add_argument("--bin", help="Input number of bins", type=int, default=15)
parser.add_argument("--phase", help="True for additional phase historgram", type=bool, default=False)
args = parser.parse_args()

bkfstable = Table.read(args.evt,hdu=1)
#print (bkfstable.info)
energy_pi = bkfstable['PI']

range_of_values = (args.min, args.max)
energy_bins = args.bin

if not args.phase:
	energy_pi = energy_pi *(1.0/100)
	n, bins, patches = plt.hist(energy_pi, range=range_of_values, bins=energy_bins)
	plt.title(str(args.evt)+ " Historgram")
	plt.xlabel("Energy [keV]")
	plt.ylabel("Number of photons")
	plt.show()
else:
	plt.subplot(211)	
	energy_pi = energy_pi *(1.0/100)
	n, bins, patches = plt.hist(energy_pi, range=range_of_values, bins=energy_bins)
	plt.title(str(args.evt)+ " Historgram")
	plt.ylabel("Number of photons")
	plt.xlabel("Energy [keV]")
	phase=bkfstable['PULSE_PHASE']
	plt.subplot(212)
	plt.hist(phase, bins = args.bin, histtype = "step")
	plt.xlabel("Phase [deg]")
	plt.title("Pulse Profile")
	plt.ylabel("Number of photons")
	plt.show()


