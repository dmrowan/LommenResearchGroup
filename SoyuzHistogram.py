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
from datetime import datetime
from datetime import timedelta
from dateutil import parser as pr

desc = """
This looks at higher energy photons as a function of time for a combined evt file
"""

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--evt", help="Input event file .evt")
parser.add_argument("--min", help="Input min eV", type=float, default=10)
parser.add_argument("--bin", help="Number of time bins", type=int, default=50)
args = parser.parse_args()

evttable  = Table.read(args.evt,hdu=1)
energy_pi = evttable['PI']
energy_pi = energy_pi * (1.0/100)

times = evttable['TIME']

firsttime = evttable['TIME'][0]
lasttime = evttable['TIME'][len(evttable['TIME'])-1]

n, timebins, patches = plt.hist(times, bins=args.bin)
plt.close()
high_energy_indicies = np.where(evttable["PI"] > (args.min * 100))[0]
n_energy_binned = []
for i in range(len(timebins)-1):
	if i == 0:
		numinbin = 0 
	else:	
		numinbin = len(np.where(evttable['TIME'][high_energy_indicies] < timebins[i])[0]) - len(np.where(evttable['TIME'][high_energy_indicies] < timebins[i-1])[0])
	
	n_energy_binned.append(numinbin)

print(n_energy_binned)

firsttime = datetime(year=2014, month=1, day=1, hour=0, minute=0, second=0)
timebins_corrected = []
for time in timebins:
	bettertime = firsttime + timedelta(seconds=time)
	timebins_corrected.append(bettertime)

print(timebins_corrected)

plt.bar(timebins_corrected[1:], n_energy_binned, fill=False)
plt.xlabel("Time?")
plt.ylabel("Number of Photons Above Threshold")
plt.show()
