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

#Dom Rowan, 2018
#This is an outdated code

desc = """
This looks at higher energy photons as a function of time for a combined evt file
"""

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--evt", help="Input event file .evt")
parser.add_argument("--min", help="Input min eV", type=float, default=10)
parser.add_argument("--bin", help="Number of time bins", type=int, default=50)
args = parser.parse_args()

#Put all the Soyzu dates in a list
soyuz_dates = ['2017-06-02', '2017-07-28', '2017-09-02', '2017-09-13', '2017-12-14', '2017-12-17', '2017-12-19']
soyuz_dates = [ pr.parse(i) for i in soyuz_dates ]

soyuz_dates_rassvet = [('2016-11-19','2017-06-02'),('2017-07-28', '2017-12-14'), ('2017-12-19', '2018-03-18')]
soyuz_dates_poisk = [('2017-04-20', '2017-09-02'), ('2017-09-13', '2018-02-27')]

soyuz_dates_rassvet = [ (pr.parse(i[0]), pr.parse(i[1])) for i in soyuz_dates_rassvet ]
soyuz_dates_poisk = [ (pr.parse(i[0]),pr.parse(i[1])) for i in soyuz_dates_poisk ]

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
plt.xlabel("Date")
plt.ylabel("Number of Photons Above Threshold")
#for i in range(len(soyuz_dates)):
#        plt.axvline(x=soyuz_dates[i])

for tup in soyuz_dates_rassvet:
        plt.axvspan(tup[0], tup[1], alpha=0.25, color='red')

for tup in soyuz_dates_poisk:
        plt.axvspan(tup[0], tup[1], alpha=0.25, color='blue')

plt.axis([pr.parse('2017-06-01'), pr.parse('2017-12-01'), 0, 1])

plt.title("Number of Photons Above Threshold for 1821-24")

plt.show()
