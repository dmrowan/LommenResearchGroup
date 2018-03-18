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
import itertools

desc = """
This looks at higher energy photons as a function of time for a group of evt files in the current directory
"""

parser = argparse.ArgumentParser(description = desc)
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

#Create an empty directory that will have our time, filename, exposure, and number of photons > threshold
files_unordered = {}

#Loop through the files in the current directory. This doesn't go sequentially. 
for filename in os.listdir(os.getcwd()):
	if os.path.isfile(filename+"/cleanfilt.evt"):
        	#print(filename)
		#First set the path to the evt file
                cleanfilt = filename + "/cleanfilt.evt"
		
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		high_energy_indicies = np.where(evttable_hdu1["PI"] > (args.min * 100))[0]
		nphotons = len(high_energy_indicies)

		#Read in HDU 2 to get the exposure info
                evttable_hdu2 = Table.read(cleanfilt, hdu=2)
                obstime = evttable_hdu2.meta['DATE-OBS']
                obstime_dt = pr.parse(obstime)
                exposure = sum(evttable_hdu2["STOP"] - evttable_hdu2["START"])
                 
                files_unordered[obstime_dt] = [filename, exposure, nphotons]

#We will fill these lists when we sort the directory
dates_ordered = []
filename_ordered = []
exposure_ordered = []
nphotons_ordered = []
#Use itertools to sort these four lists by time
for key in sorted(files_unordered.iterkeys()):
	#print(key, files_unordered[key])
        dates_ordered.append(key)
	filename_ordered.append(files_unordered[key][0])
        exposure_ordered.append(files_unordered[key][1])
	nphotons_ordered.append(files_unordered[key][2])

#Now that we have the ordered lists, we can find relative nphotons by dividing by exposure

normalized_photons = []

for i in range(len(nphotons_ordered)):
	if nphotons_ordered[i] == 0 or exposure_ordered[i] == 0:
		normalized_photons.append(0)
	else:
		normalized_photons.append(float(nphotons_ordered[i])/float(exposure_ordered[i]))

#Create a bar plot of the results
plt.bar(dates_ordered, normalized_photons, fill=False)
plt.xlabel('Date')
plt.ylabel('Relative nphotons above threshold')

#for i in range(len(soyuz_dates)): 
#	plt.axvline(x=soyuz_dates[i])

for tup in soyuz_dates_rassvet:
	plt.axvspan(tup[0], tup[1], alpha=0.25, color='red')

for tup in soyuz_dates_poisk:
	plt.axvspan(tup[0], tup[1], alpha=0.25, color='blue')

plt.axis([pr.parse('2017-06-01'), pr.parse('2017-12-01'), 0, 1])

plt.title('1821-24 High Energy Photons')
plt.savefig('1821_high_energy_photons.png')
plt.show()
