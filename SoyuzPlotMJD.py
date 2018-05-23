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
from astropy.time import Time

#Dom Rowan, 2018
#This is the most current verison of the Soyuz Plot code

desc = """
This looks at photons in a specific range as a function of time for all evt files in the current directory. Ensure that your directory only contains the files 
for the selected pulsar. This version has the option to select a time range and use hour or minute  bins instead of observation date bins
"""

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--min", help="Input min eV", type=float, default=2)
parser.add_argument("--max", help="Input max eV", type=float, default=10)
parser.add_argument("--hbin", help="Use hour bins T/F", default=False, action='store_true')
parser.add_argument("--mbin", help="Use minute bins T/F", default=False, action='store_true')
parser.add_argument("--dmin", help="Enter start date for obsID range", type=str, default=None)
parser.add_argument("--dmax", help="Enter end date for obsID range", type=str, default=None)
args = parser.parse_args()

#If we use minute and hourbins, default to minute mins
if args.mbin and args.hbin:
	args.hbin = False

#Put all the Soyzu dates in a list
soyuz_dates = ['2017-06-02', '2017-07-28', '2017-09-02', '2017-09-13', '2017-12-14', '2017-12-17', '2017-12-19'] 
soyuz_dates = [ pr.parse(i) for i in soyuz_dates ]

soyuz_dates_rassvet = [('2016-11-19 21:58:00','2017-06-02 10:47:00'),('2017-07-28 21:54:00', '2017-12-14 05:14:00'), ('2017-12-19 08:39:00', '2018-04-01')] 
soyuz_dates_poisk = [('2017-04-20 13:18:00', '2017-09-02 21:58:00'), ('2017-09-13 02:55:00', '2018-02-27 23:08:00')]

soyuz_dates_rassvet = [ (pr.parse(i[0]), pr.parse(i[1])) for i in soyuz_dates_rassvet ]
soyuz_dates_poisk = [ (pr.parse(i[0]),pr.parse(i[1])) for i in soyuz_dates_poisk ]

for i in range(len(soyuz_dates_rassvet)):
	dastropy0 = Time(soyuz_dates_rassvet[i][0], format='datetime', scale='utc')
	dastropy0.format = 'mjd'
	dastropy1 = Time(soyuz_dates_rassvet[i][1], format='datetime', scale='utc')
	dastropy1.format = 'mjd'
	soyuz_dates_rassvet[i] = (dastropy0.value, dastropy1.value)

for i in range(len(soyuz_dates_poisk)):
	dastropy0 = Time(soyuz_dates_poisk[i][0], format='datetime', scale='utc')
	dastropy0.format = 'mjd'
	dastropy1 = Time(soyuz_dates_poisk[i][1], format='datetime', scale='utc')
	dastropy1.format = 'mjd'
	soyuz_dates_poisk[i] = (dastropy0.value, dastropy1.value)

#Create an empty directory that will have our time, filename, exposure, and number of photons > threshold
files_unordered = {}

#Loop through the files in the current directory. This doesn't go sequentially. 
for filename in os.listdir(os.getcwd()):
	if os.path.isfile(filename+"/cleanfilt.evt"):
		#First set the path to the evt file
                cleanfilt = filename + "/cleanfilt.evt"
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		obstime = pr.parse(evttable_hdu1.meta['DATE-OBS'])
		
		#If we have a max and min date, select obsIDs here
		if ((args.dmin is not None) and (args.dmax is not None)):
			if (obstime >= pr.parse(args.dmin)) and (obstime <= pr.parse(args.dmax)):
				files_unordered[obstime] = [filename]
                elif ((args.dmin is not None) and (args.dmax is None)):
			if (obstime >= pr.parse(args.dmin)):
				files_unordered[obstime] = [filename]
		elif ((args.dmin is None) and (args.dmax is not None)):
			if (obstime <= pr.parse(args.dmax)):
				files_unordered[obstime] = [filename]
		else:
			files_unordered[obstime]= [filename]

#This will contain the files in order
filename_ordered = []
dates_ordered = []
#Use itertools to sort these four lists by time
for key in sorted(files_unordered.iterkeys()):
	filename_ordered.append(files_unordered[key][0])
	dates_ordered.append(key)

#The case where we use ObsID bins
if (not args.hbin) and (not args.mbin):
	print("Using obsID bins")
	nphotons_ordered = []
	exposure_ordered = []
	#For each file in order
	for filename in filename_ordered:
		cleanfilt = filename + "/cleanfilt.evt"
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		#Find high energy indicies
		high_energy_indicies = np.where((evttable_hdu1["PI"] > (args.min * 100)) & (evttable_hdu1["PI"] < (args.max * 100)))[0]
		#append number of high energy photons
                nphotons_ordered.append(len(high_energy_indicies))
		#find and append the exposure
		exposure_ordered.append(evttable_hdu1.meta['EXPOSURE'])

	#Simple photons/second conversion
	normalized_photons = []

	for i in range(len(nphotons_ordered)):
		if nphotons_ordered[i] == 0 or exposure_ordered[i] == 0:
			normalized_photons.append(0)
		else:
			normalized_photons.append(float(nphotons_ordered[i]) / exposure_ordered[i])

#Use hour or minute bins
if args.hbin or args.mbin:
	#These will be our x and y values when we go to make the plot. 
	timebins_main = []
	nphotons_main = []
	for filename in filename_ordered:
		cleanfilt = filename +"/cleanfilt.evt"
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		totalexposure = evttable_hdu1.meta['EXPOSURE']
		exp_per_hour = totalexposure / 3600.0
		exp_per_minute = totalexposure / 60.0
		high_energy_indicies = np.where((evttable_hdu1["PI"] > (args.min * 100)) & (evttable_hdu1["PI"] < (args.max * 100)))[0]
		#Put all values into a directory to sort by timebin
		d1 = {}
		#Do the time correction procedure
		firsttime = datetime(year=2014, month=1, day=1, hour=0, minute=0, second=0)
		for idx in high_energy_indicies:
			time_ref_seconds = evttable_hdu1["TIME"][idx]
			corrected_time = firsttime + timedelta(seconds=time_ref_seconds)
			photon_energy = evttable_hdu1["PI"][idx]
			d1[corrected_time] = photon_energy
		#Make a second dictionary where we will sort the dates by hour or minute
		d2 = {}
		for datekey in d1.keys():
			if args.hbin:
				if datekey.replace(minute=0,second=0, microsecond=0) in d2.keys():
					d2[datekey.replace(minute=0,second=0, microsecond=0)]+=1
				else:
					d2[datekey.replace(minute=0,second=0,microsecond=0)]=1
			if args.mbin:
				if datekey.replace(second=0, microsecond=0) in d2.keys():
					d2[datekey.replace(second=0, microsecond=0)]+=1
				else:
					d2[datekey.replace(second=0,microsecond=0)]=1

		#Now we have hour mins and nphotons in each bin for the filename
		#Append to main lists
		for datekey in d2.keys():
			timebins_main.append(datekey)
			#Divide nphotons by exp_per_hour
			if exp_per_hour == 0:
				nphotons_main.append(0)
			else:
				nphotons_main.append((d2[datekey]))


#for i in range(len(soyuz_dates)): 
#	plt.axvline(x=soyuz_dates[i])


#Make colored regions
for tup in soyuz_dates_rassvet:
	plt.axvspan(tup[0], tup[1], alpha=0.3, color='red')

for tup in soyuz_dates_poisk:
	plt.axvspan(tup[0], tup[1], alpha=0.3, color='blue')

#Set date ranges for plot
if args.dmin is not None:
	dmin = pr.parse(args.dmin)
else:
	dmin = pr.parse('2017-05-01')

dminastropy = Time(dmin, format='datetime', scale='utc')
dminastropy.format = 'mjd'
axisdatemin = dminastropy.value
print("Minimum date range is ", axisdatemin)

if args.dmax is not None:
	dmax = pr.parse(args.dmax)
else:
	dmax = pr.parse('2018-03-01')

dmaxastropy = Time(dmax, format='datetime', scale='utc')
dmaxastropy.format = 'mjd'
axisdatemax = dmaxastropy.value
print("Maximum date range is ", axisdatemax)
#Set axes for plot
if (not args.hbin) and (not args.mbin):
	majorvals = []
	iterdate = axisdatemin
	while True:
		if iterdate > axisdatemax:
			break
		else:
			majorvals.append(iterdate)
			iterdate = iterdate + 5

	plt.axis([axisdatemin, axisdatemax, 0, (1.2*max(normalized_photons))])
	labels = [ str(val) for val in majorvals ]
	plt.xticks(majorvals, labels)

if args.hbin or args.mbin:
	majorvals = []
	iterdate = axisdatemin
	if args.hbin:
		iterval = 1
	if args.mbin:
		iterval = .5
	while True:
		if iterdate > axisdatemax:
			break
		else:
			majorvals.append(iterdate)
			iterdate = iterdate + iterval


	plt.axis([axisdatemin, axisdatemax, 0, 1.2*max(nphotons_main)])
	labels = [ str(val) for val in majorvals ]
	plt.xticks(majorvals, labels)


#Create a bar plot of the results
if (not args.hbin) and (not args.mbin):
	#Conversion method into MJD
	for idx in range(len(dates_ordered)):
		dastropy = Time(dates_ordered[idx], format='datetime', scale='utc')
		dastropy.format = 'mjd'
		dates_ordered[idx] = dastropy.value

	#Width based on using markers every five days
	plt.bar(dates_ordered, normalized_photons, width=1, align='edge')

if args.hbin or args.mbin:
	#Conversion method into MJD
	for idx in range(len(timebins_main)):
		dastropy = Time(timebins_main[idx], format='datetime', scale='utc')
		dastropy.format = 'mjd'
		timebins_main[idx] = dastropy.value
	
	if args.hbin:		
		#Width based on using markers every day
		plt.bar(timebins_main, nphotons_main, color='black', width=(1.0/24.0), align='edge')

	if args.mbin:
		#Width based on using markers every 12 hours
		print(len(timebins_main))
		plt.bar(timebins_main, nphotons_main, color='black', width=(1.0/720.0), align='edge')
plt.xlabel('Date')
plt.ylabel('Photons per second between ' + str(args.min) + 'keV and ' + str(args.max)+ 'keV')


plt.title('1821-24 Photons between ' + str(args.min) + 'keV and ' + str(args.max) + 'keV')
if args.hbin:
	bintype='hbin'
elif args.mbin:
	bintype='mbin'
else:
	bintype='ObsIDbin'
plt.savefig('1821-24_photons_'+str(args.min)+'_'+str(args.max)+'_'+bintype+'.png')
plt.show()
