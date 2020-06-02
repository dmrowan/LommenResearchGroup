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

#Dom Rowan, 2018
#This code needs to be updated to fix axes and bin widths (see SoyuzPlotMJD for how this is done)

desc = """
This looks at the number of ratio rejected events as a function of time from bkf prefilt files. Ensure that your directory only contains the files 
for the selected pulsar. This version has the option to select a time range and use hour or minute bins in addition to obs  date bins
"""

parser = argparse.ArgumentParser(description = desc)
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

#Create an empty directory that will have our time and filename. 
files_unordered = {}

#Loop through the files in the current directory. This doesn't go sequentially. 
for filename in os.listdir(os.getcwd()):
	if os.path.isfile(filename+"/"+filename[:10]+"_prefilt.bkf"):
		#First set the path to the evt file
                bkf = filename + "/"+filename[:10]+"_prefilt.bkf"
		
		bkftable_hdu1 = Table.read(bkf, hdu=1)
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
	nbadratio_ordered = []
	exposure_ordered = []
	#For each file in order
	for filename in filename_ordered:
		bkf = filename + "/"+filename[:10]+"_prefilt.bkf"
		bkftable_hdu1 = Table.read(bkf, hdu=1)
		#Number of rejected ratio events
		nbadratio = sum(bkftable_hdu1['BAD_RATIO'])
		#append number of RRE
                nbadratio_ordered.append(nbadratio)
		#find and append the exposure
		cleanfilt = filename + "/cleanfilt.evt"
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		exposure_ordered.append(evttable_hdu1.meta['EXPOSURE'])

	#Simple rre/second conversion
	normalized_RRE = []

	for i in range(len(nbadratio_ordered)):
		if nbadratio_ordered[i] == 0 or exposure_ordered[i] == 0:
			normalized_RRE.append(0)
		else:
			normalized_RRE.append(float(nbadratio_ordered[i])/ exposure_ordered[i])

#Use hour or minute bins
if args.hbin or args.mbin:
	#These will be our x and y values when we go to make the plot. 
	timebins_main = []
	nrre_main = []
	for filename in filename_ordered:
		cleanfilt = filename +"/cleanfilt.evt"
		evttable_hdu1 = Table.read(cleanfilt, hdu=1)
		totalexposure = evttable_hdu1.meta['EXPOSURE']
		exp_per_hour = totalexposure / 3600.0
		exp_per_minute = totalexposure / 60.0
		
		bkf = filename + "/"+filename[:10]+"_prefilt.bkf"
		bkftable_hdu1 = Table.read(bkf, hdu=1)
		nbadratio_indicies = np.where(bkftable_hdu1['BAD_RATIO'] > 0)[0]
		#Put all values into a directory to sort by timebin
		d1 = {}
		#Do the time correction procedure
		firsttime = datetime(year=2014, month=1, day=1, hour=0, minute=0, second=0)
		for idx in nbadratio_indicies:
			time_ref_seconds = bkftable_hdu1["TIME"][idx]
			corrected_time = firsttime + timedelta(seconds=time_ref_seconds)
			rre = bkftable_hdu1["BAD_RATIO"][idx]
			d1[corrected_time] = rre
		#Make a second dictionary where we will sort the dates by hour or minute
		d2 = {}
		for datekey in d1.keys():
			if args.hbin:
				if datekey.replace(minute=0,second=0, microsecond=0) in d2.keys():
					d2[datekey.replace(minute=0,second=0, microsecond=0)] += d1[datekey]
				else:
					d2[datekey.replace(minute=0,second=0,microsecond=0)]= d1[datekey]
			if args.mbin:
				if datekey.replace(second=0, microsecond=0) in d2.keys():
					d2[datekey.replace(second=0, microsecond=0)]+= d1[datekey]
				else:
					d2[datekey.replace(second=0,microsecond=0)]= d1[datekey]

		#Now we have hour mins and n rejected events in each bin for the filename
		#Append to main lists
		for datekey in d2.keys():
			timebins_main.append(datekey)
			#Divide nrre by exp_per_hour
			if exp_per_hour == 0:
				nrre_main.append(0)
			else:
				nrre_main.append((d2[datekey])/exp_per_hour)

#Create a bar plot of the results
if (not args.hbin) and (not args.mbin):
	plt.bar(dates_ordered, normalized_RRE, fill=False, width=0.1)

if args.hbin or args.mbin:
	plt.bar(timebins_main, nrre_main, fill=False, width=0.1)

plt.xlabel('Date')
plt.ylabel('Ratio Rejected Events / S')

#for i in range(len(soyuz_dates)): 
#	plt.axvline(x=soyuz_dates[i])


#Make colored regions
for tup in soyuz_dates_rassvet:
	plt.axvspan(tup[0], tup[1], alpha=0.3, color='red')

for tup in soyuz_dates_poisk:
	plt.axvspan(tup[0], tup[1], alpha=0.3, color='blue')

#Set date ranges for plot
if args.dmin is not None:
	axisdatemin = pr.parse(args.dmin)
else:
	axisdatemin = pr.parse('2017-05-01')

if args.dmax is not None:
	axisdatemax = pr.parse(args.dmax)
else:
	axisdatemax = pr.parse('2018-03-01')

#Set axes for plot
if (not args.hbin) and (not args.mbin):
	plt.axis([axisdatemin, axisdatemax, 0, (1.2*max(normalized_RRE))])

if args.hbin or args.mbin:
	plt.axis([axisdatemin, axisdatemax, 0, 1.2*max(nrre_main)])

plt.title('1821-24 Ratio Rejected Events')
if args.hbin:
	bintype='hbin'
elif args.mbin:
	bintype='mbin'
else:
	bintype='ObsIDbin'
plt.savefig('1821-24_rre'+bintype+'.png')
plt.show()

