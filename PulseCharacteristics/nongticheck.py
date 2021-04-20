#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math

desc="""
Makes a histogram of the amplitudes of the peaks of the interpulses of pulse profiles based on time interval
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def power(x, a, b):
    return(a*(x**(-b)))

def integrationtimes(timewidth):

    fnames = pd.read_csv('crabfilenames.txt', header = None)
    fnames = list(fnames[0])
    
    filenames =  [fnames[2]]
   
    for name in filenames:
        tab = Table.read(name, hdu=1) # reads in pulse phase data
        timetab = Table.read(name, hdu=2) # reads in gti data
 
        phase = np.array(tab['PULSE_PHASE'])  # uses pulse phases in nth profile
 
        # Splits pulse phase data into profiles based on time intervals
        phases = []
        starttimes =[]
        endtimes = []
        starttime = tab['TIME'][0]
        endtime = tab['TIME'][-1]
        timediff = endtime-starttime
        nranges = int(timediff/timewidth)
        starttimes = [starttime]
        print(nranges)
        print(starttime)
        print(endtime)

        for n in range(nranges):
            endtimes.append(starttimes[n] + timewidth)
            if (endtimes[n]+timewidth) < endtime:
                starttimes.append(endtimes[n])
            else:
                break 
 
        for i in range(len(starttimes)-1):
            rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))

    #    print(starttimes)
    #    print(endtimes)
        totalphases = len(phases)
        phases = [x for x in phases if x != []]
        sections = len(phases)
        print(sections)
        print(len(phases[0]))
        countrate = []
        for x in range(len(phases)):
            countrate.append(len(phases[x])/timewidth)
        ranges = np.arange(sections)
    #    print(ranges)
        plt.plot(ranges, countrate)
        plt.ylabel('Count Rate/Second')
        plt.xlabel('Time (intervals of %s second(s))'%timewidth)
        plt.show()

times = [1]
for time in times:
    integrationtimes(time)
