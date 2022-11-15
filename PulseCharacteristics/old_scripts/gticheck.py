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
    
    filenames =  [fnames[0]]
   
    for name in filenames:
        tab = Table.read(name, hdu=1) # reads in pulse phase data
        timetab = Table.read(name, hdu=2) # reads in gti data
        
        phase = np.array(tab['PULSE_PHASE'])  # uses pulse phases in nth profile
 
        # Splits pulse phase data into profiles based on time intervals
        phases = []
        starttimes =[]
        endtimes = []
        totaltime = 0
        totalt = []
        starttime = timetab['START'][0]
        starttimes.append(starttime)
       # for i in range(len(timetab)):
        for i in range(10):
           # if (i==(len(timetab)-1)): break
            if (i==9): break
            string = ':'
            interval = timetab['STOP'][i] - starttime #amount of time in the interval
            if interval == 0:
                print('error2')
            if (timewidth < interval): # >1 intergration time fits in one interval
                string = string + '1'
                number = int(interval/timewidth) 
                for n in range(number):
                    endtime = starttime + timewidth
                    endtimes.append(endtime)
                    string =  string + 'a'
                    totaltime = endtime - starttime
                    starttime = endtime
                    starttimes.append(starttime)
                    totalt.append(totaltime)
                    totaltime = 0
                difference = timetab['STOP'][i] - starttime
                gaptime = timetab['START'][i+1] - timetab['STOP'][i]
                if (timetab['START'][i+1] + difference)<=timetab['STOP'][i+1]:
                     endtime = starttime + timewidth + gaptime
                else:
                     print('error1')
                endtimes.append(endtime)
                string = string + 'A'
                totaltime = endtime - starttime
                starttime = endtime
                starttimes.append(starttime)
                totalt.append(totaltime)
                totaltime = 0
            if (timewidth == interval): #fits exactly into interval
                string = string + '2b'
                endtime = starttime + timewidth
                endtimes.append(endtime)
                totaltime = endtime - starttime
                starttime = endtime
                starttimes.append(starttime)
                totalt.append(totaltime)
                totaltime = 0
            if (timewidth > interval): #<1 intergration time fits into interval
                string = string + '3'
                totaltime += interval
                if (totaltime >= timewidth):
                    diff = totaltime - timewidth
                    endtime = timetab['STOP'][i] - diff
                    endtimes.append(endtime)
                    string = string + 'c'
                    starttime = endtime
                    starttimes.append(starttime)
                    totalt.append(totaltime)
                    totaltime = diff
                starttime = timetab['START'][i+1]
            f = open("gtitest_%s.txt" % timewidth, "a")
            print(string)
            f.close()
        print(starttimes)
        print(endtimes)
        print(totalt)
        for i in range(len(starttimes)-1):
            rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))

        totalphases = len(phases)
        phases = [x for x in phases if x != []]
        sections = len(phases)
        print(sections)
        print(len(phases[0]))
        countrate = []
        for x in range(len(phases)):
            countrate.append(len(phases[x])/timewidth)
        ranges = np.arange(sections)
        print(ranges)
        plt.plot(ranges, countrate, '.')
        plt.ylabel('Count Rate/Second')
        plt.xlabel('Time (intervals of %s second(s))'%timewidth)
        plt.show()

times = [10]
for time in times:
    integrationtimes(time)
