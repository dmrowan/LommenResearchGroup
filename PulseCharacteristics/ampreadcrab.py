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

def power(x, a, b):
    return(a*(x**(-b)))

def integrationtimes(timewidth):

    fname = '/students/pipeline/heasoft6.27/PSR_B0531+21/1013010121_pipe/cleanfilt.evt'
    fname2 = '/students/pipeline/heasoft6.27/PSR_B0531+21/1013010127_pipe/cleanfilt.evt'
    filenames =  [fname, fname2]
   
    for name in filenames:
        tab = Table.read(name, hdu=1)
        timetab = Table.read(name, hdu=2)
 
        # Splits pulse phase data into profiles based on time intervals
        phases = []
        starttimes =[]
        endtimes = []
        totaltime = 0
        starttime = timetab['START'][0]
        starttimes.append(starttime)
        for i in range(len(timetab)):
            if (i==(len(timetab)-1)): break
            interval = timetab['STOP'][i] - starttime
            if (timewidth < interval):
                number = int(interval/timewidth)
                for n in range(number):
                    endtime = starttime + timewidth
                    endtimes.append(endtime)
                    starttime = endtime
                    starttimes.append(starttime)
                difference = timetab['STOP'][i] - starttime
                gaptime = timetab['START'][i+1] - timetab['STOP'][i]
                endtime = starttime + timewidth + gaptime
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
            if (timewidth == interval):
                endtime = starttime + timewidth
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
            if (timewidth > interval):
                totaltime += interval
                if (totaltime >= timewidth):
                    diff = totaltime - timewidth
                    endtime = timetab['STOP'][i] - diff
                    endtimes.append(endtime)
                    starttime = endtime
                    starttimes.append(starttime)
                    totaltime = diff
                starttime = timetab['START'][i+1]

        for i in range(len(starttimes)-1):
            rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))

        totalphases = len(phases)
        phases = [x for x in phases if x != []]
        sections = len(phases)
        emptyremoved = totalphases - sections

        # Makes a list of amplitude of peak in each profile
        removed = []
        number = len(phases)
        for n in range(number):
            # Makes a line plot from the histogram
            phase = np.array(phases[n])  # uses pulse phases in nth profile
            for i in range(len(phase)):
                if (phase[i] < 0.5):
                    phase[i] = 0
            phase = [x for x in phase if x!= 0]
            binnumber = 128
            yvals, xlims = np.histogram(phase,bins=binnumber) # finds heights and sides of each bin, no plot
            xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot

            # Use convolution to find the estimate for the location of the peak
            width=0.05
            x = xvals
            template = np.exp(-((x)/width)**2) # template for convolution
            convo = []
            for i in range(len(yvals)):
                convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
            m = np.max(convo) # finds peak value of convolution
            maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution

            # Does a gaussian curve fit to the histogram
            try:
                popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
            except RuntimeError:
                removed.append(n)
                continue

            intint = (popt[0]*popt[2]*np.sqrt(2*np.pi))/timewidth
            if (popt[1] >= 0.8):
                f = open("crabintdata_%s.txt" % timewidth, "a")
                print(intint, file=f)
                f.close()
            else:
                removed.append(n)

times = [30, 120, 150, 180]
for time in times:
    integrationtimes(time)
