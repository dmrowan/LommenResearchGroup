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

    fname = 'PSR_B1821-24_combined.evt'
    filenames =  [fname]
   
    for name in filenames:
        tab = Table.read(name, hdu=1)
        timetab = Table.read(name, hdu=2)

        phase = np.array(tab['PULSE_PHASE'])  # uses pulse phases in nth profile
        yvals, xlims = np.histogram(phase,bins=255) # finds heights and sides of each bin, no plot
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
        popt, pcov = curve_fit(gauss, xvals, yvals, p0 = [max(yvals), maxloc, 0.05, min(yvals)])
    
        for i in range(len(phase)):
            if ((phase[i] > popt[1]-popt[2]) & (phase[i] < popt[1]+popt[2])):
                phase[i] = 0
        phase = [x for x in phase if x!= 0]
        yvals2, xlims2 = np.histogram(phase,bins=255) # finds heights and sides of each bin, no plot
        xvals2 = xlims2[:-1] + np.diff(xlims2)/2 # finds middle of each bin, to be x values of line plot
        # Use convolution to find the estimate for the location of the peak
        width=0.05
        x = xvals2
        template = np.exp(-((x)/width)**2) # template for convolution
        convo = []
        for i in range(len(yvals2)):
            convo.append(np.sum(yvals2*np.roll(template,i))) # finds convolution
        m = np.max(convo) # finds peak value of convolution
        maxloc2 = xvals2[convo.index(m)]  # finds the location of the peak of convolution

        popt, pcov = curve_fit(gauss2, xvals, yvals, p0 = [max(yvals), maxloc, 0.05, max(yvals2), maxloc2, 0.05, min(yvals)])
        if (popt[0] > popt[3]):
            peakloc = popt[1]
            standdev = popt[2]
            intloc = popt[4]
            intstanddev = popt[5]
        if (popt[0] < popt[3]):
            peakloc = popt[4]
            standdev = popt[5]
            intloc = popt[1]
            intstanddev = popt[2]
        print(peakloc, standdev, intloc, intstanddev)

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
            a = intloc-(intstanddev*2)
            b = intloc+(intstanddev*2)
            for i in range(len(phase)):
                if ((phase[i] > a) & (phase[i] < b)):
                    phase[i] = 0
            phase = [x for x in phase if x!= 0]
           # binnumber = int(200-(200*(b-a)))
            binnumber = 200
            yvals, xlims = np.histogram(phase,bins=binnumber) # finds heights and sides of each bin, no plot
            xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot
            for i in range(len(xvals)):
                if ((xvals[i]>=a) & (xvals[i]<=b)):
                    xvals[i] = 999
                    yvals[i] = 999
            xvals = [x for x in xvals if x!= 999]
            yvals = [x for x in yvals if x!= 999]
            xvals = np.array(xvals)
            yvals = np.array(yvals)
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
                popt3, pcov3 = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) 
            except RuntimeError:
                removed.append(n)
                continue

            intint = (popt3[0]*popt3[2]*np.sqrt(2*np.pi))/timewidth
            f = open("1821intdata_%s.txt" % timewidth, "a")
            if ((popt3[1] >= peakloc-(standdev*4)) & (popt3[1] <= peakloc+(standdev*4))):
                print(intint, file=f)
            else:
                removed.append(n)
            f.close()
        print(timewidth, len(phases), len(removed))

for time in range(1800, 9000, 900):
    if (time == 1800):
        time = 2100
    integrationtimes(time)
