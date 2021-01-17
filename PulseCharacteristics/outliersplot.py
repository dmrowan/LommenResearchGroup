#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
import csv

desc="""
Plots outliers for the crab
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def outliers(timewidth):
    intint = pd.read_csv('crabintdata_%s.txt'%timewidth, header = None)
    intint = list(intint[0])
    upper = []
    lower = []
    for i in intint:
         if (i > 1.5):
              upper.append(i)
         if (i < 0.8):
              lower.append(i)
    print('The total number of data points is', len(intint))
    print('The number of outliers above integrated intensity of 1.5 is', len(upper))
    print('The number of outliers below integrated intensity of 0.8 is', len(lower))

    fnames = pd.read_csv('crabfilenames.txt', header = None)
    fnames = list(fnames[0])

    filenames =  [fnames[0], fnames[1], fnames[2], fnames[3], fnames[4], fnames[5]]

    upper = []
    upperplot = []
    lower = [] 
    lowerplot = []    
    upperintintensity = []
    lowerintintensity = []

    for name in filenames:
        log.info('Starting new int time')
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

     #   plt.hist(phase, bins=255)

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
       
        # Splits pulse phase data into profiles based on time intervals
        phases = []
        starttimes =[]
        endtimes = []
        totaltime = 0
        starttime = timetab['START'][0]
        starttimes.append(starttime)
        for i in range(4):
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
            if (a<0):
                for i in range(len(phase)):
                    if ((phase[i] >= 0) & (phase[i] < b)):
                        phase[i] = 0
                    if (phase[i] > (1+a)) & (phase[i] <=1):
                        phase[i] = 0
                binnumber = int(200-(200*(b-a)))
            if (a>=0):
                for i in range(len(phase)):
                    if ((phase[i] > a) & (phase[i] < b)):
                        phase[i] = 0
                binnumber = 200
            phase = [x for x in phase if x!= 0]
            yvals, xlims = np.histogram(phase,bins=binnumber) # finds heights and sides of each bin, no plot
            xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot
            if (a>=0):
                for i in range(len(xvals)):
                    if ((xvals[i]>=a) & (xvals[i]<=b)):
                        xvals[i] = 999
                        yvals[i] = 999
            if (a<0):
                for i in range(len(xvals)):
                    if ((xvals[i] >= 0) & (xvals[i] < b)):
                        xvals[i] = 999
                        yvals[i] = 999
                    if ((xvals[i] > (1+a)) & (xvals[i] <= 1)):
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
            try:
                for i in range(len(yvals)):
                    convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
            except ValueError:
                removed.append(n)
                continue
            m = np.max(convo) # finds peak value of convolution
            maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution

            # Does a gaussian curve fit to the histogram
            try:
                popt2, pcov2 = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf)))
               # plt.hist(phases[n], bins=200)
            except RuntimeError:
                removed.append(n)
                continue

            intint = (popt2[0]*popt2[2]*np.sqrt(2*np.pi))/timewidth
            if ((popt2[1] >= peakloc-(standdev*4)) & (popt2[1] <= peakloc+(standdev*4))):
                if ((intint < 1.5) & (intint >0.8)):
                    with open('centerintint_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow([intint])
                    with open('centercurve_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow(list(popt2))
                    with open('centerhist_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow((phases[n]))
                """ if (intint > 1.5):
                    with open('upperintint_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow([intint])
                    with open('uppercurve_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow(list(popt2))
                    with open('upperhist_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow((phases[n]))
                    
                if (intint < 0.8):
                    with open('lowerintint_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow([intint])
                    with open('lowercurve_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow(list(popt2))
                    with open('lowerhist_%s.csv'%timewidth, 'a+', newline= '') as file:
                       output = csv.writer(file, delimiter='\t')
                       output.writerow((phases[n])) """
            else:
                removed.append(n)
 
time = [10, 30, 60, 90, 120]
for t in time:
    outliers(t)

