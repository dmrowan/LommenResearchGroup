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

def gauss(x, a, b, c, d):  # defines gaussian function to use for curve fit
    return(a*np.exp(-((x-b)/c)**2)+d)

def power(x, a, b):
    return(a*(x**(-b)))

def fit(pulsarname, timewidth):

    # Reads in data and makes a table
    if (pulsarname == '1937'):
        fname = '1937_events.evt'
    if (pulsarname == '1821'):
        fname = 'PSR_B1821-24_combined.evt'
    if (pulsarname == 'crab'):
        fname = '/students/pipeline/heasoft6.27/PSR_B0531+21/1013010121_pipe/cleanfilt.evt'
    log.info('Starting new timewidth')
    tab = Table.read(fname, hdu=1)
    timetab = Table.read(fname, hdu=2)
 
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
    amplitudes = []
    number = len(phases)
    for n in range(number):
 
        # Makes a line plot from the histogram
        phase = np.array(phases[n])  # uses pulse phases in nth profile
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

         # Does a gaussian curve fit to the histogram
        try:
            popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
        except RuntimeError:
            removed.append(n)
            continue

        # Amplitude of fitted curve
        if (pulsarname == '1937'):
            amp = popt[0]/timewidth  # finds amplitude of fitted curve (the first parameter of curve fit)
            if (popt[1] <= 0.2):
                amplitudes.append(amp) # appends amplitude of peak into list of amplitudes
            else:
                removed.append(n)        
        if (pulsarname == '1821'):
            amp = popt[0]/timewidth
            if ((popt[1] >= 0.7) & (popt[1] <= 0.85)):
                amplitudes.append(amp) # appends amplitude of peak into list of amplitudes
            else:
                removed.append(n)
        if (pulsarname == 'crab'):
            amp = popt[0]/timewidth
            if (popt[1] >= 0.8):
                amplitudes.append(amp)
            else:
                removed.append(n)
   
   
    if (pulsarname == '1937'):
        binwidths = list(np.arange(0,0.015 , 0.00025))
    if (pulsarname == '1821'):
        binwidths = list(np.arange(0,0.015 , 0.00025))
    if (pulsarname == 'crab'):
        binwidths = list(np.arange(16, 20, 0.1))
    plt.hist(amplitudes, bins = binwidths) # makes histogram of amplitudes
  
    # Makes a line plot from the histogram
    amplitudes = np.array(amplitudes)  # uses pulse phases in nth profile
    yvals, xlims = np.histogram(amplitudes,bins=binwidths) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot

    # Use convolution to find the estimate for the location of the peak
    width=0.005
    if (pulsarname == 'crab'):
        width = 1
    x = xvals
    template = np.exp(-((x)/width)**2) # template for convolution
    convo = []
    for i in range(len(yvals)):
        convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution
    
    if (pulsarname == '1937'):    
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, 0.005, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
    if (pulsarname == '1821'):
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, 0.005, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) 
    if (pulsarname == 'crab'):
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, 1, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) 
    width = popt[2] #2*np.sqrt(2*(math.log(2)))*(popt[2])
    plt.plot(xvals, gauss(xvals,*popt))
    errorbar = np.absolute(pcov[2][2])**0.5
    plt.xlabel('Amplitude of Peak')
    plt.ylabel('Counts')
    plt.title('Amplitudes of Pulse Profiles')
    f = open("amphistdata.txt", "a")
    print("width = ", width, file=f)
    print("amplitude = ", popt[0], file=f)
    print("center = ", popt[1], file=f)
    f.close()
    if (pulsarname == '1937'):
        plt.savefig('1937_%s.png' % timewidth)
    if (pulsarname == '1821'):
        plt.savefig('1821_%s.png' % timewidth)
    if (pulsarname == 'crab'):
        plt.savefig('crab_%s.png' % timewidth)
    plt.clf()
    return(popt[0], width, errorbar)

amp = []
width = []
errorbars = []
timewidth=[]
for twidth in range(1800, 9000, 900):
    if (twidth == 1800):
        twidth = 2100
    a, w, e = fit('1937', twidth)
    amp.append(a)
    width.append(w)
    errorbars.append(e)
    timewidth.append(twidth)
#plt.plot(timewidth, width, 'o')
popt, pcov = curve_fit(power, timewidth, width)
plt.plot(timewidth, power(timewidth, *popt))
y =[]
for i in timewidth:
    y.append(1/(i**0.5))
#plt.plot(timewidth, y, '--')
plt.errorbar(timewidth, width, yerr = errorbars, fmt = 'o')
plt.title("Widths of Amplitude Histograms")
plt.xlabel("Integration Time (seconds)")
plt.ylabel("Width (counts/second)")
plt.show()


