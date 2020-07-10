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
    intint = []
    number = len(phases)
    for n in range(number):
 
        # Makes a line plot from the histogram
        phase = np.array(phases[n])  # uses pulse phases in nth profile
        binnumber = 255
        if (pulsarname == '1821'):
            for i in range(len(phase)):
                if (phase[i] < 0.5):
                    phase[i] = 0
            phase = [x for x in phase if x!= 0]
            bins = 128
        if (pulsarname == 'crab'):
            for i in range(len(phase)):
                if (phase[i] < 0.5):
                    phase[i] = 0
            phase = [x for x in phase if x!= 0]
            bins = 128
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

        # Amplitude and integrated intensity of fitted curve
        if (pulsarname == '1937'):
            amp = (popt[0]*popt[2]*np.sqrt(2*np.pi))/timewidth 
            if (popt[1] <= 0.2):
                intint.append(amp) # appends amplitude of peak into list of amplitudes
            else:
                removed.append(n)        
        if (pulsarname == '1821'):
            amp = (popt[0]*popt[2]*np.sqrt(2*np.pi))/timewidth
            if ((popt[1] >= 0.7) & (popt[1] <= 0.85)):
                intint.append(amp) # appends amplitude of peak into list of amplitudes
            else:
                removed.append(n)
        if (pulsarname == 'crab'):
            amp = (popt[0]*popt[2]*np.sqrt(2*np.pi))/timewidth
            if (popt[1] >= 0.8):
                intint.append(amp)
            else:
                removed.append(n)
   
   
    binwidths = list(np.arange(0, 0.0002, 0.000005))
    width = 0.00005
    if (pulsarname == 'crab'):
        binwidths = list(np.arange(0.4, 1.2, 0.01))
        width = 0.1
    plt.hist(intint, bins = binwidths) # makes histogram of amplitudes
  
    sd = np.std(intint)  # calculates standard deviation directly

    # Makes a line plot from the histogram
    intint = np.array(intint)  # uses pulse phases in nth profile
    yvals, xlims = np.histogram(intint,bins=binwidths) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot

    # Use convolution to find the estimate for the location of the peak
    x = xvals
    template = (np.exp(-((x/width)**2)/2))
    convo = []
    for i in range(len(yvals)):
        convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution
 
    popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, width, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
    plt.plot(xvals, gauss(xvals,*popt))
    errorbar = np.absolute(pcov[2][2])**0.5

    plt.xlabel('Integrated Intensity')
    plt.ylabel('Counts')
    plt.title('Integrated Intensities of Pulse Profiles')
    f = open("amphistdata.txt", "a")
    print("standard deviation1 = ", sd, file=f)
    print("standard deviation2 = ", popt[2], file=f)
    f.close()
    if (pulsarname == '1937'):
        plt.savefig('1937_%s.png' % timewidth)
    if (pulsarname == '1821'):
        plt.savefig('1821_%s.png' % timewidth)
    if (pulsarname == 'crab'):
        plt.savefig('crab_%s.png' % timewidth)
    plt.clf()
    return(sd, popt[2], errorbar)

width = []
width2 = []
errorbars = []
timewidth=[]
for twidth in range(1800, 9000, 900):
    if (twidth == 1800):
        twidth = 2100
    w, w2, e = fit('1937', twidth)
    width.append(w)
    width2.append(w2)
    errorbars.append(e)
    timewidth.append(twidth)
y =[]
for i in timewidth:
    y.append(1/(i**0.5))
plottype = 'loglog'
if (plottype == 'plot'):
    plt.plot(timewidth, width, 'o', color = 'b')
    popt, pcov = curve_fit(power, timewidth, width)
    calc = plt.plot(timewidth, power(timewidth, *popt), color = 'b', label = 'Calculated')
 #   plt.plot(timewidth, width2, 'o', color = 'g')
    popt, pcov = curve_fit(power, timewidth, width2)
    fit = plt.plot(timewidth, power(timewidth, *popt), color = 'g', label = 'Fitted')
    plt.errorbar(timewidth, width2, yerr = errorbars, fmt = 'o', color = 'g')
    plt.title("Widths of Integrated Intensity Distribution vs Integration Time")
    plt.legend()
    plt.xlabel("Integration Time (seconds)")
    plt.ylabel("Standard Deviation (counts/second)")
if (plottype == 'loglog'):
    plt.plot(np.log10(timewidth), np.log10(width), 'o', color = 'b')
    popt, pcov = curve_fit(power, timewidth, width)
    cslope = popt[1]
    cslopeerror = np.absolute(pcov[1][1])**0.5
    plt.plot(np.log10(timewidth), np.log10(power(timewidth, *popt)), color = 'b', label = 'Calculated')
  #  plt.plot(np.log10(timewidth), np.log10(width2), 'o', color = 'g')
    popt, pcov = curve_fit(power, timewidth, width2)
    fslope = popt[1]
    fslopeerror = np.absolute(pcov[1][1])**0.5
    plt.plot(np.log10(timewidth), np.log10(power(timewidth, *popt)), color = 'g', label = 'Fitted')
    plt.plot(np.log10(timewidth), np.log10(y)-2.7, '--', label = 'Constant')
    lowererror = []
    uppererror = []
    for x in range(len(errorbars)):
        lowererror.append(abs(np.log10(width2[x]-errorbars[x])-np.log10(width2[x])))
        uppererror.append(np.log10(width2[x]+errorbars[x])-np.log10(width2[x]))
    loglogerrors =  [lowererror, uppererror]
    loglogerror = np.array(loglogerrors)
    plt.errorbar(np.log10(timewidth),np.log10(width2), yerr = loglogerror , fmt = 'o', color = 'g')
    plt.title("Widths of Integrated Intensity Distribution vs Integration Time")
    plt.legend()
    plt.xlabel("log(Integration Time (seconds))")
    plt.ylabel("log(Standard Deviation (counts/second))")
    print(cslope, cslopeerror)
    print(fslope, fslopeerror)
plt.show()


