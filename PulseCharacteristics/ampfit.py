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


def fit(pulsarname, timewidth):

    # Reads in data and makes a table
    if (pulsarname == '1937'):
        fname = '1937_events.evt'
    if (pulsarname == '1821'):
        fname = 'PSR_B1821-24_combined.evt'
    log.info('Read in table')
    tab = Table.read(fname, hdu=1)
    timetab = Table.read(fname, hdu=2)
 
    # User input: choose time per profile
 #   timewidth = int(input("How much time in each pulse profile? (in seconds)"))


    # Splits pulse phase data into profiles based on time intervals
    log.info("Limit data")
    phases = []
    starttimes =[]
    endtimes = []
    totaltime = 0
    starttime = timetab['START'][0]
    starttimes.append(starttime)
    for i in range(len(timetab)):
        if (i==4503): break
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
  #  print("The number of time intervals is", totalphases)
  #  plotnumber = input("Use all profiles or first 84? (all/84)")

    # Makes a list of amplitude of peak in each profile
    log.info("Make list of amplitudes")
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
   
   
    log.info("Amplitude histogram")
  #  print("The number of profiles removed due to insufficient data is", len(removed)+emptyremoved)
    binwidths = list(np.arange(0,0.015 , 0.00001))
    plt.hist(amplitudes, bins = binwidths) # makes histogram of amplitudes
  
    # Makes a line plot from the histogram
    amplitudes = np.array(amplitudes)  # uses pulse phases in nth profile
    yvals, xlims = np.histogram(amplitudes,bins=binwidths) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot

    # Use convolution to find the estimate for the location of the peak
    width=0.005
    x = xvals
    template = np.exp(-((x)/width)**2) # template for convolution
    convo = []
    for i in range(len(yvals)):
        convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution
        
    popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, 0.005, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
    width = 2*np.sqrt(2*(math.log(2)))*(popt[2])
 #   print("The width is", width)
    plt.plot(xvals, gauss(xvals,*popt))
    plt.xlabel('Amplitude of Peak')
    plt.ylabel('Counts')
    plt.title('Amplitudes of Pulse Profiles')
  #  plt.legend(['width =', width], ['amplitude =', popt[0]], ['center =', popt[1]])
    f = open("amphistdata.txt", "a")
    print("width = ", width, file=f)
    print("amplitude = ", popt[0], file=f)
    print("center = ", popt[1], file=f)
    f.close()
    plt.savefig('%s.png' % timewidth)
    plt.clf()
    return(popt[0], width)

amp = []
width = []
timewidth=[]
for twidth in range(1800, 7200, 1800):
    a, w = fit('1937', twidth)
    amp.append(a)
    width.append(w)
    timewidth.append(twidth)
plt.plot(timewidth, width)
