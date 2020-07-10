#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

desc="""
Plots all histograms and curve fits for given time intervals
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def main():    

    # Reads in data and makes a table
    fname = 'PSR_B1821-24_combined.evt'
    log.info('Read in table')
    tab = Table.read('PSR_B1821-24_combined.evt', hdu=1)
    timetab = Table.read('PSR_B1821-24_combined.evt', hdu=2)
    
    # User input: choose time per profile
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

    # Splits pulse phase data into profiles based on time intervals
    log.info("Limit data")
    phases = []
    starttimes =[]
    endtimes = []
    totaltime = 0
    starttime = timetab['START'][0]
    starttimes.append(starttime)
    for i in range(len(timetab)):
        if (i==(len(timetab) -1)): break
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
        phase = list(phase)
        phases.append(phase)

    sections = len(phases)
    print("The number of time intervals is", sections)
  

   # phases = [x for x in phases if x != []] # gets rid of all empty lists in phases


    # Makes a list of xvalues and fits for each profile

    log.info("Make list of curve fits")
    xvalues = []
    curves = []
    for n in range(len(phases)):
 
        # Makes a line plot from the histogram
        phase = np.array(phases[n])  # uses pulse phases in nth profile
        for i in range(len(phase)):
            if (phase[i] < 0.5):
                phase[i] = 0
        phase = [x for x in phase if x!= 0]
        yvals, xlims = np.histogram(phase,bins=128) # finds heights and sides of each bin, no plot
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
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram 
        if ((popt[1] >= 0.7)&(popt[1]<= 0.85)):
            xvalues.append(xvals)
            curves.append(popt)
        if (popt[1] < 0.7):
            phases[n] = []
        if (popt[1] > 0.85):
            phases[n] = []

    phases = [x for x in phases if x != []]

    notempty = len(phases)
    print("The number of time intervals with data is", notempty)
    row = int(input("How many rows of subplots?"))
    col = int(input("How many columns of subplots?"))

    log.info("Plot all curve fits")
    fig, ax = plt.subplots(row, col, sharex = 'col')
    i = 0
    j = 0
    for n in range(len(phases)):
        if (j > (col-1)):
            j = 0
            i += 1
        ax[i, j].hist(np.array(phases[n]), bins = 255)
        ax[i, j].plot(xvalues[n], gauss(xvalues[n], *curves[n]))
        j += 1 
    
    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical') 
    
    plt.show()
        
if __name__ == '__main__':
    main()

