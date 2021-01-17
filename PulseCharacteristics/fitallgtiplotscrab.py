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

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def main():    

    fnames = pd.read_csv('crabfilenames.txt', header = None)
    fnames = list(fnames[0])
    fname =  fnames[1]

    tab = Table.read(fname, hdu=1)
    timetab = Table.read(fname, hdu=2)
    
    # User input: choose time per profile
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

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
    log.info("Limit data")
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

    sections = len(phases)
    print("The number of time intervals is", sections)
  

   # phases = [x for x in phases if x != []] # gets rid of all empty lists in phases


    # Makes a list of xvalues and fits for each profile

    log.info("Make list of curve fits")
    xvalues = []
    curves = []
    intints = []
    removed = []
    for n in range(len(phases)):
 
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
            phases[n] = []
            continue
        m = np.max(convo) # finds peak value of convolution
        maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution

        # Does a gaussian curve fit to the histogram
        popt2, pcov2 = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram 
        intint = (popt2[0]*popt2[2]*np.sqrt(2*np.pi))/timewidth
        if ((popt2[1] >= peakloc-(standdev*4)) & (popt2[1] <= peakloc+(standdev*4))):    
            xvalues.append(xvals)
            curves.append(popt2)         
            intints.append(intint)  
        else:
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
        ax[i, j].hist(np.array(phases[n]), bins = 200)
        ax[i, j].plot(xvalues[n], gauss(xvalues[n], *curves[n]))
        j += 1 
    
    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical') 
    
    plt.show()
        
if __name__ == '__main__':
    main()

