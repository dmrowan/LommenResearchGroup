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

def gauss(x, a, b, c, d):  # defines gaussian function to use for curve fit
    return(a*np.exp(-((x-b)/c)**2)+d)

def main():    

    # Reads in data and makes a table
    fname = '1937_events.evt'
    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)
 
    # User input: choose time per profile
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

    # Splits pulse phase data into profiles based on time intervals
    ranges = np.arange(tab['TIME'][0], tab['TIME'][-1], timewidth) # makes a numpy array that splits the     total time into intervals based on the time interval given
    sections = len(ranges)
    print("The number of time intervals is", sections)

    log.info("Limit data")
    phases = []
    for starttime in ranges: # goes through each time interval
        rows = np.where((tab['TIME'] > starttime) & (tab['TIME'] <= starttime + timewidth)) # finds rows in each time interval
        newphase = list(tab['PULSE_PHASE'][rows[0]]) # makes a list of phases in each time interval
        phases.append(newphase)
    phases = [x for x in phases if x != []] # gets rid of all empty lists in phases


    # Makes a list of xvalues and fits for each profile

    log.info("Make list of curve fits")
    xvalues = []
    curves = []
    for n in range(len(phases)):
 
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
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram 
        if (popt[1] <= 0.2):
            xvalues.append(xvals)
            curves.append(popt)           
        if (popt[1] > 0.2):
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

