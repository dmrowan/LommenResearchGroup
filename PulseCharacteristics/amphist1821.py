#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

desc="""
Makes a histogram of the amplitudes of the peaks of pulse profiles based on time interval
"""

def gauss(x, a, b, c, d):  # defines gaussian function to use for curve fit
    return(a*np.exp(-((x-b)/c)**2)+d)

def main():    

    # Reads in data and makes a table
    fname = 'PSR_B1821-24_combined.evt'
    log.info('Read in table')
    tab = Table.read('PSR_B1821-24_combined.evt', hdu=1)
 
    # User input: choose time per profile
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))


    # Splits pulse phase data into profiles based on time intervals
    log.info("Limit data")
    ranges = np.arange(tab['TIME'][0], tab['TIME'][-1], timewidth) # makes a numpy array that splits the     total time into intervals based on the time interval given

    phases = []
    for starttime in ranges: # goes through each time interval
        rows = np.where((tab['TIME'] > starttime) & (tab['TIME'] <= starttime + timewidth)) # finds rows in each time interval
        newphase = list(tab['PULSE_PHASE'][rows[0]]) # makes a list of phases in each time interval
        phases.append(newphase)
    phases = [x for x in phases if x != []] # gets rid of all empty lists in phases

    # Makes a list of amplitude of peak in each profile
    log.info("Make list of amplitudes")
    amplitudes = []
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
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
       # plt.plot(xvals, gauss(xvals, *popt))

        # Amplitude of fitted curve
        amp = popt[0]  # finds amplitude of fitted curve (the first parameter of curve fit)
        if ((popt[1] >= 0.7) & (popt[1] <= 0.85)): 
            amplitudes.append(amp) # appends amplitude of peak into list of amplitudes
   
    log.info("Amplitude histogram")
    plt.hist(amplitudes, bins = int(len(amplitudes)/2)) # makes histogram of amplitudes
    plt.xlabel('Amplitude of Peak')
    plt.ylabel('Counts')
    plt.title('Amplitudes of Pulse Profiles')
    plt.show()
        
if __name__ == '__main__':
    main()

