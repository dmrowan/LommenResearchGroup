#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

desc="""
Does a gaussian curve fit to the pulse profile in time interval, plots it, and returns peak amplitude
"""

def gauss(x, a, b, c, d):  # defines gaussian function to use for curve fit
    return(a*np.exp(-((x-b)/c)**2)+d)

def main():

    # Reads in data and makes a table
    fname = '1937_events.evt' # for pulsar 1973
    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)


    # User input: choose time per profile and if to show plot
    plot = input("Show the curve fit plot? (y/n)")
    histogram = input("Show the histogram? (y/n)")
    interpulse = input("Show the interpulse curve fit? (y/n)")
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))
    b = int(input("How many bins?"))


    # Splits pulse phase data into profiles based on time intervals
    log.info("Limit data")
    ranges = np.arange(tab['TIME'][0], tab['TIME'][-1], timewidth) # makes a numpy array that splits the     total time into intervals based on the time interval given
    sections = len(ranges)
    print("The number of time intervals is", sections)
    n = int(input("Which time interval do you want?")) # choose which profile to plot

    phases = []
    for starttime in ranges: # goes through each time interval
        rows = np.where((tab['TIME'] > starttime) & (tab['TIME'] <= starttime + timewidth)) # finds rows in each time interval
        newphase = list(tab['PULSE_PHASE'][rows[0]]) # makes a list of phases in each time interval
        phases.append(newphase)
    phases = [x for x in phases if x != []]


    # Makes a line plot from the histogram 
    phase = np.array(phases[n])  # uses pulse phases in nth profile
    if (histogram == 'y'):
        yvals, xlims, _ = plt.hist(phase,bins=b) # finds heights and sides of each bin and plots histogram
    if (histogram == 'n'):
        yvals, xlims = np.histogram(phase,bins=b) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of line plot
    #  plt.plot(xvals, yvals)   # option to plot the pulse profile as a line


    # Use convolution to find the estimate for the location of the peak
    log.info("Convolution")
    width=0.05
    x = xvals
    template = np.exp(-((x)/width)**2) # template for convolution
    convo = []
    for i in range(len(yvals)):  
        convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution 
    #  plt.plot(xvals, convo)   # option to plot the convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc = xvals[convo.index(m)]  # finds the location of the peak of convoltion


    # Does a gaussian curve fit to the histogram
    log.info("Curve fit")
    popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)]) # uses gaussian    function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location 
    plt.plot(xvals, gauss(xvals, *popt)) # plots the fitted curve
    if (interpulse == "y"):
        popt2, pcov2 = curve_fit(gauss, xvals, yvals, bounds = ([min(yvals)-100, 0.5, 0.025, min(yvals)-100], [max(yvals), 0.7, 0.1, max(yvals)]))
        plt.plot(xvals, gauss(xvals, *popt2), color = 'orange') # plots the fitted curve to interpulse
    plt.xlabel('Pulse Phase')
    plt.ylabel('Counts')
    plt.title('Gaussian Fit to Pulse Profile')
 

    # Amplitude of fitted curve
    amp = popt[0] # finds amplitude of fitted curve (the first parameter of curve fit)
    print(amp)

    if (plot == "y"): # based on user input, plots fitted curve
        plt.show()

if __name__ == '__main__':
    main()

