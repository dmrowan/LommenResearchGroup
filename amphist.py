#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

desc="""
Plots the convolution of the peak, then does a gaussian curve fit, plots it and returns peak amplitude
"""

def gauss(x, a, b, c, d):
    return(a*np.exp(-((x-b)/c)**2)+d)

def main():    

    # Read in data and made a table
    fname = '1937_events.evt'
    tab = Table.read('1937_events.evt', hdu=1)
 
    #User input
    ntype = input("Do you want to use a fraction of data or seconds of data? (valid input: fraction OR seconds)")
    totalt = tab['TIME'][-1] - tab['TIME'][0]
    if (ntype == "fraction"):
        N = int(input("What fraction of data do you want? (nth of data)"))
        print("Each profile is", totalt/N , "seconds of the data.")
    if (ntype == "seconds"):
        s = int(input("How many seconds of data do you want?"))
        N = totalt/s
  #  plot = input("Show the curve fit plot? (y/n)")

    limitedtab = np.array_split(tab['PULSE_PHASE'], N)
    amplitudes = []

    for n in range(N):
        phase = np.array(limitedtab[n])
        yvals, xlims = np.histogram(phase,bins=255)
        xvals = xlims[:-1] + np.diff(xlims)/2

        #Use convolution to find the estimate for the location of the peak
        width=0.05
        x = xvals
        template = np.exp(-((x)/width)**2)
        convo = []
        for i in range(len(yvals)):
            convo.append(np.sum(yvals*np.roll(template,i)))

        m = np.max(convo)
        maxloc = xvals[convo.index(m)]  # finds the location of the peak

        #Does a gaussian curve fit to the histogram
        popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)],bounds=([min(yvals),0,0,min(yvals)-100],[max(yvals)+100, 0.2, 0.2, max(yvals)]))
       # plt.plot(xvals, gauss(xvals, *popt))

        #Amplitude of fitted curve
        amp = popt[0]
        amplitudes.append(amp)
   
    print(amplitudes)
    plt.show(plt.hist(amplitudes, bins = int(N/2)))
        
if __name__ == '__main__':
    main()

