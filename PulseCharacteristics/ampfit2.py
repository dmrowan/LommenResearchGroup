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

def gauss(x, a, m, s):
    return(a*np.exp(-(((x - m)/s)**2)/2))

def power(x, a, b):
    return(a*(x**(-b)))

def fit(pulsarname, timewidth):
 
    if (pulsarname == '1937'):
        intint = pd.read_csv('1937intdata_%s.txt' %timewidth, header = None)
    if (pulsarname == '1821'):
        intint = pd.read_csv('1821intdata_%s.txt' %timewidth, header = None)
    intint = list(intint[0])

    if (pulsarname == '1821'): 
        binwidths = list(np.arange(0, 0.0004, 0.000005))
    if (pulsarname == '1937'):
        binwidths = list(np.arange(0, 0.0002, 0.000005))
    width = 0.00005
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
 
    popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, width], bounds = ((0, 0, 0), (np.inf, np.inf, np.inf))) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
    plt.plot(xvals, gauss(xvals,*popt))
    errorbar = np.absolute(pcov[2][2])**0.5

    plt.xlabel('Integrated Intensity')
    plt.ylabel('Counts')
    plt.title('Integrated Intensities of Pulse Profiles')
    if (pulsarname == '1937'):
        plt.savefig('1937_%s.png' % timewidth)
    if (pulsarname == '1821'):
        plt.savefig('1821_%s.png' % timewidth)
    plt.clf()
    return(sd, popt[2], errorbar)

pname = '1821'
plottype = 'plot'
width = []
width2 = []
errorbars = []
timewidth=[]
for twidth in range(1800, 9000, 900): 
    if (twidth == 1800):
        twidth = 2100
    w, w2, e = fit(pname, twidth)
    width.append(w)
    width2.append(w2)
    errorbars.append(e)
    timewidth.append(twidth)
y =[]
for i in timewidth:
    y.append(1/(i**0.5))
if (plottype == 'plot'):
    plt.plot(timewidth, width, 'o', color = 'b')
    popt, pcov = curve_fit(power, timewidth, width)
    calc = plt.plot(timewidth, power(timewidth, *popt), color = 'b', label = 'Standard deviation')
 #   plt.plot(timewidth, width2, 'o', color = 'g')
    popt, pcov = curve_fit(power, timewidth, width2)
    fit = plt.plot(timewidth, power(timewidth, *popt), color = 'g', label = 'Gaussian curve fit')
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
    plt.plot(np.log10(timewidth), np.log10(power(timewidth, *popt)), color = 'b', label = 'Standard deviation')
  #  plt.plot(np.log10(timewidth), np.log10(width2), 'o', color = 'g')
    popt, pcov = curve_fit(power, timewidth, width2)
    fslope = popt[1]
    fslopeerror = np.absolute(pcov[1][1])**0.5
    plt.plot(np.log10(timewidth), np.log10(power(timewidth, *popt)), color = 'g', label = 'Gaussian curve fit')
    if (pname == '1937'):
        shift = 2.9
    if (pname == '1821'):
        shift = 2.6
    plt.plot(np.log10(timewidth), np.log10(y)-shift, '--', label = 'Gaussian distibution (slope = -1/2)')
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


