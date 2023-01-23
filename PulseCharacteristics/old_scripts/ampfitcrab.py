#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from functions import *

desc="""
Makes a histogram of the integrated intensities of pulse profiles for each integration time
"""

def fit(timewidth):
 
    intint = pd.read_csv('intdata/crabintdata_%s.txt' %timewidth, header = None)
    intint = list(intint[0])
    intint = [x for x in intint if x > 0.001]
    intint = [x for x in intint if x < 3]
    intint = np.array(intint)

    print('The total number of profiles in the %s pulse histogram is '%timewidth, len(intint))

    binwidths = list(np.linspace(0, 0.12, 100))

    plt.hist(intint, binwidths)

    width = 0.05
    sd = np.std(intint)  # calculates standard deviation directly
    
    # Makes a line plot from the histogrm
    intint = np.array(intint)  # uses pulse phases in nth profile
    xvals, yvals = hist_to_curve(intint, binwidths)

    # Use convolution to find the estimate for the location of the peak
    #x = xvals
    #template = (np.exp(-((x/width)**2)/2))
    #convo = []
    #for i in range(len(yvals)):
    #    convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
    #m = np.max(convo) # finds peak value of convolution
    #maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution

    #popt, pcov = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc, width], bounds = ((0, 0, 0), (np.inf, np.inf, np.inf))) # uses gaussian function to do a curve fit to the line version fo the histogram; uses maxloc for the guess for location
   
    #popt, pcov = curve_fit(lognormal, xvals, yvals, p0 = [100, 0, 2])
    #print(popt)
    #plt.plot(xvals, lognormal(xvals, *popt))
    
    plt.hist(intint, bins=binwidths)
    #plt.plot(xvals, gauss(xvals,*popt))
    #errorbar = np.absolute(pcov[2][2])**0.5

    plt.xlabel('Integrated intensity (counts/pulse)')
    plt.ylabel('# of Profiles')
    plt.title('Integrated intensity distribution for %s pulses/profile'%timewidth)
    plt.savefig('crab_%s.png' % timewidth)
    plt.clf()
    return()

plottype = 'loglog'
width = []
width2 = []
errorbars = []
timewidth=[]
times = [15, 20, 30, 90, 150, 300, 900]
for twidth in times: 
    if (twidth == 0):
        twidth = 10
    fit(twidth)
    #width.append(w)
    #width2.append(w2)
    #errorbars.append(e)
    #timewidth.append(twidth)
#y =[]
"""
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
    plt.title("Width of Integrated Intensity Distribution vs Integration Time")
    plt.legend()
    plt.xlabel("Integration Time (seconds)")
    plt.ylabel("Standard Deviation (counts/second)")


if (plottype == 'loglog'):
    plt.plot(np.log10(timewidth), np.log10(width), 'o', color = 'b')
    #popt, pcov = curve_fit(power, timewidth, width)
    fit, cov =  np.polyfit(np.log10(timewidth), np.log10(width), 1, cov=True)
    cslope = fit[0]
    csloperror = np.absolute(cov[0][0])**0.5
    #cslopeerror = np.absolute(pcov[1][1])**0.5
    plt.plot(np.log10(timewidth), np.log10(timewidth)*fit[0]+fit[1], color = 'b', label = 'Standard deviation, %s$\pm$%s'%(float('%.2g' % cslope), float('%.1g' % csloperror)))
    plt.plot(np.log10(timewidth), np.log10(width2), 'o', color = 'g')
    
    #popt, pcov = curve_fit(power, timewidth, width2)
    fit, cov =  np.polyfit(np.log10(timewidth), np.log10(width2), 1, cov=True)
    fslope = fit[0]
    #fslopeerror = np.absolute(pcov[1][1])**0.5
    fsloperror = np.absolute(cov[0][0])**0.5
    plt.plot(np.log10(timewidth), np.log10(timewidth)*fit[0]+fit[1], color = 'g', label = 'Gaussian curve fit, %s$\pm$%s'%(float('%.2g' % fslope), float('%.1g' % fsloperror)))
    
    shift = 1.3
    plt.plot(np.log10(timewidth), np.log10(y)-shift, '--', label = 'Gaussian distribution (slope = -1/2)')
    
    lowererror = []
    uppererror = []
    for x in range(len(errorbars)):
        lowererror.append(abs(np.log10(width2[x]-errorbars[x])-np.log10(width2[x])))
        uppererror.append(np.log10(width2[x]+errorbars[x])-np.log10(width2[x]))
    loglogerrors =  [lowererror, uppererror]
    loglogerror = np.array(loglogerrors)
    plt.errorbar(np.log10(timewidth),np.log10(width2), yerr = loglogerror , fmt = 'o', color = 'g')
    plt.title("Width of Integrated Intensity Distribution vs Number of Pulses")
    plt.legend()
    plt.xlabel("log(Number of Pulses per Profile)")
    plt.ylabel("log(Width (counts/second))")
    
    print(cslope, csloperror)
    print(fslope, fsloperror)
plt.show()

"""
