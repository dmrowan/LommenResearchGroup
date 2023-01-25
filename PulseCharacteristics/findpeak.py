#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
import os

"""
Function used by ampread.py; finds the location and width of both the main pulse and the interpulse in the full profile for one ObsID using a single and double Gaussian curve fit and convolutions
"""

def gauss(x, a, m, s, d): #same as in ampread.py
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d): #same as in functions.py, see gauss2 there
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def findpeak(phase, showplot = False, savedata = False): #finds the location and width of the main pulse and interpulse for one ObsID

    # Gaussian curve fit to only the main pulse
    yvals, xlims = np.histogram(phase,bins=255) # finds heights and sides of each bin, no plot
    xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be the x values of curve 
    
    # Use convolution to find the estimate for the location (pulse phase) of the center of the main pulse (same as in ampread.py)
    width=0.05
    x = xvals
    template = np.exp(-((x)/width)**2) # template for convolution: very basic Gaussian
    convo = []
    for i in range(len(yvals)):
        convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc = xvals[convo.index(m)]  # finds the location (pulse phase) of the peak of convolution; this is the center of the main pulse
    popt, pcov = curve_fit(gauss, xvals, yvals, p0 = [max(yvals), maxloc, 0.05, min(yvals)]) #Gaussian curve fit, with guesses

    if showplot == True: #can use plot to check if fit looks right
        plt.hist(phase, bins=255)

    # Temporarily remove the main pulse in order to find the location (pulse phase) of the interpulse
    phase2 = []
    for i in phase: #initially append all same values from phase (want a deep copy, since we'll modify phase2 but want phase to be the same)
        phase2.append(i)
    for i in range(len(phase2)): #remove the main pulse from phase2 by first setting those values to zero
        if ((phase2[i] > popt[1]-popt[2]) & (phase2[i] < popt[1]+popt[2])):
            phase2[i] = 0
    phase2 = [x for x in phase2 if x!= 0] #remove all zeros
    yvals2, xlims2 = np.histogram(phase2,bins=255) # finds heights and sides of each bin, no plot
    xvals2 = xlims2[:-1] + np.diff(xlims2)/2 # finds middle of each bin, to be x values of surve
    
    # Use convolution to find the estimate for the location of the center of the interpulse
    width=0.05
    x = xvals2
    template = np.exp(-((x)/width)**2) # template for convolution
    convo = []
    for i in range(len(yvals2)):
        convo.append(np.sum(yvals2*np.roll(template,i))) # finds convolution
    m = np.max(convo) # finds peak value of convolution
    maxloc2 = xvals2[convo.index(m)]  # finds the location of the peak of convolution; this is the center of the interpulse
    
    #Use the FULL profile curve (phase, not phase2) to curve fit to BOTH peaks simulataneously
    popt, pcov = curve_fit(gauss2, xvals, yvals, p0 = [max(yvals), maxloc, 0.05, max(yvals2), maxloc2, 0.05, min(yvals)])

    #Make sure main pulse and interpulse are assigned correctly; if not, swap
    #Note we only care about locations and widths here, not max/min values so we don't worry about bridge emission here)
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

    if savedata == True: #can save fit data to a file if needed
        f = open("pulsefitdata.txt", "a")
        print(name, 'main pulse, interpulse: ', peakloc, intloc, file=f)
        f.close()

    if showplot == True: #can plot the final curve fit to both pulses
        plt.plot(xvals, gauss2(xvals, *popt))
        plt.show()

    #Returns the location and width of the main pulse (peakloc, standdev) and interpulse (intloc, instanddev)
    return(peakloc, standdev, intloc, intstanddev)
