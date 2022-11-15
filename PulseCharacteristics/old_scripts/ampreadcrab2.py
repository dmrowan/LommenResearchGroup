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
from gti import gti
from splitprofiles import splitprofile
from findpeak import findpeak

desc="""
Makes a histogram of the amplitudes of the peaks of the interpulses of pulse profiles based on time interval
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def power(x, a, b):
    return(a*(x**(-b)))

def intint(n_rotations):

    fnames = pd.read_csv('validobsids.txt', header = None)
    fnames = list(fnames[0])
    filenames =  fnames
    saveplot = True
    
    for name in filenames:

        outliers = []
        histogram = []
        xvalues = []

        log.info('Starting ObsID %s'%name)
        path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe/cleanfilt2.evt'%name
        tab = Table.read(path, hdu=1)
        phase = np.array(tab['PULSE_PHASE'])
        time = np.array(tab['TIME'])
        timetab = Table.read(path, hdu=2)
        
        try: 
            peakloc, standdev, intloc, intstanddev = findpeak(phase)
        except RuntimeError: 
            #plt.hist(phase, bins=255)
            #plt.show()
            shift = 0.5
            for p in range(len(phase)):
                phase[p] = phase[p] - shift
                if phase[p] < 0:
                    phase[p] = phase[p] + 1
                if phase[p] > 1:
                    phase[p] = phase[p] - 1 
            peakloc, standdev, intloc, intstanddev = findpeak(phase)
            plt.hist(phase, bins=255)
 
        # Splits pulse phase data into profiles based on time intervals
        """
        starttimes, endtimes = gti(timewidth, timetab, phase) 
        phases = []
        phase = np.array(phase)
        time = np.array(time)
        print(len(phase), len(time))
        for i in range(len(starttimes)-1):
            rows = np.where((time >= starttimes[i]) & (time <= endtimes[i]))
            phasevals = phase[rows[0]]
            phases.append(list(phasevals))

        totalphases = len(phases)
        phases = [x for x in phases if x != []]
        sections = len(phases)
        emptyremoved = totalphases - sections
        """

        phases = splitprofile(n_rotations, timetab, phase, time)
        phases = [x for x in phases if x != []]

        # Shifts pulses so that they are always at the same phase
        shift = peakloc - 0.8
        for n in range(len(phases)):
            for p in range(len(phases[n])):
                phases[n][p] = phases[n][p] - shift
                if phases[n][p] < 0:
                    phases[n][p] = phases[n][p] + 1
                if phases[n][p] > 1:
                    phases[n][p] = phases[n][p] - 1
        
        peakloc = peakloc - shift
        intloc = intloc - shift
        if intloc < 0:
            intloc = intloc +1
        if intloc > 1:
            intloc = intloc -1 
        #print(peakloc, intloc)

        # Makes a list of amplitude of peak in each profile
        removed = []
        number = len(phases)
        for n in range(number):
            # Makes a line plot from the histogram
            phase = np.array(phases[n])  # uses pulse phases in nth profile
            a = intloc-(intstanddev*4)
            if a < 0:
                a = a + 1
            if a > 1:
                a = a - 1
            b = intloc+(intstanddev*4)
            if b < 0:
                b = b + 1
            if b > 1:
                b = b - 1
            c = peakloc+(standdev*3)
            diff = c-b
            if diff > 0:
                for i in range(len(phase)):
                    if (phase[i] < b):
                        phase[i]=999
                    if (phase[i] > c):
                        phase[i]=999
                binnumber = int(200-(200*(1-c+b)))
            """
            if diff < 0:
                for i in range(len(phase)):
                    if (phase[i] > c)&(phase[i] < b):
                        phase[i]=999
                binnumber = 200
            """
            phase = [x for x in phase if x!= 999]
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
                m = np.max(convo) # finds peak value of convolution
            except ValueError:
                removed.append(n)
                continue
            maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution
            
            m = peakloc 
            s = standdev
            # Does a gaussian curve fit to the histogram
            try:
                popt2, pcov2 = curve_fit(lambda xvals, a, d: gauss(xvals, a, m, s, d), xvals, yvals) 
                #print(popt2)
                if saveplot == True:
                    #print(binnumber)
                    #print(a, b, c)
                    plt.clf()
                    plt.hist(phases[n], bins=200)
                    #plt.hist(phase, bins=binnumber)
                    yval, xlims = np.histogram(phase,bins=binnumber)
                    xval = xlims[:-1] + np.diff(xlims)/2 
                    #plt.plot(xval, yval, '.') 
                    plt.plot(xvals, gauss(xvals, popt2[0], m, s, popt2[1]))
                    plt.savefig('profile%s.png'%n_rotations)
                    #plt.show()
                    plt.clf()
                    saveplot = False
            except RuntimeError:
                removed.append(n)
                continue

            intint = (popt2[0]*s*np.sqrt(2*np.pi))/n_rotations
            f = open("intdata/crabintdata2_%s.txt" %n_rotations, "a")
            print(intint, file=f)
            if intint <= 0.015:
                print(name)
                outliers.append(popt2)
                histogram.append(phases[n])
                xvalues.append(xvals)

                """
                if intint <= 1:
                    print(name)
                    outliers.append(popt2)
                    histogram.append(phases[n])
                    xvalues.append(xvals)
                if intint >= 2:
                    print(name)
                    outliers.append(popt2)
                    histogram.append(phases[n])
                    xvalues.append(xvals)
                """
            else:
                removed.append(n)
            f.close()
        del(phase)
        del(time)
        
        """
        row = 3
        col = 3
        fig, ax = plt.subplots(row, col, sharex = 'col', figsize = (13, 7))
        i = 0
        j = 0
        for k in range(9):
            if (j > (col-1)):
                j = 0
                i += 1
            ax[i, j].hist(histogram[k], bins=200)
            ax[i, j].plot(xvalues[k], gauss(xvalues[k], *outliers[k]))
            j += 1
        fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
        fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')
        plt.show()
        """

rotations = [900, 300, 150, 30, 15] #rotations
for t in rotations:
    intint(t)

