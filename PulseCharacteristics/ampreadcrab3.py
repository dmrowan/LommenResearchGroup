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

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def power(x, a, b):
    return(a*(x**(-b)))

def integrationtimes(timewidth):

    fnames = pd.read_csv('crabfilenames.txt', header = None)
    fnames = list(fnames[0])
    
    filenames =  [fnames[0]]
   
    for name in filenames:
        log.info('Starting new int time')
        tab = Table.read(name, hdu=1)
        timetab = Table.read(name, hdu=2)
        
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
    #    print(peakloc, standdev, intloc, intstanddev)        
    #    plt.plot(xvals, gauss2(xvals, *popt))
    #    plt.show()    
 
        # Splits pulse phase data into profiles based on time intervals
        phases = []
        starttimes = [] #start time of range
        endtimes = [] #end time of range
        totalt = [] #in each range
        starttime = timetab['START'][0]
        starttimes.append(starttime)
        totaltime = 0 #counter
        for i in range(len(timetab)):
            if (i==(len(timetab)-1)): break
            if (starttime >= timetab['START'][i]) & (starttime <= timetab['STOP'][i]): #assumption: starttime in the same interval
                interval = timetab['STOP'][i] - starttime
            if (starttime < timetab['START'][i]):
                gaptime = timetab['START'][i] - timetab['STOP'][i-1]
                interval = timetab['STOP'][i] - starttime - gaptime
                starttime = timetab['START'][i]
            if (starttime > timetab['STOP'][i]):
                continue
            if ((timewidth - totaltime) < interval):
                if totaltime != 0: 
                    headstart = totaltime
                else:
                    headstart = 0
                number = int(interval/timewidth) #number of ranges that will fully fit into one interval
                if ((interval/timewidth) >= 1):
                    for n in range(number): #fills as many ranges as possible
                        endtime = starttime + (timewidth - headstart)
                        endtimes.append(endtime)
                        totaltime += (endtime - starttime) 
                        totalt.append(totaltime)
                        starttime = endtime
                        starttimes.append(starttime)
                        totaltime = 0 #resets totaltime for next range
                        headstart = 0
                    if (starttime + timewidth) < timetab['STOP'][i]:
                        endtime = starttime + (timewidth)
                        endtimes.append(endtime)
                        totaltime += (endtime - starttime) 
                        totalt.append(totaltime)
                        starttime = endtime
                        starttimes.append(starttime)
                        totaltime = 0 #resets totaltime for next range
                else:
                    if (interval > (timewidth - headstart)):
                        endtime = starttime + (timewidth - headstart)
                        endtimes.append(endtime)
                        totaltime += (endtime - starttime) 
                        totalt.append(totaltime)
                        starttime = endtime
                        starttimes.append(starttime)
                        totaltime = 0 #resets totaltime for next range
                        headstart = 0
                difference = timetab['STOP'][i] - starttime
                if (difference != 0):
                    gaptime = timetab['START'][i+1] - timetab['STOP'][i]
                    if (timetab['START'][i+1] + (timewidth-difference))<=timetab['STOP'][i+1]:
                        endtime = starttime + timewidth + gaptime
                        endtimes.append(endtime)
                        totaltime += (endtime - starttime - gaptime)
                        totalt.append(totaltime)
                        totaltime = 0
                        starttime = endtime
                        starttimes.append(starttime)
                    else:
                        try:
                            counter = 1
                            totaltime1 = difference
                            gaps = []
                            while totaltime1 < timewidth:
                                g = timetab['START'][i+counter] - timetab['STOP'][i+counter-1]
                                gaps.append(g)
                                m = (timetab['STOP'][i+counter] - timetab['START'][i+counter])
                                totaltime1 += m
                                counter += 1
                                print(counter)
                            diff = totaltime1 - totaltime
                            endtime = starttime + timewidth
                            for value in gaps:
                                endtime += value
                            endtimes.append(endtime)
                            totaltime += (endtime - starttime)
                            for value in gaps:
                                totaltime = totaltime - value
                            totalt.append(totaltime)
                            totaltime = 0
                            starttime = endtime
                            starttimes.append(starttime)                       
                        except IndexError:
                            break 
                else: print(' ')
                continue 
            if ((timewidth - totaltime) == interval):
                endtime = starttime + timewidth
                endtimes.append(endtime)
                totaltime += (endtime - starttime)
                totalt.append(totaltime)
                starttime = endtime
                starttimes.append(starttime)
                totaltime = 0
                continue
            if ((timewidth - totaltime) > interval): 
                if (totaltime + interval) < timewidth:
                    #need this interval, plus next one
                    totaltime += interval
                    starttime = timetab['START'][i+1]
                continue
        if starttimes[-1] == endtimes[-1]:
            starttimes = starttimes[:-1]
   
        for i in range(len(starttimes)-1):
            rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))

        totalphases = len(phases)
        phases = [x for x in phases if x != []]
        sections = len(phases)
        emptyremoved = totalphases - sections

        # Makes a list of amplitude of peak in each profile
        removed = []
        number = len(phases)
        for n in range(number):
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
                continue
            m = np.max(convo) # finds peak value of convolution
            maxloc = xvals[convo.index(m)]  # finds the location of the peak of convolution

            # Does a gaussian curve fit to the histogram
            try:
                popt2, pcov2 = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) 
               # plt.hist(phases[n], bins=200)
               # plt.hist(phase, bins=binnumber)
               # yval, xlims = np.histogram(phase,bins=binnumber)
               # xval = xlims[:-1] + np.diff(xlims)/2 
               # plt.plot(xval, yval, '.') 
               # plt.plot(xvals, gauss(xvals, *popt2))
               # plt.show()
            except RuntimeError:
                removed.append(n)
                continue

            intint = (popt2[0]*popt2[2]*np.sqrt(2*np.pi))/timewidth
            f = open("crabintdata3_%s.txt" % timewidth, "a")
            if ((popt2[1] >= peakloc-(standdev*4)) & (popt2[1] <= peakloc+(standdev*4))):
                print(intint, file=f)
            else:
                removed.append(n)
            f.close()
        print(timewidth, len(phases), len(removed))

times = [10]
for time in times:
    integrationtimes(time)
    print(time)
