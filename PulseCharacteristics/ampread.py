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
#from gti import gti
from splitprofiles import splitprofile
from findpeak import findpeak

"""
Reads in all data sequentially by ObsID and calculates the integrated intensity and saves integrated intensity data in directory intdata
Has options to plot each ObsID with the curve fit, or to plot multiple ObsIDs as subplots

NOTE!! Do not append more data into existing output files unless you are certain it is not redundant, because this will affect the integrated intensity histogram. Please change the name of the output file before running
"""

def gauss(x, a, m, s, d): #single peak Gaussian WITH shift (different from gauss in functions.py)
    return((a*np.exp(-(((x - m)/s)**2)/2))+d) #a=amplitude; m=mean; s=standard deviation/width; d=shift along vertical axis

def intensity(n_rotations): #reads in all data, splits into profiles of n_rotations pulses each, calculates and saves integrated intensity data

    # Reads ObsID names from event file into a list
    fnames = pd.read_csv('validobsids.txt', header = None) 
    fnames = list(fnames[0]) 
    filenames = fnames #can select specific ObsIDs by indexing 
    
    saveplot = True #saves one example pulse profile plot for each value of n_rotations 

    for name in filenames: #for each ObsID

        # Only needed if plotting subplots
        #outliers = []
        #histogram = []
        #xvalues = []

        # Read in data
        log.info('Starting ObsID %s'%name)
        path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe/cleanfilt2.evt'%name
        tab = Table.read(path, hdu=1) #contains all data
        phase = np.array(tab['PULSE_PHASE']) #pulse phases for each photon
        time = np.array(tab['TIME']) #time corresponding to each pulse phase value
        timetab = Table.read(path, hdu=2) #table of GTIs
        
        # Find locations and width of main pulse (peakloc, standdev) and interpulse (intloc, instanddev)
        # Used later to improve the curve fitting and make sure that the MP and IP are always found correctly
        try: 
            peakloc, standdev, intloc, intstanddev = findpeak(phase)
        except RuntimeError: #if cannot find the peak because it is too close to phase=0=1 
            shift = 0.5
            for p in range(len(phase)): #shifts all pulse values by 0.5, makes sure none are negative or greater than 1
                phase[p] = phase[p] - shift
                if phase[p] < 0:
                    phase[p] = phase[p] + 1
                if phase[p] > 1:
                    phase[p] = phase[p] - 1 
            peakloc, standdev, intloc, intstanddev = findpeak(phase) #try again
            plt.hist(phase, bins=255) #plots pulse profile to make sure it looks right
 
        # Splits pulse phase data into profiles containing n_rotations pulses each
        phases = splitprofile(n_rotations, timetab, phase, time) #phases is a list of lists of phase values, split into profiles
        phases = [x for x in phases if x != []] #remove any empty profiles, just in case

        # Shifts pulses so that they are always at the same phase (makes the curve fitting easier)
        shift = peakloc - 0.8 #shift main pulse to phase=0.8
        for n in range(len(phases)): #for each profile
            for p in range(len(phases[n])): #for each pulse phase value in each profile, shift and make sure none are negative or >1
                phases[n][p] = phases[n][p] - shift
                if phases[n][p] < 0:
                    phases[n][p] = phases[n][p] + 1
                if phases[n][p] > 1:
                    phases[n][p] = phases[n][p] - 1
        
        # Shifts peakloc and intloc values calculated earlier to match the shifted phase values
        peakloc = peakloc - shift
        intloc = intloc - shift
        if intloc < 0:
            intloc = intloc +1
        if intloc > 1:
            intloc = intloc -1  

        # Converts pulse profile histogram into a curve, fits a Gaussian to the main pulse, calculates integrated intensity for each profile
        removed = [] #the index of the pulse profiles where the curve fit cannot find the main pulse will be added here
        number = len(phases) 
        for n in range(number): #for each pulse profile in that ObsID
            
            # Converts the histogram into a curve
            phase = np.array(phases[n])  #uses pulse phases in nth profile
            
            # As the code fits a single peak Gaussian, need to remove the interpulse AND the bridge emission (between MP and IP)
            # Approximate the left and right boundaries of the interpulse, and the right of main pulse, to remove IP and bridge
            a = intloc-(intstanddev*4) #left boundary of interpulse
            if a < 0:
                a = a + 1
            if a > 1:
                a = a - 1
            b = intloc+(intstanddev*4) #right boundary of interpulse
            if b < 0:
                b = b + 1
            if b > 1:
                b = b - 1
            c = peakloc+(standdev*3) #right boundary of interpulse
            diff = c-b #width of bridge emission PLUS interpulse (width of phase range that must be cut)
            if diff > 0:
                for i in range(len(phase)): #go through each pulse phase and set any in bridge+IP to 999 (to be easily removed later)
                    if (phase[i] < b):
                        phase[i]=999
                    if (phase[i] > c):
                        phase[i]=999
                binnumber = int(200-(200*(1-c+b))) #adjusts the number of bins so that they are distributed the same as without removing bridge+IP
            
            # Convert histogram to a curve
            phase = [x for x in phase if x!= 999] #actually remove the bridge+IP phases
            yvals, xlims = np.histogram(phase,bins=binnumber) # finds heights and sides of each bin, no plot unlike using plt.hist
            xvals = xlims[:-1] + np.diff(xlims)/2 # finds middle of each bin, to be x values of curve
            xvals = np.array(xvals)
            yvals = np.array(yvals)

            # Use convolution to find the estimate for the location of the peak
            width=0.05 #approximate width of the main pulse, better to overestimate
            x = xvals
            template = np.exp(-((x)/width)**2) # template for convolution
            convo = []
            try: 
                for i in range(len(yvals)): #use convolution -> sum of yvals and template as template is shifted along x axis
                    convo.append(np.sum(yvals*np.roll(template,i))) # finds convolution
                m = np.max(convo) # finds peak value of convolution, the corresponding x value is the location of the center of the pulse
            except ValueError: 
                removed.append(n) #if can't find a peak, skip and move onto next profile
                continue
            maxloc = xvals[convo.index(m)]  # finds the location (pulse phase) of the peak of the convolution
            
            # Fits a Gaussian to the main pulse; use of "guesses" and bounds helps the curve fit find the peak correctly
            try:
                popt2, pcov2 = curve_fit(gauss, xvals, yvals, p0= [max(yvals),maxloc,0.05, min(yvals)], bounds = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))) 
                if saveplot == True: #use plotting to check if the curve fit looks right; can save, or show plot here
                    plt.clf() #clear any previous plots
                    plt.hist(phases[n], bins=200) #plots full pulse profile as histogram 
                    #plt.hist(phase, bins=binnumber) #plots pulse profile as histogram WITHOUT bridge+interpulse
                    yval, xlims = np.histogram(phase,bins=binnumber) #use curve fit to plot fitted curve
                    xval = xlims[:-1] + np.diff(xlims)/2 
                    #plt.plot(xval, yval, '.') #plots pulse profile as points (no bridge+interpulse)
                    plt.plot(xvals, gauss(xvals, *popt2)) #plots FITTED CURVE to the pulse profile; make sure this looks ok
                    plt.savefig('profile%s.png'%n_rotations) #saves plot
                    #plt.show() #can show plot too
                    plt.clf()
                    saveplot = False #only shows a plot ONCE for each n_rotations, not for each profile; can modify if needed
            except Exception: 
                removed.append(n) #if curve fit cannot find the peak, skip and move onto next profile
                continue

            #Calculate integrated intensity
            intint = (popt2[0]*popt2[2]*np.sqrt(2*np.pi))/n_rotations # amplitude*width*sqrt(2*pi)/number of pulses per profile
            
            #Append integrated intensity to file 
            #PLEASE MAKE SURE you are saving to the correct file. The code APPENDS data, not overwrites the existing file
            #Please do not create redundant data in the file you will use for further analysis or it WILL affect your results
            #Because otherwise you will append the same profile into the output file multiple times, and that will affect the histogram 
            #If testing things, make a temporary file or comment the whole section, and only append if you're sure that data is not in there already
            #I usually add a 2 when testing ("crabintdata2_%s.txt" in the line below) or something like that
            f = open("intdata/crabintdata_%s.txt" %n_rotations, "a") #opens file; appends, not overwrites
            if ((popt2[1] >= peakloc-(standdev*4)) & (popt2[1] <= peakloc+(standdev*4))): #make sure peak was found correctly
                print(intint, file=f) 
                
                """
                #Can use to find outliers/plot profiles which have a specific integrated intensity value; I used these for finding outliers, hence the name
                if intint <= 0.015:
                    print(name)
                    outliers.append(popt2)
                    histogram.append(phases[n])
                    xvalues.append(xvals)
                """
            else:
                removed.append(n) #if the peak is not found correctly (not at correct phase value based on ALL data in that ObsID), skip, move onto next
            f.close()
        
        #Reset values for next ObsID (otherwise had issues)
        del(phase)
        del(time)
        
        """
        #Can use to plot outliers (or any chosen integrated intensity range) as subplots
        #Can change number of rows/columns as needed
        row = 2
        col = 2
        fig, ax = plt.subplots(row, col, sharex = 'col', figsize = (13, 7))
        i = 0
        j = 0
        for k in range(2):
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
        
#Runs the function intensity on N pulses per profile to calculate and append integrated intensity data into output file 
#15 pulses per profile is the lowest you can go and still have sufficient data per profile for a reliable fit (to calculate integrated intensity from)
rotations = [15] #N pulses  per profile; note there are about 30 pulses/second for the Crab
for N in rotations:
    intensity(N) #run the function intensity for N pulses per profile

