#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

desc="""
Reads in table and plots pulse profiles (histograms) based on time in each profile 
"""

def main():

    # Reads in data and makes a table
    fname = '/students/pipeline/heasoft6.27/PSR_B0531+21/1013010121_pipe/cleanfilt.evt'  
    log.info('Read in table')
    tab = Table.read(fname, hdu=1)
    timetab = Table.read(fname, hdu=2)
    timetab = timetab.to_pandas()
    print(tab)


    # User input: choose time per profile, how many plots to show, bins
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

    log.info("Limit data")
    phases = []
    starttimes =[]
    endtimes = []
    totaltime = 0
    starttime = timetab['START'][0]
    starttimes.append(starttime)
    for i in range(len(timetab)):
        if (i==(len(timetab)-1)): break
        interval = timetab['STOP'][i] - starttime
        if (timewidth < interval):
            number = int(interval/timewidth)
            for n in range(number):
                endtime = starttime + timewidth
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
            difference = timetab['STOP'][i] - starttime
            gaptime = timetab['START'][i+1] - timetab['STOP'][i]
            endtime = starttime + timewidth + gaptime
            endtimes.append(endtime)
            starttime = endtime
            starttimes.append(starttime)
        if (timewidth == interval):
            endtime = starttime + timewidth
            endtimes.append(endtime)
            starttime = endtime
            starttimes.append(starttime)
        if (timewidth > interval):
            totaltime += interval
            if (totaltime >= timewidth):
                diff = totaltime - timewidth
                endtime = timetab['STOP'][i] - diff
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
                totaltime = diff
            starttime = timetab['START'][i+1]
    
    for i in range(len(starttimes)-1):
        rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
        phase = tab['PULSE_PHASE'][rows[0]]
        phases.append(list(phase))
  
    sections = len(phases)
    print("The number of time intervals is", sections)
    plotnumber = input("Plot all profiles or first 84? (all/84)")
    if (plotnumber == 'all'):
        row = int(input("How many rows of subplots?"))
        col = int(input("How many columns of subplots?"))
    if (plotnumber == '84'):
        row = 12
        col = 7

    # Plots histograms of the profiles 
    log.info("Making Pulse Profile")
   
    fig, ax = plt.subplots(row, col, sharex = 'col')
    i = 0
    j = 0
    if (plotnumber == 'all'):
        for n in range(len(phases)):
            if (j > (col-1)):
                j = 0
                i += 1
            ax[i, j].hist(phases[n], bins = 255)
            j += 1    
  
    if (plotnumber == '84'):  
        for n in range(84):
            if (j > (col-1)):
                j = 0
                i += 1
            ax[i, j].hist(phases[n], bins = 255)
            j += 1
  
  
    for axs in ax.flat:
        axs.label_outer() # hides x labels and ticks for top plots and y ticks for plots on right

    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')
    

    plt.show() 


if __name__ == '__main__':
    main()
