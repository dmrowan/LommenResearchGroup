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
    fname = '1937_events.evt'  # for pulsar 1973
    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)
    timetab = Table.read('1937_events.evt', hdu=2)
    timetab = timetab.to_pandas()
    print(tab)


    # User input: choose time per profile, how many plots to show, bins
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

    log.info("Limit data")
    phases = []
    totaltime = 0
    counter = 0
    diff = 0
    starttime = timetab['START'][0]
    timetab['timedif'] = timetab['STOP'] - timetab['START']
    while (counter <= len(timetab)):
        while ((totaltime + timetab['timedif'][counter]) < timewidth):
            totaltime += (timetab['timedif'][counter] - diff)
            diff = 0
            counter += 1
            if (counter == len(timetab)-1):
                endtime = timetab['STOP'][counter]
                rows = np.where((tab['TIME'] >= starttime) & (tab['TIME'] <= endtime))
                phase = tab['PULSE_PHASE'][rows[0]]
                phases.append(list(phase))
                break
        if (totaltime < timewidth):
            if (counter == len(timetab)-1): break
            counter += 1
            diff = timewidth - totaltime
            totaltime += diff
            endtime = timetab['STOP'][counter] + diff
            rows = np.where((tab['TIME'] >= starttime) & (tab['TIME'] <= endtime))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))
            starttime = endtime
            totaltime = 0
     
    sections = len(phases)
    print("The number of time intervals is", sections)
    row = int(input("How many rows of subplots?"))
    col = int(input("How many columns of subplots?"))


    # Plots histograms of the profiles 
    log.info("Making Pulse Profile")
   
    fig, ax = plt.subplots(row, col, sharex = 'col')
    i = 0
    j = 0
    for n in range(len(phases)):
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
