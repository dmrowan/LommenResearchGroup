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
    print(tab)


    # User input: choose time per profile, how many plots to show, bins
    timewidth = int(input("How much time in each pulse profile? (in seconds)"))

    # Splits pulse phase data into profiles based on time intervals
    ranges = np.arange(tab['TIME'][0], tab['TIME'][-1], timewidth) # makes a numpy array that splits the     total time into intervals based on the time interval given
    sections = len(ranges)
    print("The number of time intervals is", sections)
    row = int(input("How many rows of subplots?"))
    col = int(input("How many columns of subplots?"))

    log.info("Limit data")
    phases = []
    for starttime in ranges: # goes through each time interval 
        rows = np.where((tab['TIME'] > starttime) & (tab['TIME'] <= starttime + timewidth)) # finds rows in each time interval
        newphase = list(tab['PULSE_PHASE'][rows[0]]) # makes a list of phases in each time interval
        phases.append(newphase)
 
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
