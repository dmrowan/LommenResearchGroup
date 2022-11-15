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
    p = int(input("How many plots do you want?"))
    b = int(input("How many bins?"))
    random = input("Plot random plots? (y/n)")
 

    # Splits pulse phase data into profiles based on time intervals
    log.info("Limit data")
    ranges = np.arange(tab['TIME'][0], tab['TIME'][-1], timewidth) # makes a numpy array that splits the     total time into intervals based on the time interval given

    phases = []
    for starttime in ranges: # goes through each time interval 
        rows = np.where((tab['TIME'] > starttime) & (tab['TIME'] <= starttime + timewidth)) # finds rows in each time interval
        newphase = list(tab['PULSE_PHASE'][rows[0]]) # makes a list of phases in each time interval
        phases.append(newphase)
    phases = [x for x in phases if x != []]
  
    # Plots histograms of the profiles 
    log.info("Making Pulse Profile")
    fig, axs = plt.subplots(p,1, figsize=(8,6)) 
    fig.suptitle('Histograms')
   
    if (random == 'n'): 
        for i in range(p): # based on number of plots from input
            axs[i].hist(phases[i], bins=b, color='xkcd:darkblue') # plots p histograms
 
    if (random == 'y'):
        numbers = np.random.randint(low=1, high=len(phases), size=p)
        c = 0
        for i in numbers:
            axs[c].hist(phases[i], bins=b, color='xkcd:darkblue') # plots p random histograms
            print(i, c)
            c += 1

    for axs in ax.flat:
        axs.set(xlabel='Pulse Phase', ylabel='Counts') # labels plots

    for axs in ax.flat:
        axs.label_outer() # hides x labels and ticks for top plots and y ticks for plots on right

    plt.show() 


if __name__ == '__main__':
    main()
