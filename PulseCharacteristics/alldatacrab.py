#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os.path
import os
import pandas as pd

"""
Plots all data for each ObsID of the Crab as subplots
"""

def main():

    #Reads ObsID names from event file into a list
    fnames = pd.read_csv('/homes/alevina/research2020/LommenResearchGroup/PulseCharacteristics/validobsids.txt', header = None)
    filenames = list(fnames[0])

    #Reads in data from all ObsIDs into a list of arrays
    phaselist = []
    for f in filenames:
        path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe'%f
        path = path + '/cleanfilt2.evt' #an event file is a FITS file

        log.info('Read in table for ObsID %s'%f)
        tab = Table.read(path, hdu=1) #reads all data from event file into a table
        phases = np.array(tab['PULSE_PHASE']) #array of pulse phases
   
        phaselist.append(phases) #appends array of pulse phases into list
        del(phases) 
    
    #Plots pulse profiles in subplots, uses user input questions to allow for any number ObsIDs to be plotted
    #If the number of ObsIDs cannot evenly be distributed in a grid, overestimate number of rows or columns
    #Or else it will give an error
    log.info("Plotting pulse profiles")
    n = len(filenames)
    row = int(input("There are %s ObsIDs; how many rows of subplots? "%n))
    col = int(input("How many columns? "))

    fig, ax = plt.subplots(row, col, sharex = 'col')
    i = 0
    j = 0
    for n in range(len(phaselist)):
        if (j > (col-1)):
            j = 0
            i += 1
        ax[i, j].hist(phaselist[n], bins = 255) #uses 255 bins in each histogram
        j += 1

    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')

    plt.show()

if __name__ == '__main__':
    main()
