#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os.path
import os
import pandas as pd

desc="""
Plots all data of each ObsID of the Crab
"""

def main():

    fnames = pd.read_csv('/homes/alevina/research2020/LommenResearchGroup/PulseCharacteristics/validobsids.txt', header = None)
    fnames = list(fnames[0])
    filenames =  fnames

    phaselist = []
    for f in filenames:
        path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe'% f
        path = path + '/cleanfilt2.evt'

        log.info('Read in table for ObsID %s'%f)
        tab = Table.read(path, hdu=1)
        phases = np.array(tab['PULSE_PHASE'])
   
        phaselist.append(phases)
        del(phases)
    
    log.info("Making Pulse Profile")
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
        ax[i, j].hist(phaselist[n], bins = 255)
        j += 1

    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')

    plt.show()

if __name__ == '__main__':
    main()
