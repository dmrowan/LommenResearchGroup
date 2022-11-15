#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os.path
import os

desc="""
Plots full data of ONE ObsID for the Crab
"""

def main():

    psr = str(input('ObsID: '))
    path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe/cleanfilt2.evt'% psr
    tab = Table.read(path, hdu=1)
    phases = np.array(tab['PULSE_PHASE'])
   
    print("The total number of photons is", len(phases))
    
    log.info("Making Pulse Profile")
    fig, ax = plt.subplots(1, 1, figsize=(12,6))
    ax.set_xlabel("Pulse Phase", fontsize=20)
    ax.set_ylabel("Counts", fontsize=20)
    ax.hist(phases, bins=255, color='xkcd:darkblue')

    plt.show()

if __name__ == '__main__':
    main()
