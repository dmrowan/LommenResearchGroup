#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


desc="""
Plots full data set in one plot for Crab
"""

def main():
    fnames = np.loadtxt('validobsids.txt', str)
    fnames = list(fnames)

    fname =  fnames[0]
    path = '/homes/alevina/research2020/PSR_B0531+21/%s_pipe/cleanfilt2.evt'%fname
    
    log.info('Read in table')
    tab = Table.read(path, hdu=1)

    log.info("Making Pulse Profile")
    plt.hist(tab['PULSE_PHASE'], bins=255)
    plt.title("Pulse Profile of ObsID %s"%fname)
    plt.xlabel("Pulse Phase")
    plt.ylabel("Counts")
    plt.show()

if __name__ == '__main__':
    main()
