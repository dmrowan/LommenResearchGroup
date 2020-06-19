#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt


#Modified from code by Dom Rowan, 2020

desc="""
Plots full data set in one plot for 1937
"""

def main():
    fname = '1937_events.evt'

    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)
    print(tab)

    log.info("Making Pulse Profile")
    fig, ax = plt.subplots(1, 1, figsize=(12,6))
    ax.set_xlabel("Pulse Phase", fontsize=20)
    ax.set_ylabel("Counts", fontsize=20)
    ax.hist(tab['PULSE_PHASE'], bins=255, color='xkcd:darkblue')

    plt.show()

if __name__ == '__main__':
    main()
