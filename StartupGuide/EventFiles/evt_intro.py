#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt


#Dom Rowan, 2020

desc="""
Introduction to EVT files for NICER data analysis
"""

def convert_time(time):
    timezero = datetime.datetime(year=2014, month=1,
                                 day=1, hour=0, minute=0, second=0)
    new_time = timezero+datetime.timedelta(seconds=time)

    return new_time


def main():
    fname = '1937_events.evt'

    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)
    print(tab)

    log.info("Listing columns")
    print(tab.colnames)

    log.info("Making Pulse Profile")
    fig, ax = plt.subplots(1, 1, figsize=(12,6))
    ax.set_xlabel("Pulse Phase", fontsize=20)
    ax.set_ylabel("Counts", fontsize=20)
    ax.hist(tab['PULSE_PHASE'], bins=255, color='xkcd:darkblue')

    plt.show()

    log.info("Example of time conversion")
    example_time = 142039018.640167803
    fixed_time = convert_time(example_time)
    print(f'The time {round(example_time)} is equivalent to {fixed_time}')

if __name__ == '__main__':
    main()
