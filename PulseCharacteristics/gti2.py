#!/usr/bin/env python
from astropy.table import Table
from astropy import log
import numpy as np
import pandas as pd

desc="""
Creates bins of "good times" of an arbitrary length N (as determined by the user)
"""
def main(): 

    fname = '1937_events.evt'
    log.info('Read in text file')
    timetab = Table.read('1937_events.evt', hdu=2) 
    tab = Table.read('1937_events.evt', hdu=1)
    timetab = timetab.to_pandas()
    timewidth = int(input("What time interval?"))
    phases = []
    totaltime = 0
    timetab['timedif']= timetab['STOP'] - timetab['START']
    totaltime=0
    startime = timetab['START']
    for counter in range(np.size(timetab[0])):
        totaltime += timetab['timedif'][counter]
        if (totaltime >= timewidth):
            diff = totaltime-timewidth
            totaltime -= diff
            endtime = timetab['STOP'][counter] - diff
            # make profiles
            goodtimes = np.where(time>startime && time < endtime)
            make the profile
            fit it
            write it to a file
            startime=endtime
            totaltime=0

    print(timetab)
    print(counter)
    print(totaltime)

if __name__ == '__main__':
    main()
