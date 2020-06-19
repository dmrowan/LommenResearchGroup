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
    counter = 0
    starttime = timetab['START'][0] 
    timetab['timedif'] = timetab['STOP'] - timetab['START'] 
    while ((totaltime + timetab['timedif'][counter]) < timewidth):
        totaltime += timetab['timedif'][counter]
        counter += 1    
    if (totaltime < timewidth):
        counter += 1
        diff = timewidth - totaltime
        totaltime += diff
        endtime = timetab['STOP'][counter] + diff
        rows = np.where((tab['TIME'] >= starttime) & (tab['TIME'] <= endtime))
        phase = tab['PULSE_PHASE'][rows[0]]
        phases.append(list(phase))
        starttime = endtime
        totaltime = 0
    print(timetab)
    print(counter)
    print(phases)
  
    """
    ranges = []
    rowcount = 0
    rangecount = 0
    shift = 0
    while rowcount <= len(timetab):
        for i in timetab['timedif']:
            while (totaltime + i  < (timewidth*rangecount)):
                totaltime += i
                counter += 1
        if totaltime < timewidth:
            counter += 1
            diff = timewidth - totaltime
            totaltime += diff
            ranges.append(totaltime)
            shift = timetab['timedif'][counter] - diff
            totaltime = shift
        rowcount += counter
        rangecount +=1 
      
    print(counter)
    print(totaltime)
    print(ranges)
    print(rowcount)"""
#   log.info("Limit data")
#   ranges = np.arange(tab[0], tab[length], tab[length])
         

if __name__ == '__main__':
    main()
