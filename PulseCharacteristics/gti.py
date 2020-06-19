#!/usr/bin/env python
from astropy.table import Table
from astropy import log
import numpy as np
import pandas as pd

desc="""
Creates bins of "good times" of an arbitrary length (as determined by the user)
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
    diff = 0
    starttime = timetab['START'][0] 
    timetab['timedif'] = timetab['STOP'] - timetab['START'] 
    while (counter <= len(timetab)):
        while ((totaltime + timetab['timedif'][counter]) < timewidth):
            totaltime += (timetab['timedif'][counter] - diff)
            diff = 0
            counter += 1  
            if (counter == len(timetab)-1):
                endtime = timetab['STOP'][counter]
                rows = np.where((tab['TIME'] >= starttime) & (tab['TIME'] <= endtime))
                phase = tab['PULSE_PHASE'][rows[0]]
                phases.append(list(phase))
                break
        if (totaltime < timewidth):
            if (counter == len(timetab)-1): break
            counter += 1
            diff = timewidth - totaltime
            totaltime += diff
            endtime = timetab['STOP'][counter] + diff
            rows = np.where((tab['TIME'] >= starttime) & (tab['TIME'] <= endtime))
            phase = tab['PULSE_PHASE'][rows[0]]
            phases.append(list(phase))
            starttime = endtime
            totaltime = 0

    print(counter)
    print("The number of profiles is", len(phases))  


if __name__ == '__main__':
    main()
