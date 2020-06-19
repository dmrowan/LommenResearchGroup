#!/usr/bin/env python
from astropy.table import Table
from astropy import log
from astropy.io import ascii
import numpy as np
import pandas as pd

desc="""
Creates bins of "good times" of an arbitrary length N (as determined by the user)
"""
def main(): 

    fname = '1937_events.evt'
    log.info('Read in text file')
    timetab = Table.read('1937_events.evt', hdu=2)
    tab = Table.read('1937_events.evt', hdu =1)
 #   tab.remove_column('EVENT_FLAGS')
 #   timetab =  timetab.to_pandas()
 #   tab = tab.to_pandas()
    timewidth = int(input("What time interval?"))
    starttime = []
    stoptime = []
    timedif =  timetab['STOP'][0] - timetab['START'][0]
    starttime.append(timetab['STOP'][0])
    if (timewidth > timedif):
        difference = timewidth - timedif
        nexttime =  timetab['STOP'][1] + difference
        stoptime.append(nexttime)
        print(stoptime)
        for i in range(1, len(timetab)):
            timedif = timetab['STOP'][i] - stoptime[i-1]
            count = i + 1
            difference = timewidth - timedif
            if (difference > 0):
                nexttime = timetab['START'][count] + difference
                if (nexttime <= timetab['STOP'][count]):
                    stoptime.append(nexttime)
                else:
                    difference  = nexttime - timetab['STOP'][count]
                    count += 1
    print(starttime)
    print(stoptime)        
    

if __name__ == '__main__':
    main()
