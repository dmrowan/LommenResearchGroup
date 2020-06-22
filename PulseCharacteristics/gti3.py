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
    starttimes =[]
    endtimes = []
    totaltime = 0
    starttime = timetab['START'][0]
    starttimes.append(starttime)
    for i in range(len(timetab)):
        if (i==4503): break
        interval = timetab['STOP'][i] - starttime
        if (timewidth < interval):
            number = int(interval/timewidth)
            for n in range(number):
                endtime = starttime + timewidth
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
            difference = timetab['STOP'][i] - starttime
            gaptime = timetab['START'][i+1] - timetab['STOP'][i]
            endtime = starttime + timewidth + gaptime
            endtimes.append(endtime)
            starttime = endtime
            starttimes.append(starttime)
        if (timewidth == interval):
            endtime = starttime + timewidth
            endtimes.append(endtime)
            starttime = endtime
            starttimes.append(starttime)
        if (timewidth > interval):
            totaltime += interval
            if (totaltime >= timewidth):
                diff = totaltime - timewidth
                endtime = timetab['STOP'][i] - diff
                endtimes.append(endtime)
                starttime = endtime
                starttimes.append(starttime)
                totaltime = diff
            starttime = timetab['START'][i+1]

    for i in range(len(starttimes)-1):
        rows = np.where((tab['TIME'] >= starttimes[i]) & (tab['TIME'] <= endtimes[i]))
        phase = tab['PULSE_PHASE'][rows[0]]
        phases.append(list(phase))


    print(len(phases))
    print(starttimes[0:20])
    print(endtimes[0:20])

if __name__ == '__main__':
    main()
