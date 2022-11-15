#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math

def gti(timewidth, timetab, phase):
    phases = []
    starttimes = [] #start time of range
    endtimes = [] #end time of range
    totalt = [] #in each range
    starttime = timetab['START'][0]
    starttimes.append(starttime)
    totaltime = 0 #counter

    for i in range(len(timetab)):
        if (i==(len(timetab)-1)): break
        if (starttime >= timetab['START'][i]) & (starttime <= timetab['STOP'][i]): #assumption: starttime in the same interval
            interval = timetab['STOP'][i] - starttime
        if (starttime < timetab['START'][i]):
            gaptime = timetab['START'][i] - timetab['STOP'][i-1]
            interval = timetab['STOP'][i] - starttime - gaptime
            starttime = timetab['START'][i]
        if (starttime > timetab['STOP'][i]):
            continue
        if ((timewidth - totaltime) < interval):
            if totaltime != 0: 
                headstart = totaltime
            else:
                headstart = 0
            number = int(interval/timewidth) #number of ranges that will fully fit into one interval
            if ((interval/timewidth) >= 1):
                for n in range(number): #fills as many ranges as possible
                    endtime = starttime + (timewidth - headstart)
                    endtimes.append(endtime)
                    totaltime += (endtime - starttime) 
                    totalt.append(totaltime)
                    starttime = endtime
                    starttimes.append(starttime)
                    totaltime = 0 #resets totaltime for next range
                    headstart = 0
                if (starttime + timewidth) < timetab['STOP'][i]:
                    endtime = starttime + (timewidth)
                    endtimes.append(endtime)
                    totaltime += (endtime - starttime) 
                    totalt.append(totaltime)
                    starttime = endtime
                    starttimes.append(starttime)
                    totaltime = 0 #resets totaltime for next range
            else: 
                if interval > (timewidth - headstart):
                    endtime = starttime + (timewidth - headstart)
                    endtimes.append(endtime)
                    totaltime += (endtime - starttime) 
                    totalt.append(totaltime)
                    starttime = endtime
                    starttimes.append(starttime)
                    totaltime = 0 #resets totaltime for next range
                    headstart = 0
            difference = timetab['STOP'][i] - starttime
            if (difference != 0):
                gaptime = timetab['START'][i+1] - timetab['STOP'][i]
                if (timetab['START'][i+1] + (timewidth - difference))<=timetab['STOP'][i+1]:
                    endtime = starttime + timewidth + gaptime
                    endtimes.append(endtime)
                    totaltime += (endtime - starttime - gaptime)
                    totalt.append(totaltime)
                    totaltime = 0
                    starttime = endtime
                    starttimes.append(starttime)
                else:
                    try:
                        counter = 1
                        totaltime1 = difference
                        gaps = []
                        while totaltime1 < timewidth:
                            g = timetab['START'][i+counter] - timetab['STOP'][i+counter-1]
                            gaps.append(g)
                            m = (timetab['STOP'][i+counter] - timetab['START'][i+counter])
                            totaltime1 += m
                            counter += 1
                        diff = totaltime1 - totaltime
                        endtime = starttime + timewidth
                        for value in gaps:
                            endtime += value
                        endtimes.append(endtime)
                        totaltime += (endtime - starttime)
                        for value in gaps:
                            totaltime = totaltime - value
                        totalt.append(totaltime)
                        totaltime = 0
                        starttime = endtime
                        starttimes.append(starttime)
                    except IndexError:
                        break      
            else: print('difference = 0')
            continue 
        if ((timewidth - totaltime) == interval):
            endtime = starttime + timewidth
            endtimes.append(endtime)
            totaltime += (endtime - starttime)
            totalt.append(totaltime)
            starttime = endtime
            starttimes.append(starttime)
            totaltime = 0
            continue
        if ((timewidth - totaltime) > interval): 
            if (totaltime + interval) < timewidth:
                #need this interval, plus next one
                totaltime += interval
                starttime = timetab['START'][i+1]
            continue
    if starttimes[-1] == endtimes[-1]:
        starttimes = starttimes[:-1]
    return(starttimes, endtimes)
