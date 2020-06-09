#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt

def convert_time(time):

    timezero = datetime.datetime(year=2014, month=1,
                                 day=1, hour=0, minute=0, second=0)
    new_time = timezero+datetime.timedelta(seconds=time)

    return new_time

#info from feb3_energies_fitting.py
startTime = 390+1.92224*10**8
##uncertainties from Noah's data
lowE_frac_unc = 0.8580811012375814
midE_frac_unc = 0.16515273968538902
highE_frac_unc = 0.12636021430062486

#x-intercepts found on mathematica
lowE_xint = 9.57002
lowE_time = startTime + lowE_xint
midE_xint = 4.54309
midE_time = startTime + midE_xint
highE_xint = 2.13372
highE_time = startTime + highE_xint

dateTime_lowE = convert_time(lowE_time)
dateTime_midE = convert_time(midE_time)
dateTime_highE = convert_time(highE_time)

print(f'The dateTime for the high energy crossing is: {dateTime_highE}, with an uncertainty of {highE_frac_unc*highE_xint} seconds')
print(f'The dateTime for the mid energy crossing is: {dateTime_midE}, with an uncertainty of {midE_frac_unc*midE_xint} seconds')
print(f'The dateTime for the low energy crossing is: {dateTime_lowE}, with an uncertainty of {lowE_frac_unc*lowE_xint} seconds')
