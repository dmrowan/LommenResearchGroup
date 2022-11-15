#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
import csv
import random

desc="""
Plots outliers for the crab as subplots
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def plot(timewidth):
    intints = pd.read_csv('2pulseintint_%s.csv'%(timewidth), delimiter = '\t', header=None)
    curves = pd.read_csv('2pulsecurve_%s.csv'%(timewidth), delimiter = '\t', header=None)
    largest_column_count = 0
    data_file_delimiter = '\t'
    with open('2pulsehist_%s.csv'%(timewidth), 'r') as temp_f:
    # Read the lines
        lines = temp_f.readlines()
        for l in lines:
            # Count the column count for the current line
            column_count = len(l.split(data_file_delimiter)) + 1
            # Set the new most column count
            largest_column_count = column_count if largest_column_count < column_count else largest_column_count
   # Close file
    temp_f.close()
   # Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
    column_names = [i for i in range(0, largest_column_count)]
    hists = pd.read_csv('2pulsehist_%s.csv'%(timewidth), delimiter = '\t', header=None, names = column_names)
    for n in range(len(intints)):
        if len(np.array(hists.loc[n].dropna())) < 100:
            hists = hists.drop(n)      
    #plot subplots of all outliers
    row = 4
    col = 4

    log.info("Plot all curve fits")
    fig, ax = plt.subplots(row, col, sharex = 'col', figsize = (13, 7))
    i = 0
    j = 0
    index = hists.index.tolist()
   # for n in random.sample(index, 9):
    for n in index[:16]:
        if (j > (col-1)):
            j = 0
            i += 1
        yvals, xlims = np.histogram(np.array(hists.loc[n].dropna()),bins=200) # finds heights and sides of each bin, no plot
        xvals = xlims[:-1] + np.diff(xlims)/2
        ax[i, j].plot(xvals, yvals)
        ax[i, j].plot(xvals, gauss2(xvals, *curves.loc[n]))
        ax[i, j].set_title('%f, %f'%(intints[0].loc[n], intints[1].loc[n]), fontsize = 8)
        j += 1

    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')
    fig.suptitle('Fitted Profiles for Integration Time of %s s'%(timewidth))


    plt.savefig('fittedprofiles_%s.png' %(timewidth))

    plt.clf()
    

time =  [90]
for t in time:
    plot(t)
