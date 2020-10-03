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

desc="""
Plots outliers for the crab
"""

def gauss(x, a, m, s, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+d)

def gauss2(x, a, m, s, b, c, e, d):
    return((a*np.exp(-(((x - m)/s)**2)/2))+(b*np.exp(-(((x - c)/e)**2)/2))+d)

def outliersplot(scale, timewidth):
    data = pd.read_csv('upperoutliers_10.csv', delimiter = '\t', header=None)
    data.columns = ['intint', 'curve', 'hist']
    print(data.head(10))
    intints = list(data['intint'])
    print(intints)
    curves = list(data['curve'])
    print(curves)
    hists = list(data['hist'])

    #plot subplots of all outliers
    row = 5
    col = 5

    log.info("Plot all curve fits")
    fig, ax = plt.subplots(row, col, sharex = 'col')
    i = 0
    j = 0
    for n in range(len(intints)):
        if (j > (col-1)):
            j = 0
            i += 1
        ax[i, j].hist(np.array(hists[n]), bins = 200)
        #ax[i, j].plot(upperplot[n], gauss(upperplot[n], *curves[n]))
        #ax[i, j].legend(intints[n], loc='upper right')
        j += 1

    fig.text(0.5, 0.04, 'Pulse Phase', ha='center')
    fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')

    plt.savefig('%s_outliers_%s.png' % (scale, timewidth))
    plt.clf()


outliersplot('upperoutliers', '10')

