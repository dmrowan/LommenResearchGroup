#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

desc="""
Plots the convolution of the peak, then does a gaussian curve fit, plots it and returns peak amplitude
"""
def gauss(x, a, b, c, d):
    return(a*np.exp(-((x-b)/c)**2)+d)

def main():
    fname = '1937_events.evt'

    log.info('Read in table')
    tab = Table.read('1937_events.evt', hdu=1)

    log.info("Limit data")
    n = int(input("What fraction of data do you want? (nth of data)"))
    limitedtab = np.array_split(tab['PULSE_PHASE'], n)
    phase = np.array(limitedtab[0])
    yvals, xlims, _ = plt.hist(phase,bins=255)
    xvals = xlims[:-1] + np.diff(xlims)/2
    
    width=0.05
    x = xvals
    template = np.exp(-((x)/width)**2)

    convo = []
    for i in range(len(yvals)):
        convo.append(np.sum(yvals*np.roll(template,i)))

  #  plt.plot(xvals, convo)
    m = np.max(convo)
    maxloc = xvals[convo.index(m)]
    print(m, maxloc)

    popt, pcov = curve_fit(gauss, xvals, yvals, p0= [500,0.037,0.05, 200],bounds=([0,0,0,100],[600, 0.2, 0.2, 400]))
    print(popt)
    plt.plot(xvals, gauss(xvals, *popt))

 #   p =	int(input("Which profile do you want? (valid input: 0-(n-1))"))
 #   b = int(input("How many bins?"))


    plt.show()

if __name__ == '__main__':
    main()

