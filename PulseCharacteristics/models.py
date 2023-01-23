#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from functions import *
from scipy import stats, interpolate
from copy import deepcopy

desc="""

Runs tests on model for the distribution of integrated intensities. 

"""

def model(n_pulses, model, p, plot=False, modelplot=False, split=None):
 
    if split is None:
        intint = pd.read_csv('intdata/crabintdata_%s.txt' %n_pulses, header = None)
        intint = list(intint[0])
        intint = [x for x in intint if x > 0]
        intint = [x for x in intint if x < 0.2]
        intint = np.array(intint)
        binwidths = list(np.linspace(0, 0.12, 100))

    else:
        intint = split[0]
        if split[1] < 5:
            s = 100
        elif split[1] < 16:
            s = 50
        else:
            s = 25
        binwidths = list(np.linspace(0, 0.12, s))
    #plt.hist(intint, binwidths)

    width = 0.05
    intint = np.array(intint) 
    xvals, yvals = hist_to_curve(intint, binwidths)
    yvals_full = yvals

    y_error = np.sqrt(yvals)

    #Models
    if model == 'Gaussian':
        #create model, convolve it with itself n_pulses times
        a = p[0]
        b = p[1]
        c = p[2]
        xvalues = np.arange(0, 0.2, 0.001) #arbitrary range of x values
        yvalues = gauss(xvalues, a, b, c)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c))
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()

        #scale so that peak has a value of 1
        #yvals = yvals/max(yvals)
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'Log normal':
        #create model, convolve it with itself n_pulses times
        a = p[0]
        b = p[1]
        c = p[2]
        xvalues = np.arange(0.01, 3, 0.01) #arbitrary range of x values
        yvalues = lognormal(xvalues, a, b, c)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c)) 
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()
        
        #scale so that peak has a value of 1
        #yvals = yvals/max(yvals)
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'Power law':
        #create model, convolve it with itself n_pulses times
        a = p[0]
        b = p[1]
        xvalues = np.arange(0.01, 0.04, 0.0001) #arbitrary range of x values
        yvalues = power(xvalues, a, b)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c))
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()
    
        #scale so that peak has a value of 1
        #yvals = yvals/max(yvals)
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'double power law':
        #create model, convolve it with itself n_pulses times
        a = p[0] #use best from single power law
        b = p[1] #use best from single power law
        c = p[2]
        d = p[3]
        xvalues = np.arange(0.01, 0.04, 0.0001) #arbitrary range of x values
        yvalues = power2(xvalues, a, b, c, d)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c))
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()

    #scale and shift x axis
    xvalues_c = matchtodata(xvals, yvals, xvalues_c, yvalues_c)
    print(max(yvals), max(yvalues_c))

    #resample model curve at x values of data
    try:
        f = interpolate.interp1d(xvalues_c, yvalues_c)
        yvalues_r = f(xvals)
    except ValueError:
        if max(xvals) > max(xvalues_c):
            xvals_temp = []
            for i in xvals:
                if i < max(xvalues_c):
                    xvals_temp.append(i)
        f = interpolate.interp1d(xvalues_c, yvalues_c)
        yvalues_r = f(xvals_temp)
        xvals = xvals_temp
        yvals = yvals[:len(xvals)]

    if plot == True:
        plt.plot(xvals, yvals, label='Data')
        plt.plot(xvals, yvalues_r, label='Model')
        plt.xlabel('Integrated intensity (counts/pulse)')
        plt.ylabel('# of Profiles')
        plt.title('Integrated intensity distribution, %s'%model)
        plt.legend()
        plt.show()

    return(xvals, yvals, yvalues_r, y_error, yvals_full)

def matchtodata(xvals, yvals, xvalues_c, yvalues_c):
    #scale the x axis using the FWHM of the data and the model
    fwhm_d = fwhm(xvals, yvals)
    fwhm_m = fwhm(xvalues_c, yvalues_c)
    scale = fwhm_d/fwhm_m
    xvalues_c = xvalues_c*scale

    try:
        f = interpolate.interp1d(xvals, yvals)
        yvals_temp = f(xvalues_c)
    except ValueError:
        xvals_temp = deepcopy(xvals)
        yvals_temp = deepcopy(yvals)
        tempx = []
        tempy = []
        if max(xvalues_c) > max(xvals):
            for i in xvalues_c:
                if i > max(xvals):
                    tempx.append(i)
                    tempy.append(0)
            xvals_temp = np.concatenate((xvals_temp, np.array(tempx)), axis=None)
            yvals_temp = np.concatenate((yvals_temp, np.array(tempy)), axis=None)
        tempx = []
        tempy = []
        if min(xvalues_c) < min(xvals):
            for i in xvalues_c:
                if i < min(xvals):
                    tempx.append(i)
                    tempy.append(0)
            xvals_temp = np.concatenate((np.array(tempx), xvals_temp), axis=None)
            yvals_temp = np.concatenate((np.array(tempy), yvals_temp), axis=None)
        
        f = interpolate.interp1d(xvals_temp, yvals_temp)
        yvals_new = f(xvalues_c)

    #find the x shift by convolving data and model plots
    m = convolve1(yvalues_c, yvals_new)
    shift = xvalues_c[m] - xvalues_c[0]
    xvalues_c = xvalues_c-shift
    
    #return shifted and scaled x values
    return(xvalues_c)

def chisq(yvals, yvalues_r, y_error):
   
    #normalize
    #total =  np.trapz(yvals)
    y_error = np.sqrt(yvals)
    #yvals = yvals/total
    #yvalues_r = yvalues_r/total
    #y_error = y_error/total

    chisqr = 0
    for i in range(len(yvals)):
        chisqr = chisqr + ((yvals[i]-yvalues_r[i])**2)/(y_error[i]**2) 
  
    csq = stats.chisquare(yvals, yvalues_r)
    return(csq[0], chisqr)

def likelihood(yvals, yvalues_c, y_error, yvals_full):

    #normalize
    total =  np.trapz(yvals)
    yvals_p = yvals/total
    model_p = yvalues_c/total
    
    #y_error = y_error/total
    #coef = 1/((2*np.pi*np.prod(y_error**2))**(len(yvals)/2))
    #lh = np.exp(-sum(((yvals - yvalues_c)/y_error)**2)/2)

    #print(yvals_full)
    #print(model_p)

    lh = 0
    for i in range(len(yvals_full)):
        lh += yvals_full[i]*np.log10(model_p[i])
    
    return(lh, 10**(lh))

def residuals(yvals, yvalues_c):
    
    res = yvals - yvalues_c
    plt.plot(res)
    plt.show()

def testmodel(n_pulses, model_type, parameters):

    print("Testing %s for %s pulses:"%(model_type, n_pulses))
    x, y, y_r, y_err, y_full = model(n_pulses, model_type, parameters, plot=True, modelplot=False)
    
    print("The chi squared value is", chisq(y, y_r, y_err))
    print("The likelihood is", likelihood(y, y_r, y_err, y_full))



N = 15
#models = ["Log normal", "Gaussian", "Power law"]
#models = [["Log normal", [1, -1.3, 0.8]]]
#models = [['Gaussian', [1, 0.1, 0.01]]]
#models = [['Power law', [1, 5.25]]]
models = [["Log normal", [1, -1.3, 0.8]], ['Gaussian', [1, 0.1, 0.01]], ['Power law', [1, 5.25]]]
#models = [['power law', [1, 5.25]], ['double power law', [1, 5.25, 1, 6]]]
for m in models:
    testmodel(N, m[0], m[1])
