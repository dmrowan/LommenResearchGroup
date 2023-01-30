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

"""
A set of functions for fitting models to an intensity histogram and doing statistics on it
Note that here I will refer to the calculated Crab integrated intensity values as the "data" and the model being fitted to it as the "model"

Also note that because the pulse profiles are the sum of N pulses, cannot directly fit a model to it, as we want the underlying distribution
This means we cannot use a basic curve fit as we need to convolve the with itself N times first, and then fit that convolved model to the data
There is no clean way to do this (that we know of), so we make a template model, convolve it N times, and scale its height and scale and shift it along the horizontal axis until it fits the data. The input model parameters are for the underlying, not convolved model


"""

def model(n_pulses, model, p, plot=False, modelplot=False, split=None): 
    #fits the given model to the data; p=model parameters; can plot; 
    #split is the number of sections you are splitting the data into (e.g. to see any variation) 
    #note that if using split, read in data and split data into sections: split = [[data], N sections]
 
    # Read in data and set up bins
    if split is None: #if using full data set of intensities, read it in, remove outliers, set bins for histogram (check it looks right beforehand)
        intint = pd.read_csv('intdata/crabintdata_%s.txt' %n_pulses, header = None)
        intint = list(intint[0])
        intint = [x for x in intint if x > 0]
        intint = [x for x in intint if x < 0.2]
        intint = np.array(intint)
        binwidths = list(np.linspace(0, 0.12, 100))

    else: #if splitting data, split = [[intensity values for section of data set as an array], N sections you split data set into]
        intint = split[0] #get data values 
        if split[1] < 5: #make sure bins are not too wide/narrow; these are approximate ranges of number of bins that worked for me
            s = 100
        elif split[1] < 16:
            s = 50
        else:
            s = 25
        binwidths = list(np.linspace(0, 0.12, s)) #set up bins for histogram (check that these look right)

    # Convert the data histogram to a curve so that a curve can be fit to it
    width = 0.05
    intint = np.array(intint) 
    xvals, yvals = hist_to_curve(intint, binwidths)
    yvals_full = yvals

    y_error = np.sqrt(yvals) #calculate error in y values; Poisson error 

    # Set up the template models
    # If adding any new, copy any as a template and modify for number of parameters 
    if model == 'Gaussian':
        #create model, convolve it with itself n_pulses times
        a = p[0]
        b = p[1]
        c = p[2]
        xvalues = np.arange(0, 0.2, 0.001) #arbitrary range of x values, just make sure no major features of model are being cut off
        yvalues = gauss(xvalues, a, b, c) #get y values using function gauss (in functions.py)
        yvalues_c = convolve(yvalues, n_pulses) #convolved y values; convolve model with itself n_pulses times
        xvalues_c = np.arange(len(yvalues_c)) #get a new set of x values to go with the new y values (len(yvalues) != len(yvalues_c))
        if modelplot == True: #can show a plot of the underlying model
            plt.plot(xvalues, yvalues)
            plt.show()

        #Scale so that peak has the same amplitude as the data (scale model along y axis)
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'Log normal':
        a = p[0]
        b = p[1]
        c = p[2]
        xvalues = np.arange(0.01, 3, 0.01) 
        yvalues = lognormal(xvalues, a, b, c)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c)) 
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()
        
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'Power law':
        a = p[0]
        b = p[1]
        xvalues = np.arange(0.01, 0.04, 0.0001) 
        yvalues = power(xvalues, a, b)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c))
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()
    
        yvalues_c = yvalues_c/max(yvalues_c)*max(yvals)

    if model == 'double power law':
        a = p[0] #use best from single power law
        b = p[1] #use best from single power law
        c = p[2]
        d = p[3]
        xvalues = np.arange(0.01, 0.04, 0.0001) 
        yvalues = power2(xvalues, a, b, c, d)
        yvalues_c = convolve(yvalues, n_pulses)
        xvalues_c = np.arange(len(yvalues_c))
        if modelplot == True:
            plt.plot(xvalues, yvalues)
            plt.show()

    # Scale and shift the model along the x axis
    xvalues_c = matchtodata(xvals, yvals, xvalues_c, yvalues_c)

    # Resample model at x values of data (because currently len(yvalues_c) > len(yvals))
    try:
        f = interpolate.interp1d(xvalues_c, yvalues_c)
        yvalues_r = f(xvals) #resampled model y values; resample at xvals so that len(yvalues_r) = len(yvals)
    except ValueError: #if interpolation fails, it means some values do not overlap, get rid of those
        if max(xvals) > max(xvalues_c):
            xvals_temp = [] 
            for i in xvals:
                if i < max(xvalues_c):
                    xvals_temp.append(i)
        f = interpolate.interp1d(xvalues_c, yvalues_c) 
        yvalues_r = f(xvals_temp) #resample again at the x values that DO overlap fully with model
        xvals = xvals_temp
        yvals = yvals[:len(xvals)] #make sure number of yvals and xvals match up

    if plot == True: #can plot the model overlaid over the data; very useful to making sure everything worked
        plt.plot(xvals, yvals, label='Data')
        plt.plot(xvals, yvalues_r, label='Model')
        plt.xlabel('Integrated intensity (counts/pulse)')
        plt.ylabel('# of Profiles')
        plt.title('Integrated intensity distribution, %s'%model)
        plt.legend()
        plt.show()

    # Returns the x values for both, the data y values, the fitted model y values, the error in yvals, and yvals without any scaling(?)
    return(xvals, yvals, yvalues_r, y_error, yvals_full)

def matchtodata(xvals, yvals, xvalues_c, yvalues_c):
    #scale the x axis using the full width half max of the data and the model; shift model along x axis
    
    fwhm_d = fwhm(xvals, yvals) #fwhm of data
    fwhm_m = fwhm(xvalues_c, yvalues_c) #fwhm of model
    scale = fwhm_d/fwhm_m 
    xvalues_c = xvalues_c*scale #scale x values by ratio of fwhm of data and model

    #resample data at x values of the convolved model (note len(yvalues_c) > len(yvals))
    try:
        f = interpolate.interp1d(xvals, yvals)
        yvals_temp = f(xvalues_c)
    except ValueError: #error if model xvalues_c wider than data xvals (there is a region of no overlap)
        
        #make temporary deep copies of the data
        xvals_temp = deepcopy(xvals)
        yvals_temp = deepcopy(yvals)
        tempx = []
        tempy = []
        if max(xvalues_c) > max(xvals): #if xvalues_c extends further right than xvals, append values into xvals and zeros into yvals
            for i in xvalues_c: #for each value in xvalues_c 
                if i > max(xvals): #if that value is greater than the highest x values the data goes to
                    tempx.append(i) #append i into a temporary list of x values
                    tempy.append(0) #append 0 into a temporary list of y values 
            xvals_temp = np.concatenate((xvals_temp, np.array(tempx)), axis=None) #add tempx onto end of xvals_temp
            yvals_temp = np.concatenate((yvals_temp, np.array(tempy)), axis=None) #add tempy onto end of yvals_temp
        
        #repeat, but for any values of xvalues_c that extend further left than xvals 
        tempx = [] 
        tempy = []
        if min(xvalues_c) < min(xvals):
            for i in xvalues_c:
                if i < min(xvals):
                    tempx.append(i)
                    tempy.append(0)
            xvals_temp = np.concatenate((np.array(tempx), xvals_temp), axis=None)
            yvals_temp = np.concatenate((np.array(tempy), yvals_temp), axis=None)
        
        #resample data at x values of the convolved model 
        f = interpolate.interp1d(xvals_temp, yvals_temp)
        yvals_new = f(xvalues_c) 

    #Shift the model along the x axis so that the peaks match up
    #find the x shift by convolving data and model plots
    m = convolve1(yvalues_c, yvals_new) 
    shift = xvalues_c[m] - xvalues_c[0] #the value by which the model needs to be shifted
    xvalues_c = xvalues_c-shift
    
    #return shifted and scaled x values
    return(xvalues_c)



def chisq(yvals, yvalues_r, y_error): #calculate the chi squared given the data, model, and error in data (but that might be recalculated here)
   
    #Normalize (is this necessary?)
    #total =  np.trapz(yvals)
    y_error = np.sqrt(yvals)
    #yvals = yvals/total
    #yvalues_r = yvalues_r/total
    #y_error = y_error/total

    #Calculate the chi squared directly
    chisqr = 0
    for i in range(len(yvals)):
        chisqr = chisqr + ((yvals[i]-yvalues_r[i])**2)/(y_error[i]**2) 
  
    #Calculate the chi squared using the stats function
    csq = stats.chisquare(yvals, yvalues_r)
    return(csq[0], chisqr) #returns both chi squares calculated both ways; should match up

def likelihood(yvals, yvalues_c, y_error, yvals_full): #calculate the log likelihood

    #Normalize
    total =  np.trapz(yvals)
    yvals_p = yvals/total
    model_p = yvalues_c/total
    
    #Using the equation (the Andrea way to calculate it)
    #y_error = y_error/total
    #coef = 1/((2*np.pi*np.prod(y_error**2))**(len(yvals)/2))
    #lh = np.exp(-sum(((yvals - yvalues_c)/y_error)**2)/2)

    #print(yvals_full)
    #print(model_p)

    #Using Kent's way of calculating it (in the emails)
    lh = 0
    for i in range(len(yvals_full)):
        lh += yvals_full[i]*np.log10(model_p[i])
    
    return(lh, 10**(lh)) #returns log likelihood and likelihood

def residuals(yvals, yvalues_c): #plot the residuals of data - model
    
    res = yvals - yvalues_c
    plt.plot(res)
    plt.show()

def testmodel(n_pulses, model_type, parameters): #fit model to the data, given N pulses per profile, name of model, and model parameters

    print("Testing %s for %s pulses:"%(model_type, n_pulses))
    #runs function model using N pulses, name of model, parameters, and can also show plots
    #outputs are x values, data y values, fitted model y values, error on model y values, and data y values without any scaling(?)
    x, y, y_r, y_err, y_full = model(n_pulses, model_type, parameters, plot=True, modelplot=False) 
    
    print("The chi squared value is", chisq(y, y_r, y_err)) #chi squared for the model
    print("The likelihood is", likelihood(y, y_r, y_err, y_full)) #likelihood of the model


# For N pulses per profile, test models with the specified parameters
# The current parameters are the best I've found using file modelstest.py 
N = 15
#models = [["Log normal", [1, -1.3, 0.8]]]
#models = [['Gaussian', [1, 0.1, 0.01]]]
#models = [['Power law', [1, 5.25]]]
models = [["Log normal", [1, -1.3, 0.8]], ['Gaussian', [1, 0.1, 0.01]], ['Power law', [1, 5.25]]]
#models = [['power law', [1, 5.25]], ['double power law', [1, 5.25, 1, 6]]]
for m in models:
    testmodel(N, m[0], m[1])
