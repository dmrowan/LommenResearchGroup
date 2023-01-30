#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate
from functions import *
from models import *
import seaborn as sns

"""
Given range of parameters, uses a heatmap or curve to find the best fit parameters for a convolved model fit to an intensity histogram using chi squared

As we use convolved models we cannot directly use curve fit to find best fit values of the model to the data, so with the curve or heatmap, I selected the best models by eye. Make sure to always check the plot to make sure significant features of the model are not being cut off, as you may get false results for a set of parameters

"""

#Select model you are testing
modelname = 'log normal' 


#Model parameters - try to minimize number of free parameters as much as possible
#Remember, model will be arbitrarily scaled along both axes, and shifted along x axis 
#Refer to models.py and functions.py for functions
if modelname == 'log normal': 
    a = 1
    #m = np.linspace(-3, -1.5, 10) 
    m = -1.3
    s = np.linspace(0.1, 1.5, 10)
if modelname == 'gaussian':
    a = 1
    m = 0.1
    s = np.linspace(0.01, 0.1, 20)
if modelname == 'power law':
    a = 1
    s = np.linspace(2, 7, 20) 
if modelname == 'double power law':
    a = 1
    b = 5.25
    m = np.linspace(0, 3, 10)
    s = np.linspace(5, 7, 10)

#If testing log normal with two free parameters, uncomment section and heatmap to select best parameters
#And comment out the "if log normal" if statement below 
"""
if modelname != 'power law' and modelname != 'gaussian': #power law and gaussian only have one free parameter each
    c = []
    for i in s: #for one of the parameters
        row = [] #create rows of the heatmap by varying parameters
        for j in m: #the other parameter
            if modelname != 'double power law': #if double power law, only need the two free parameters
                x, y, y_r, y_err = model(15, modelname, [a, j, i], plot=False)
            else: #if log normal, need two set and two free
                x, y, y_r, y_err = model(15, modelname, [a, b, j, i], plot=False) #run model from models.py on it
            chi2 = chisq(y, y_r) #calculate chi squared using chisq in models.py
            row.append(chi2)
            #print(i, j, chi2)
            #if chi2<2:
            #    print(i, j, chi2)
        c.append(row)

    #Plot as heatmap using seaborn 
    m = [round(i, 2) for i in m]
    s = [round(i, 2) for i in s]
    ax = sns.heatmap(c, cmap=sns.cm.rocket_r, xticklabels = m, yticklabels = s)
    ax.set(xlabel='mu', ylabel='sigma')
    plt.show()
"""

#If using one free parameter plus two set for log normal, do the same but plot as a curve using matplotlib
if modelname == 'log normal':
    c = []
    for i in s:
        x,y,y_r,y_err = model(15, modelname, [a, m, i], plot=False)
        chi2 = chisq(y,y_r)
        c.append(chi2)
    plt.plot(s, c, '.')
    plt.xlabel('s')
    plt.ylabel('chi squared')
    plt.show()
else: #if one free parameter, one set, use this
    c = []
    for i in s:
        x,y,y_r,y_err = model(15, modelname, [a, i], plot=False)
        chi2 = chisq(y,y_r)
        c.append(chi2)
    plt.plot(s, c, '.')
    plt.xlabel('b')
    plt.ylabel('chi squared')
    plt.show()
