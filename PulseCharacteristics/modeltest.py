#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate
from functions import *
from models import *
import seaborn as sns

modelname = 'log normal' #model being tested

if modelname == 'log normal': 
    a = 1
    #m = np.linspace(-3, -1.5, 10)
    m = -1.3
    s = np.linspace(0.1, 1.5, 20)
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
"""
if modelname != 'power law' and modelname != 'gaussian':
    c = []
    for i in s:
        row = []
        for j in m:
            if modelname != 'double power law':
                x, y, y_r, y_err = model(15, modelname, [a, j, i], plot=False)
            else:
                x, y, y_r, y_err = model(15, modelname, [a, b, j, i], plot=False)
            chi2 = chisq(y, y_r)
            row.append(chi2)
            #print(i, j, chi2)
            #if chi2<2:
            #    print(i, j, chi2)
        c.append(row)

    m = [round(i, 2) for i in m]
    s = [round(i, 2) for i in s]
    ax = sns.heatmap(c, cmap=sns.cm.rocket_r, xticklabels = m, yticklabels = s)
    ax.set(xlabel='mu', ylabel='sigma')
    plt.show()
"""
if modelname == 'log normal':
    c = []
    for i in s:
        x,y,y_r,y_err = model(15, modelname, [a, m, i], plot=False)
        chi2 = chisq(y,y_r)
        c.append(chi2)
    plt.plot(s, c, '.')
    plt.xlabel('s')
    plt.ylabel('chi square')
    plt.show()
else:
    c = []
    for i in s:
        x,y,y_r,y_err = model(15, modelname, [a, i], plot=False)
        chi2 = chisq(y,y_r)
        c.append(chi2)
    plt.plot(s, c, '.')
    plt.xlabel('b')
    plt.ylabel('chi square')
    plt.show()
