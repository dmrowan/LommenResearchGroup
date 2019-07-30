#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import exp

def round_sigfigs(val, err):
    idx = len(str(err))
    skip1=True
    for i in range(len(str(err))):
        if str(err)[i] in ['-', '0', '.'] and skip1:
            continue
        elif str(err)[i] == '1' and skip1:
            skip1=False
            continue
        else:
            idx=i
            break
    err_rounded = round(err, idx-1)
    val_rounded = round(val, idx-1)
    return val_rounded, err_rounded


#Fix phasetups that overlap 0
def phase_correction(phasetup):
    if type(phasetup) != tuple:
        raise TypeError
    if phasetup[1] <= phasetup[0]:
        phasetup = (phasetup[0], phasetup[1]+1)
    return phasetup

#Basic ploting parameters
def plotparams(ax):
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=15)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
    return ax


def gaus(x, a, x0, sigma, b):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + b


def two_gaus(x, a_0, x0_0, sigma_0,
                a_1, x0_1, sigma_1, b):

    return a_0*exp(-(x-x0_0)**2/(2*sigma_0**2)) + a_1*exp(-(x-x0_1)**2/(2*sigma_1**2)) + b

def two_lorentzians(x, a_0, x0_0, gamma_0,
                       a_1, x0_1, gamma_1, b):

    return (a_0*(gamma_0)**2 * (1 / (( gamma_0)**2 + (x -x0_0)**2 )) + 
            a_1*(gamma_1)**2 * (1 / (( gamma_1)**2 + (x -x0_1)**2 )) + b)

