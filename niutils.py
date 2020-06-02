#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import exp
import matplotlib.patches as patches

#Set colors for all plots
def get_colors():
    d = {'br_colors': ['#ffa600', '#ff6361', '#bc5090',  
                       '#58508d',  '#003f5c'],
         'br_colors_dark': ['#805300', '#803230', '#3D1A2F',
                            '#302C4D', '#000000'],
         'environ_color': '#217525',
         '3c50_color': '#009FE8'}
                       
    return d

#Map colors to ranges
def map_colors(ranges=['--40', '40--60', '60--80', '80--180', '180--'],
               dark=False):
    if dark:
        br_colors = get_colors()['br_colors_dark']
    else:
        br_colors = get_colors()['br_colors']
    assert(len(ranges) == len(br_colors))
    map_colors = {}
    for i in range(len(ranges)):
        map_colors[ranges[i]] = br_colors[i]

    return map_colors

#Check if list tuple, np array, etc
def check_iter(var):
    if type(var) == str:
        return False
    elif hasattr(var, '__iter__'):
        return True
    else:
        return False

#Format bright earth range for legend or txt
def br_label(r):
    split = r.split('--')
    assert(len(split)==2)
    if split[0] == '':
        return r'BR\_EARTH $<{0}^\circ$'.format(split[1])
    elif split[1] == '':
        if split[0] in ['180', 180]:
            return r'BR\_EARTH$=200^\circ$'
        else:
            return r'BR\_EARTH $>{0}^\circ$'.format(split[0])
    else:
        return r'${0}^\circ<$BR\_EARTH$<{1}^\circ$'.format(*split)

#Format sun angle range for legend or txt
def sa_label(r):
    split = r.split('--')
    assert(len(split)==2)

    if split[0] == '':
        return r'SUN\_ANGLE $<{0}^\circ$'.format(split[1])
    else:
        return r'SUN\_ANGLE $>{0}^\circ$'.format(split[0])

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

    if 1 <= phasetup[0] <= phasetup[1]:
        phasetup = (phasetup[0]-1, phasetup[1]-1)
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

def DoubleY(ax, colors=('black', 'black')):
    '''
    Create a double y axis with two seperate colors

    :param ax: axes to modify

    :type ax: matplotlib axes object

    :param colors: 2-tuple of axes colors

    :type colors: tuple length 2

    :returns: two axes, modified original and new y scale
    '''
    if (type(colors) != tuple) or (len(colors) != 2):
        raise TypeError("colors must be 2-tuple")
    ax2 = ax.twinx()
    ax.minorticks_on()
    ax.xaxis.set_ticks_position('both')
    for a in [ax, ax2]:
        a.minorticks_on()
        a.tick_params(direction='in', which='both', labelsize=15)
        a.tick_params('both', length=8, width=1.8, which='major')
        a.tick_params('both', length=4, width=1, which='minor')
    ax.tick_params('y', colors=colors[0], which='both')
    ax2.tick_params('y', colors=colors[1], which='both')
    return ax, ax2

def add_CharErrBar(ax, counts, xfrac, yfrac):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    char_err = np.median(counts**(0.5))

    xpos = [ ((xmax-xmin)*xfrac)+xmin ]
    ypos = [ ((ymax-ymin)*yfrac)+ymin ]

    ax.errorbar(xpos, ypos, yerr=char_err, marker='.', 
                color='xkcd:violet', 
                capsize=4, xerr=0, zorder=8)

    #xpad = (xmax-xmin)*(.05)
    xpad = .1
    ypad = (ymax-ymin)*(.075) + ((char_err) * 1.5)
    """
    patch = patches.FancyBboxPatch(xy=(1, ypos[0]-.5*ypad),
                                   width=.00001, height=ypad,
                                   boxstyle='round',
                                   edgecolor='black', 
                                   facecolor='white', alpha=.6,
                                   zorder=1)
    """
    patch = patches.Rectangle(xy=(xpos[0]-.5*xpad, ypos[0]-.5*ypad),
                              width=xpad, height=ypad,
                              edgecolor='black', facecolor='white',
                              alpha=.6, zorder=6, lw=2)


    ax.add_patch(patch)

    bounds = [patch.get_x(), patch.get_y(),
              patch.get_x()+patch.get_width(),
              patch.get_y()+patch.get_height()]



    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    return ax, bounds


def gaus(x, a, x0, sigma, b):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + b


def two_gaus(x, a_0, x0_0, sigma_0,
                a_1, x0_1, sigma_1, b):

    return a_0*exp(-(x-x0_0)**2/(2*sigma_0**2)) + a_1*exp(-(x-x0_1)**2/(2*sigma_1**2)) + b

def two_lorentzians(x, a_0, x0_0, gamma_0,
                       a_1, x0_1, gamma_1, b):

    return (a_0*(gamma_0)**2 * (1 / (( gamma_0)**2 + (x -x0_0)**2 )) + 
            a_1*(gamma_1)**2 * (1 / (( gamma_1)**2 + (x -x0_1)**2 )) + b)



def component_label(popt_tup, p1_shift, p2_shift):
    if hasattr(popt_tup, 'primary_fwhm'):
        p1_coords = (popt_tup.primary_position+p1_shift[0]*
                     popt_tup.primary_fwhm,
                     popt_tup.primary_amplitude*p1_shift[1]+
                     popt_tup.vertical_shift)
        p2_coords = (popt_tup.secondary_position+p2_shift[0]*
                     popt_tup.secondary_fwhm,
                     popt_tup.secondary_amplitude*p2_shift[1]+
                     popt_tup.vertical_shift)
    else:                                                                       
        p1_coords = (popt_tup.primary_position+p1_shift[0]*
                     popt_tup.primary_sigma,
                     popt_tup.primary_amplitude*p1_shift[1]+
                     popt_tup.vertical_shift)
        p2_coords = (popt_tup.secondary_position+p2_shift[0]*
                     popt_tup.secondary_sigma,
                     popt_tup.secondary_amplitude*p2_shift[1]+
                     popt_tup.vertical_shift)

    return p1_coords, p2_coords                   
