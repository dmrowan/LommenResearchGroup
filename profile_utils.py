#!/usr/bin/env python

import numpy as np
import collections
import pandas as pd
import matplotlib.pyplot as plt
from LCClass import LightCurve
import argparse
from fuzzywuzzy import process
import os
import spectraplots

desc="""
Various profile tools to use for generating pulse profiles
"""

# Produce energy profile over selected energy range in keV
def energy_filtered_profile(lc_input, energy_min, energy_max, ax=None,
                            phase_min=None, phase_max=None, label=True):


    for var_in in [energy_min, energy_max, phase_min, phase_max]:
        if var_in is not None:
            if type(var_in) not in [int, float]:
                raise ValueError(f"{var_in} must be int or float")

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8,4))
        created_fig=True
    else:
        created_fig=False

    pi_min = energy_min * 100
    pi_max = energy_max * 100

    lc = parse_lc_input(lc_input)
    lc.mask(lower_pi=pi_min, upper_pi=pi_max)
    lc.generate()

    ax = lc.plot(ax=ax, label=label)
    plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
    if phase_min is not None and phase_max is not None:
        ax.axvspan(phase_min, phase_max, color='gray', alpha=.2)

    if created_fig:
        plt.show()
    else:
        return ax

#Generate profile and zoom on phase region
def zoom_profile(evt, phase_min, phase_max):

    if type(evt) != str:
        raise ValueError("filename must be string")
    if any( [type(v) not in [float, int] for v in [phase_min, phase_max] ]):
        raise ValueError("phase must be int or float")
    if not os.path.isfile(evt):
        raise FileNotFoundError

    fig, ax = plt.subplots(1, 1, figsize=(8,4))

    lc = LightCurve(evt)
    pi_min = energy_min * 100
    pi_max = energy_max * 100
    lc.mask(lower_pi=pi_min, upper_pi=pi_max)
    lc.generate()

    ax = lc.plot(ax=ax)
    ax.set_xlim(left=phase_min, right=phase_max)
    plt.subplots_adjust(bottom=.08, top=.98, right=.98, left=.15)
    plt.show()

def parse_lc_input(lc_input):
    if type(lc_input) == str:
        filename = lc_input
        if not os.path.isfile(filename):
            raise FileNotFoundError
        lc = LightCurve(filename)
        lc.generate()
    elif type(lc_input) == LightCurve:
        lc = lc_input
    else:
        raise TypeError

    return lc


# Find primary and interpulse phase ranges with user
# defined offpeak phase ranges
def find_phase_ranges(lc_input, off1, off2, nsigma, verbose=True):

    lc = parse_lc_input(lc_input)

    if any( [type(v) not in [float, int] for v in [off1, off2, nsigma] ]):
        raise ValueError("value must be int or float")

    cutofftup = lc.peak_cutoff(off1, off2, nsigma=nsigma)
    if verbose:
        print(f"Min phase primary: {cutofftup.min_phase}")
        print(f"Max phase primary: {cutofftup.max_phase}")
        print(f"Min phase interpulse: {cutofftup.min_phase_ip}")
        print(f"Max phase interpulse: {cutofftup.max_phase_ip}")
    return cutofftup

#Find trailing edge phase ranges
def find_edge(lc_input, off1, off2, nsigma):

    lc = parse_lc_input(lc_input)

    if any( [type(v) not in [float, int] for v in [off1, off2, nsigma] ]):
        raise ValueError("value must be int or float")

    cutofftup = find_phase_ranges(lc, off1, off2, nsigma, verbose=False)
    edge_tuple = collections.namedtuple('edge_tuple', ['peak', 'min', 'max'])
    if cutofftup.min_phase < cutofftup.max_phase:
        tup = edge_tuple(lc.peak_center()[0], cutofftup.min_phase, 
                         cutofftup.max_phase)
    else:
        tup = edge_tuple(lc.peak_center()[0]+1, cutofftup.min_phase, 
                         cutofftup.max_phase+1)
    return tup

#Produce multiple profiles with different energy ranges
def multiple_profiles(evt, energy_ranges, fit_one=False, 
                      fit_two=False, component=None,
                      output=None, show=True):
    if type(evt) != str:
        raise ValueError("filename must be string")
    if type(energy_ranges) not in [list, tuple]:
        raise ValueError("ranges must be entered in lists or tuples")
    for range_pair in energy_ranges:
        if type(range_pair) not in [list, tuple]:
            raise ValueError("ranges must be entered in lists or tuples")
        if len(range_pair) != 2:
            raise ValueError("range must have length 2")
        if any( [type(v) not in [float, int] for v in range_pair] ):
            raise ValueError("value must be int or float")

    sourcename = process.extract(evt, ['PSR B1821-24', 'PSR B1937+21'],
                                 limit=1)[0][0]
    if fit_one and component is None:
        raise ValueError("component must be defined for fit")

    fig, ax = plt.subplots(len(energy_ranges), 1, figsize=(8, len(energy_ranges)*3.5))
    ratios = []
    popts  = []
    for i in range(len(energy_ranges)):
        a = ax.reshape(-1)[i]
        if fit_one or fit_two:
            lc = LightCurve(evt)
            lc.mask(lower_pi=energy_ranges[i][0]*100, upper_pi=energy_ranges[i][1]*100)
            lc.generate(bs=.005)
            name = lc.name
            if fit_one:
                a, popt, ratio = lc.fit_gauss(component, ax=a)
                ratios.append(ratio)
            else:
                a, popt = lc.fit_two_gauss(ax=a, annotate=False)
            popts.append(popt)
        else:
            lc = LightCurve(evt)
            name = lc.name
            a = energy_filtered_profile(evt, energy_ranges[i][0], 
                                        energy_ranges[i][1], ax=a,
                                        label=False)
        a.text(.95, .95, 
               f"{energy_ranges[i][0]}"+r'$-$'+ f"{energy_ranges[i][1]} keV",
               ha='right', va='top', fontsize=20, transform=a.transAxes)
        a.set_xlabel("")
        if i != len(energy_ranges)-1:
            a.tick_params(labelbottom=False)

    ax.reshape(-1)[0].text(.05, .95, f"{name}", ha='left',
                           va='top', transform=ax.reshape(-1)[0].transAxes,
                           fontsize=20)
    ax.reshape(-1)[len(energy_ranges)-1].set_xlabel("Phase", fontsize=20)
    fig.text(.02, .5, "Photon Counts", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    plt.subplots_adjust(hspace=0, bottom=.08, top=.94, right=.98)
    if output is None:
        if show:
            plt.show()
            plt.close()
    else:
        fig.savefig(output, dpi=2000)

    if fit_one:
        return popts, ratios
    elif fit_two:
        return popts

#Test different energy range fits to find optimal spectra range
def find_min_energy(evt, component):
    if type(evt) != str:
        raise ValueError("filename must be string")
    if type(component) != str:
        raise ValuerError ("component must be string")

    sourcename = process.extract(evt, ['PSR B1821-24', 'PSR B1937+21'], 
                                 limit=1)[0][0]

    attempted_energies = [ round(a, 2) for a in np.arange(0, 2, 0.1) ]
    successful_energies = []
    fig, ax = plt.subplots(len(attempted_energies), 1, figsize=(8, len(attempted_energies)*4))
    for i in range(len(attempted_energies)):
        eng = attempted_energies[i]
        lc = LightCurve(evt)
        a = ax.reshape(-1)[i]
        print(eng)
        lc.mask(lower_pi=0, upper_pi=eng*100)
        lc.generate()
        a, chisq, p, ratio = lc.fit_gauss(component, ax=a)
        if ratio > 1:
            successful_energies.append(eng)
    plt.show()
    return successful_energies

def shifting_profile(evt, output=None):
    if type(evt) != str:
        raise ValueError("filename must be string")

    sourcename = process.extract(evt, ['PSR B1821-24', 'PSR B1937+21'],
                                 limit=1)[0][0]

    en = 0.5
    energy_ranges = []
    while en < 9.5:
        energy_ranges.append( (en, en+0.5) )
        en += 0.5

    popts = multiple_profiles(evt, energy_ranges, fit_two=True, show=False)
    print(popts)
    primary_loc = [popt[1]-1 for popt in popts]
    primary_fwhm = [ popt[2]*2.355 for popt in popts ]
    interpulse_loc = [ popt[4]-1 for popt in popts ]
    interpulse_fwhm = [ popt[5]*2.355 for popt in popts ]

    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    plt.subplots_adjust(top=.98, right=.98, hspace=0)
    ax[0].scatter([r[0] for r in energy_ranges], primary_loc, color='xkcd:violet')
    ax[1].scatter([r[0] for r in energy_ranges], interpulse_loc, color='xkcd:violet')
    ax[0] = spectraplots.plotparams(ax[0])
    ax[1] = spectraplots.plotparams(ax[1])

    #ax[0].set_xlabel("Energy (keV)", fontsize=20)
    ax[1].set_xlabel("Energy (keV)", fontsize=20)

    ax[0].set_ylabel("Center Phase", fontsize=20)
    ax[1].set_ylabel("Center Phase", fontsize=20)
    
    ax[0].text(.95, .95, "Primary Pulse", fontsize=20, 
               ha='right', va='top', transform=ax[0].transAxes)
    ax[0].text(.95, .82, sourcename, fontsize=20, 
               ha='right', va='top', transform=ax[0].transAxes)
    ax[1].text(.95, .95, "Interpulse", fontsize=20,
               ha='right', va='top', transform=ax[1].transAxes)
    ax[1].text(.95, .82, sourcename, fontsize=20,
               ha='right', va='top', transform=ax[1].transAxes)

    for i in [0,1]:
        text = ax[i].text(.2, .1, "Example Plot", fontsize=50, 
                   ha='left', va='bottom', rotation=25, color='red',
                   transform=ax[i].transAxes)
        text.set_alpha(0.6)

    if output is None:
        plt.show()
    else:
        fig.savefig(output)



