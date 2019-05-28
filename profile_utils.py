#!/usr/bin/env python

import numpy as np
import collections
import pandas as pd
import matplotlib.pyplot as plt
from LCClass import LightCurve
import argparse
import os

desc="""
Various profile tools to use for generating pulse profiles
"""

# Produce energy profile over selected energy range in keV
def energy_filtered_profile(evt, energy_min, energy_max, ax=None,
                            phase_min=None, phase_max=None):
    if type(evt) != str:
        raise ValueError("filename must be string")
    for var_in in [energy_min, energy_max, phase_min, phase_max]:
        if var_in is not None:
            if type(var_in) not in [int, float]:
                raise ValueError(f"{var_in} must be int or float")
    if not os.path.isfile(evt):
        raise FileNotFoundError

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8,4))
        created_fig=True
    else:
        created_fig=False

    lc = LightCurve(evt)
    pi_min = energy_min * 100
    pi_max = energy_max * 100
    lc.mask(lower_pi=pi_min, upper_pi=pi_max)
    lc.generate()

    ax = lc.plot(ax=ax)
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
    tup = edge_tuple(lc.peak_center()[0], cutofftup.min_phase, 
                     cutofftup.max_phase)
    return tup

def multiple_profiles(evt, energy_ranges):
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

    fig, ax = plt.subplots(len(energy_ranges), 1, figsize=(8, len(energy_ranges)*4))
    for i in range(len(energy_ranges)):
        a = ax.reshape(-1)[i]
        a = energy_filtered_profile(evt, energy_ranges[i][0], 
                                    energy_ranges[i][1], ax=a)
        a.text(.05, .95, f"{energy_ranges[i][0]}"+r'$-$'+f"{energy_ranges[i][1]} keV", 
               ha='left', va='top', fontsize=20, transform=a.transAxes)
        a.set_xlabel("")

    ax.reshape(-1)[len(energy_ranges)-1].set_xlabel("Phase", fontsize=20)
    plt.subplots_adjust(hspace=.3, bottom=.08, top=.98, right=.98)

    plt.show()



