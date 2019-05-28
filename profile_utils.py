#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from LCClass import LightCurve
import argparse
import os

desc="""
Various profile tools to use for generating pulse profiles
"""

# Produce energy profile over selected energy range in keV
def energy_filtered_profile(evt, energy_min, energy_max, 
                            phase_min=None, phase_max=None):
    if type(evt) != str:
        raise ValueError("filename must be string")
    for var_in in [energy_min, energy_max, phase_min, phase_max]:
        if var_in is not None:
            if type(var_in) not in [int, float]:
                raise ValueError(f"{var_in} must be int or float")
    if not os.path.isfile(evt):
        raise FileNotFoundError

    fig, ax = plt.subplots(1, 1, figsize=(8,4))

    lc = LightCurve(evt)
    pi_min = energy_min * 100
    pi_max = energy_max * 100
    lc.mask(lower_pi=pi_min, upper_pi=pi_max)
    lc.generate()

    ax = lc.plot(ax=ax)
    plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
    if phase_min is not None and phase_max is not None:
        ax.axvspan(phase_min, phase_max, color='gray', alpha=.2)
    plt.show()

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
    plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
    plt.show()


# Find primary and interpulse phase ranges with user
# defined offpeak phase ranges
def find_phase_ranges(evt, off1, off2, nsigma):

    if type(evt) != str:
        raise ValueError("filename must be string")
    if any( [type(v) not in [float, int] for v in [off1, off2, nsigma] ]):
        raise ValueError("value must be int or float")
    if not os.path.isfile(evt):
        raise FileNotFoundError

    lc = LCClass.LightCurve(evt)
    lc.generate()
    cutofftup = lc.peak_cutoff(off1, off2, nsigma=nsigma)
    print(f"Min phase primary: {cutofftup.min_phase}")
    print(f"Max phase primary: {cutofftup.max_phase}")
    print(f"Min phase interpulse: {cutofftup.min_phase_ip}")
    print(f"Max phase interpulse: {cutofftup.max_phase_ip}")
    return cutofftup

