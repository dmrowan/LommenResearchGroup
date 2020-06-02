#!/usr/bin/env python

import argparse
from astropy import log
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from spectraplots import plotparams
import profile_utils
import multispectra

"""
isolate_errorbars.py

Purpose: Determine if a spectra has bleeding error bars and try 
         to deduce the cause

"""

#Locate bleeding bars and color them in spectra
def flag_energies(xd, ax, sigma_cutoff, color):
    
    log.info("Flagging energies")
    med = np.median(xd.data['counts_err'])
    std = np.std(xd.data['counts_err'])
    flagged_energies = []
    for i in range(len(xd.data['energy'])):
        if xd.data['counts_err'][i] <= med - sigma_cutoff*std:
            flagged_energy = xd.data['energy'][i]
            flagged_energies.append(flagged_energy)
            ax.axvline(flagged_energy, color=color, alpha=.2, lw=4)
    
    return ax

#Make histograms corresponding to the errors on counts
def error_hists(txt_list, phase_ranges, mincounts_list, colors=None):
    for f in txt_list:
        assert(os.path.isfile)

    xd_list = [ multispectra.xspecdata(f) for f in txt_list ]
    if colors is None:
        colors = ['xkcd:azure'] * len(xd_list)

    fig, ax = plt.subplots(len(xd_list), 1, figsize=(8, 4*len(xd_list)))
    for i in range(len(xd_list)):
        xd_list[i].set_phaserange(phase_ranges[i][0], phase_ranges[i][1])
        xd_list[i].set_counts(mincounts_list[i])
        a = ax.reshape(-1)[i]
        a = plotparams(a)
        a.hist(xd_list[i].data['counts_err'], bins=20, 
               color=colors[i], edgecolor='black')
        med = np.median(xd_list[i].data['counts_err'])
        std = np.std(xd_list[i].data['counts_err'])
        for s in [1, 2, 3]:
            a.axvline(med+s*std, color='gray', ls=':', lw=4)
            a.axvline(med-s*std, color='gray', ls=':', lw=4)

        a.set_xlabel("Error on counts/sec", fontsize=20)
        a.set_ylabel("Number of energy bins", fontsize=20)

    plt.subplots_adjust(top=.98, right=.98, bottom=.1, left=.1)
    fig.savefig("CountsErrors.pdf", dpi=500)

