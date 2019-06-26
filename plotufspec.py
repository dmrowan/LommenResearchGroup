#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import numpy as np
import os
import sys
import subprocess
from astropy.table import Table
import genspectra
from spectraplots import plotparams
from fuzzywuzzy import process
#Dom Rowan 2019

desc = """
Plot unfolded spectra produced with xspec included one model component
Uses both primary and interpulse text files as input
"""

def plotufspec(sourcename, primarytxt, interpulsetxt):

    #Init figure
    fig = plt.figure(figsize=(10, 11))
    plt.subplots_adjust(top=.98, right=.98, hspace=.15, left=.15)
    outer = gridspec.GridSpec(2, 1, height_ratios=[1,1])

    inner_p = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                               hspace=0, height_ratios=[3, 1])
    inner_i = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],
                                               hspace=0, height_ratios=[3,1])
    axp1 = plt.Subplot(fig, inner_p[1])
    axp0 = plt.Subplot(fig, inner_p[0], sharex=axp1)
    axi1 = plt.Subplot(fig, inner_i[1])
    axi0 = plt.Subplot(fig, inner_i[0], sharex=axi1)
    
    #Read in text files
    df_list = []
    for f in [primarytxt, interpulsetxt]:
        df = pd.read_csv(f, skiprows=3, delimiter=" ", header=None)
        assert(len(df.columns) == 5)
        df.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
        if df['counts'].dtype not in [int, float]:
            for i in range(len(df['counts'])):
                if df['counts'][i] == '0':
                    cutoff = i
                    break;
            df = df.head(cutoff)
            df = df.astype(float)
        df_list.append(df)

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR B1821-24', 'PSR B1937+21', 'PSR J0218+4232'],
                                limit=1)[0][0]


    errorbarparams = dict(ls=' ', color='#28145b')
    labels=["Primary Pulse", "Interpulse"]

    #Plot data
    for i, ax in enumerate([axp0, axi0]):
        ax.errorbar(df_list[i]['energy'], df_list[i]['counts'],
                    xerr=df_list[i]['energy_err'], 
                    yerr=df_list[i]['counts_err'],
                    **errorbarparams, marker='o', label='Data', zorder=2)
        ax = plotparams(ax)
        if sourcename == 'PSR_B1821-24':
            ax.text(.98, .95, labels[i], transform=ax.transAxes, fontsize=20, #for 1821
                    ha='right', va='top')
        else:
            ax.text(.98, .35, labels[i], transform=ax.transAxes, fontsize=20, #for 1937
                    ha='right', va='top')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        ax.plot(df_list[i]['energy'], df_list[i]['model'],
                ls='-', lw=3, color='#0da0ff', zorder=1, 
                label="Powerlaw*tbabs")
        if sourcename == 'PSR_B1821-24':
            ax.legend(loc=(.05, .75), fontsize=15, edgecolor='black') #for 1821
        else:
            ax.legend(loc=(.68, .05), fontsize=15, edgecolor='black') #for 1937
        fig.add_subplot(ax)

    #Plot residuals
    residuals = [ np.subtract(df_list[i]['counts'], df_list[i]['model'])
                  for i in range(len(df_list)) ]

    for i, ax in enumerate([axp1, axi1]):
        ax.errorbar(df_list[i]['energy'], residuals[i],
                    yerr=df_list[i]['counts_err'], ls=' ', marker='.',
                    color='#0da0ff')
        ax.axhline(0, ls=':', color='0.8', lw=4)
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        if sourcename == 'PSR_B1821-24':
            ax.set_ylim(bottom=-0.0001, top=0.0001) #for 1821
        else:
            ax.set_ylim(bottom=-.00003, top=.00003) #for 1937
        fig.add_subplot(ax)

    plt.setp(axi0.get_xticklabels(), visible=False)
    plt.setp(axp0.get_xticklabels(), visible=False)

    fig.text(.03, .55, "Normalized Flux", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axi1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig("plotunfolded.jpeg", dpi=300)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("sourcename", help='psr sourcename',
                        type=str)
    parser.add_argument("primary", help='Primary txt file', 
                        type=str)
    parser.add_argument("interpulse", help='Interpulse txt file', 
                        type=str)
    args = parser.parse_args()
    plotufspec(args.sourcename, args.primary, args.interpulse)
