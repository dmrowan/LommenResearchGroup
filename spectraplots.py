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
#Dom Rowan 2019

desc = """
Collection of Spectra Plotting options
"""

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

#Compare spectra between two lc segments
def comparsionplot(pha1, wi1, pha2, wi2, evt, output="SpectraCompare.pdf"):
    s_peak = pulsarspectra.Spectra(pha1)
    s_peak.set_phase_fromfile(wi1)
    
    s_offpeak = pulsarspectraSpectra(pha2)
    s_offpeak.set_phase_fromfile(wi2)

    lc = pulsarspectra.LightCurve(evt)
    lc.mask(lower_pi=50, upper_pi=200)
    lc.generate(2, bs=.01)

    fig = plt.figure(figsize=(16,6))
    mygs = gs.GridSpec(2, 2)
    plt.subplots_adjust(top=.98, right=.98, wspace=0.1, hspace=.32, left=.08)
    fig.text(.02, .5, 'Counts', fontsize=25, ha='left', va='center', 
             rotation=90)

    colors = ["#8229b8", "#33bb00", "#0069ef", "#ff4319"]
    
    axLC = plt.subplot2grid((2,2), (0,0), colspan=2, rowspan=1)
    axLC.plot(lc.phasebins_extended, lc.counts_extended, marker='o', 
              ls='-', color=colors[0])
    axLC.set_xlabel('Phase', fontsize=25, zorder=5, 
                    bbox=dict(boxstyle='square', fc='w', ec='w', alpha=.4))
    axLC.set_ylim(ymax=max(lc.counts_extended)*1.1)

    ax_peak = plt.subplot2grid((2,2), (1,0), colspan=1, rowspan=1)
    ax_peak.plot(s_peak.keV, s_peak.counts, marker='o', ls='-', 
                 color=colors[1])
    ax_peak.set_xlabel('Energy (keV)', fontsize=25)
    ax_offpeak = plt.subplot2grid((2,2), (1,1), colspan=1, rowspan=1)
    ax_offpeak.plot(s_offpeak.keV, s_offpeak.counts, marker='o', ls='-', 
                 color=colors[2])
    ax_offpeak.set_xlabel('Energy(keV)', fontsize=25)

    for a in [axLC, ax_peak, ax_offpeak]:
        a = plotparams(a)
    
    phase_ranges = [s_peak.phase_range, s_offpeak.phase_range]
    axes_list = [ax_peak, ax_offpeak]
    for i in range(len(axes_list)):
        phase_lower = phase_ranges[i][0] 
        phase_upper = phase_ranges[i][1]
        condition = ((lc.phasebins >= phase_lower) & 
                    (lc.phasebins <= phase_upper))
        idx_phaserange = np.where(condition)[0]
        countvalues = [ lc.counts[ii] for ii in idx_phaserange ]
        counts_lower = min(countvalues) - 10
        counts_upper = max(countvalues) + 50
        width = phase_upper - phase_lower
        height = counts_upper - counts_lower
        phase_lower += i
        phase_upper += i
        p = mpl.patches.Rectangle(
                (phase_lower, counts_lower), width=width, height=height, 
                edgecolor=colors[i+1], lw=2, facecolor='#E6E6E6', alpha=.5)
        axLC.add_patch(p)

        con1 = mpl.patches.ConnectionPatch(
                xyA=(phase_lower, counts_lower),
                xyB=(axes_list[i].get_xlim()[0], axes_list[i].get_ylim()[1]),
                coordsA="data", coordsB="data", axesA=axLC, 
                axesB=axes_list[i], color='black', alpha=.4, lw=2, zorder=1)
        con2 = mpl.patches.ConnectionPatch(
                xyA=(phase_upper, counts_lower),
                xyB=(axes_list[i].get_xlim()[1], axes_list[i].get_ylim()[1]),
                coordsA="data", coordsB="data", axesA=axLC, 
                axesB=axes_list[i], color='black', alpha=.4, lw=2, zorder=1)
        axLC.add_artist(con1)
        axLC.add_artist(con2)
                
            
    fig.savefig(output)


#Show 1 plot with residuals
def xspec_plot(datafile):
    #Read in xspec output with pandas
    df = pd.read_csv(datafile, skiprows=3, delimiter=" ", header=None)
    df.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
    
    #Default plot params
    errorbarparams = dict(ls=' ', marker='.', color='#58508d')
    modelplotparams = dict(ls='--', color='#ffa600', lw=4)

    fig = plt.figure(figsize=(16,8))
    gs1 = gridspec.GridSpec(2, 1, height_ratios=[2,1])
    ax0 = plt.subplot(gs1[0])
    ax1 = plt.subplot(gs1[1])

    ax0.errorbar(df['energy'], df['counts'], 
                   yerr=df['counts_err'], **errorbarparams)
    ax0.plot(df['energy'], df['model'], **modelplotparams)

    #Default plot params
    ax0 = plotparams(ax0)

    residuals = np.subtract(df['counts'], df['model'])
    ax1.errorbar(df['energy'], residuals, yerr=df['counts_err'],
                   **errorbarparams)
    ax1.axhline(0, ls='--', color='0.8', lw=4)
    ax1 = plotparams(ax1)

    fig.savefig("ResidualPlot.pdf")


#Show residuals on all three panels
def xspec_triple(datafile_on, datafile_off, datafile_sub):
    fig = plt.figure(figsize=(15, 20))
    plt.subplots_adjust(top=.98, right=.98)
    outer = gridspec.GridSpec(3, 1)

    fnames = [datafile_on, datafile_off, datafile_sub]
    errorbarparams = dict(ls=' ', marker='.', color='#58508d')
    modelplotparams = dict(ls='--', color='#ffa600', lw=4)
    labels = ["On-Pulse", "Off-Pulse", "On-Pulse Background Subtracted"]
    for i, f in enumerate(fnames):
        df = pd.read_csv(f, skiprows=3, delimiter=" ", header=None)
        df.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
        inner = gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=outer[i], hspace=0, height_ratios=[3,1])
        ax0 = plt.Subplot(fig, inner[0])
        ax1 = plt.Subplot(fig, inner[1])

        ax0.errorbar(df['energy'], df['counts'], 
                       yerr=df['counts_err'],
                       **errorbarparams)
        ax0.plot(df['energy'], df['model'],  **modelplotparams)

        ax0.set_ylim(bottom=df['counts'].min(), top=df['counts'].max())
        
        residuals = np.subtract(df['counts'], df['model'])
        ax1.errorbar(df['energy'], residuals,
                     yerr=df['counts_err'], **errorbarparams)

        ax1.axhline(0, ls='--', color='0.8', lw=4)

        ax0.set_ylabel("Normalized Flux", fontsize=20)
        ax1.set_xlabel("Channel", fontsize=20)

        ax0.text(.95, .95, labels[i], transform=ax0.transAxes, fontsize=20, 
                 ha='right', va='top')

        for a in [ax0, ax1]:
            a = plotparams(a)
            fig.add_subplot(a)


    fig.savefig("threepanelresiduals.pdf")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--pha", help="Spectra file path for spec", type=str, 
                        default=None)
    parser.add_argument("--evt", help="Event file path for lc", type=str,
                        default=None)
    args = parser.parse_args()

    xspec_triple('onpeak_data.txt', 
                       'offpeak_data.txt', '1821_subtracted.txt')
