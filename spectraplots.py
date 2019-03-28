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
def plot_xspec_three_model(datafile_on, datafile_off, datafile_sub, output):
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

    fig.savefig(output)
 
def plot_xspec_subtracted(datafile_on, datafile_off, datafile_sub,
                          output, nchan=1, mode='p',
                          datafile_sub2=None, pha_on=None, pha_off=None, 
                          breakenergy=None):
    #Initialize figure
    fig = plt.figure(figsize=(9.75, 12))
    plt.subplots_adjust(top=.98, right=.98, hspace=.15)
    #Define outer grispec with different height ratios
    outer = gridspec.GridSpec(3, 1, height_ratios=[1/2, 1/2, 1])

    #Default plot params
    errorbarparams = dict(ls=' ', color='#28145b')

    #Generate lists of files/labels
    fnames = [datafile_on, datafile_off, datafile_sub]
    #Use mode p for primary, i for interpulse
    if mode in ['p', 'P', 'primary', 'Primary']:
        labels = ["On-Pulse", "Off-Pulse", "On-Pulse Background Subtracted"]
    elif mode in ['i', 'I', 'interpulse', 'Interpulse']:
        labels = ["Interpulse", "Off-Pulse", "Interpulse Background Subtracted"]
    else:
        raise ValueError("mode must be primary or interpulse")

    if datafile_sub2 is not None:
        fnames.append(datafile_sub2)

    #Find exposures for annotation
    if pha_on is not None and pha_off is not None:
        exposures = [genspectra.get_exposure(pha_on), 
                     genspectra.get_exposure(pha_off)]
    else:
        exposures=[]

    #Read in all txt files as pandas df
    df_list = []
    for f in fnames:
        df = pd.read_csv(f, skiprows=3, delimiter=" ", header=None)
        df.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
        df_list.append(df)

    #Locate each axis
    ax_on = plt.Subplot(fig, outer[0])
    ax_off = plt.Subplot(fig, outer[1])
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[2], 
                                             hspace=0, height_ratios=[3, 1])
    ax0 = plt.Subplot(fig, inner[0])
    ax1 = plt.Subplot(fig, inner[1])

    #Iterate through each axis
    for i, ax in enumerate([ax_on, ax_off, ax0]):
        #Plot data with errorbars
        ax.errorbar(df_list[i]['energy']*nchan/100, df_list[i]['counts'],
                    yerr=df_list[i]['counts_err'],
                    **errorbarparams, marker='o', label='Data', zorder=1)
        #Default plotting params
        ax = plotparams(ax)
        #Slightly increase axis thickness (big plot)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.9)
        #Add label annotation
        ax.text(.95, .95, labels[i], transform=ax.transAxes, fontsize=20, 
                 ha='right', va='top')
        #Add subplot to figure
        fig.add_subplot(ax)

    #Plot the model for subtracted panel
    ax0.plot(df_list[2]['energy']*nchan/100, df_list[2]['model'],  
             ls='-', lw=5, color='#0da0ff', zorder=2, 
             label='Bknpower Model')

    #Compute residuals
    residuals1 = np.subtract(df_list[2]['counts'], df_list[2]['model'])

    #Plot residuals
    ax1.errorbar(df_list[2]['energy']*nchan/100, residuals1,
                 yerr=df_list[2]['counts_err'], 
                 ls=' ', marker='.', color='#0da0ff')

    #Do same for second model if it exists
    if len(df_list) == 4:
        ax0.plot(df_list[3]['energy']*nchan/100, df_list[3]['model'],
                 ls='-', lw=5, color='#480ebd', zorder=3, 
                 label='Powerlaw Model')
        residuals2 = np.subtract(df_list[3]['counts'], df_list[3]['model'])
        ax1.errorbar(df_list[3]['energy']*nchan/100, residuals2,
                     yerr=df_list[3]['counts_err'],
                     ls=' ', marker='.', color='#480ebd')

    #Add line at 0 for residual panel
    ax1.axhline(0, ls='--', color='0.8', lw=4)
    #Default plot params
    ax1 = plotparams(ax1)
    #Increase axis thickness
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.9)
    fig.add_subplot(ax1)

    #Add exposure annotation
    if len(exposures) != 0:
        ax_on.text(.95, .83, f"Exposure: {round(exposures[0])} s", 
                   transform=ax_on.transAxes, fontsize=20,
                   ha='right', va='top')
        ax_off.text(.95, .83, f"Exposure: {round(exposures[1])} s",
                    transform=ax_off.transAxes, fontsize=20,
                    ha='right', va='top')

    #Add vertical line at break energy
    if breakenergy is not None:
        for a in [ax0, ax1]:
            a.axvline(breakenergy, ls='--', color='black', label="Break Energy")
        legend_y = .50
    else:
        legend_y = .60

    if len(df_list) == 4:
        legend_y -= .12


    #Subtracted panel legend
    ax0.legend(loc=(.55, legend_y), fontsize=20, edgecolor='black')

    #Remove ticks to precent overlapping with residuals
    ax0.tick_params(labelbottom=False)

    #Axes labels
    fig.text(.03, .55, "Normalized Flux", ha='center', va='center',
             rotation='vertical', fontsize=30)
    ax1.set_xlabel("Energy (keV)", fontsize=30)

    #Save figure
    fig.savefig(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--on_data", help="On peak data path", type=str,
                        default=None)
    parser.add_argument("--off_data", help="Off peak data path", type=str,
                        default=None)
    parser.add_argument("--sub_data", help="Subtracted data path", type=str,
                        default=None)
    parser.add_argument("--sub_data2", help="Second subtracted data path",
                        type=str, default=None)
    parser.add_argument("--output", help="Output filename", type=str,
                        default=None)
    parser.add_argument("--nchan", help="nchan used for grppha",
                        default=1, type=int)
    parser.add_argument("--m", help="Primary or interpulse mode", 
                        default='p', type=str)
    parser.add_argument("--pha_on", help="Onpeak pha file",
                        default=None, type=str)
    parser.add_argument("--pha_off", help="Off peak pha file",  
                        default=None, type=str)
    parser.add_argument("--be", help="break energy for bknpower", 
                        default=None, type=float)

    args = parser.parse_args()

    plot_xspec_subtracted(args.on_data, args.off_data,
                          args.sub_data, args.output, nchan=args.nchan,
                          mode=args.m, datafile_sub2=args.sub_data2,
                          pha_on=args.pha_on,
                          pha_off=args.pha_off, breakenergy=args.be)
    """
    ####Example Function Calls####

    plot_xspec_subtracted('onpeak_data.txt', 'offpeak_data.txt',
                          'subtracted.txt', 'ModelPlot.pdf', 
                          nchan=5, mode='p', 
                          pha_on='onpeak.pha', pha_off='offpeak.pha',
                          breakenergy=2.33)

    plot_xspec_subtracted('interpulse_data.txt', 'offpeak_data.txt'
                          'interpulse_subtracted.txt', 
                          'InterpulseModel.pdf', nchan=5, mode='i', 
                          pha_on='interpulse.pha', pha_off='offpeak.pha', 
                          breakenergy=1.59 
    """
