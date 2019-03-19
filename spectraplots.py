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

def filter_compare(datadir='./', eventfile='combined.evt'):

    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    plt.subplots_adjust(top=.95, right=.95, hspace=0)
    colors = ["#8229b8", "#33bb00", "#0069ef", "#ff4319"]
    default_plot = dict(marker='o', ls='-', color=colors[0])
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_weight('light')
    default_text = dict(ha='left', va='center', fontsize=15, fontproperties=font)

    on_peak_lower = .5
    on_peak_upper = .6
    off_peak_lower = .75
    off_peak_upper = .85

    ############# PXEX #############
    s_pxex0 = pulsarspectra.gen_spectra(lower_e=50, upper_e=200, 
                                        lower_phase=on_peak_lower, upper_phase=on_peak_upper,
                                        datadir=datadir, eventfile=eventfile,
                                        session='autopython', show=False)
    ax[0][0].plot(s_pxex0.keV, s_pxex0.counts, **default_plot)
    ax[0][0] = plotparams(ax[0][0])
    ax[0][0].text(.5, .75, "Energy: xselect\nPhase: xselect", 
                  transform=ax[0][0].transAxes, **default_text)

    s_pxex1 = pulsarspectra.gen_spectra(lower_e=50, upper_e=200, 
                                        lower_phase=off_peak_lower, upper_phase=off_peak_upper,
                                        datadir=datadir, eventfile=eventfile,
                                        session='autopython2', show=False)
    ax[0][1].plot(s_pxex1.keV, s_pxex1.counts, **default_plot)
    ax[0][1] = plotparams(ax[0][1])
    ax[0][1].text(.5, .75, "Energy: xselect\nPhase: xselect", 
                  transform=ax[0][1].transAxes, **default_text)

    ############# PXEM #############
    s_pxem0 = pulsarspectra.gen_spectra(lower_e=0, upper_e=1200, 
                                        lower_phase=on_peak_lower, upper_phase=on_peak_upper, 
                                        datadir=datadir, eventfile=eventfile,
                                        session='autopython3', show=False)
    s_pxem0.filter_energy_pha(50, 200)
    ax[1][0].plot(s_pxem0.keV, s_pxem0.counts, **default_plot)
    ax[1][0] = plotparams(ax[1][0])
    ax[1][0].text(.5, .75, "Energy: Python\nPhase: xselect", 
                  transform=ax[1][0].transAxes, **default_text)

    s_pxem1 = pulsarspectra.gen_spectra(lower_e=0, upper_e=1200, 
                                       lower_phase=off_peak_lower, upper_phase=off_peak_upper,
                                       datadir=datadir, eventfile=eventfile,
                                       session='autopython4', show=False)
    s_pxem1.filter_energy_pha(50, 200)
    ax[1][1].plot(s_pxem1.keV, s_pxem1.counts, **default_plot)
    ax[1][1] = plotparams(ax[1][1])
    ax[1][1].text(.5, .75, "Energy: Python\nPhase: xselect", 
                  transform=ax[1][1].transAxes, **default_text)

    ############# PMEX #############
    tab1 = Table.read(eventfile, hdu=1)
    tab1 = tab1[np.where( (tab1['PULSE_PHASE'] >= on_peak_lower) & 
                          (tab1['PULSE_PHASE'] <= on_peak_upper))[0]]
    print(tab1['PULSE_PHASE'].min(), tab1['PULSE_PHASE'].max())
    tab1.write('autocompare1.fits', overwrite=True)
    s_pmex0 = pulsarspectra.gen_spectra(lower_e=50, upper_e=200,
                                        lower_phase=0, upper_phase=1,
                                        datadir=datadir, eventfile='autocompare1.fits',
                                        session='autopython5', show=False)
    ax[2][0].plot(s_pmex0.keV, s_pmex0.counts, **default_plot)
    ax[2][0] = plotparams(ax[2][0])
    ax[2][0].text(.5, .75, "Energy: xselect\nPhase: Python", 
                  transform=ax[2][0].transAxes, **default_text)

    tab2 = Table.read(eventfile, hdu=1)
    tab2 = tab2[np.where( (tab2['PULSE_PHASE'] >= off_peak_lower) & 
                          (tab2['PULSE_PHASE'] <= off_peak_upper))[0]]
    tab2.write('autocompare2.fits', overwrite=True)
    print(tab2['PULSE_PHASE'].min(), tab2['PULSE_PHASE'].max())
    s_pmex1 = pulsarspectra.gen_spectra(lower_e=50, upper_e=200, 
                                        lower_phase=0, upper_phase=1,
                                        datadir=datadir, eventfile='autocompare2.fits',
                                        session='autopython6', show=False)
    ax[2][1].plot(s_pmex1.keV, s_pmex1.counts, **default_plot)
    ax[2][1] = plotparams(ax[2][1])
    ax[2][1].text(.5, .75, "Energy: xselect\nPhase: Python", 
                  transform=ax[2][1].transAxes, **default_text)


    ############# PMEM #############
    s_pmem0 = pulsarspectra.Spectra(eventfile, evt=True)
    s_pmem0.filter_energy_and_phase_evt(50, 200, on_peak_lower, on_peak_upper)
    ax[3][0].plot(s_pmem0.keV, s_pmem0.counts, **default_plot)
    ax[3][0] = plotparams(ax[3][0])
    ax[3][0].text(.5, .75, "Energy: Python\nPhase: Python", 
                  transform=ax[3][0].transAxes, **default_text)

    s_pmem1 = pulsarspectra.Spectra(eventfile, evt=True)
    s_pmem1.filter_energy_and_phase_evt(50, 200, off_peak_lower, off_peak_upper)
    ax[3][1].plot(s_pmem1.keV, s_pmem1.counts, **default_plot)
    ax[3][1] = plotparams(ax[3][1])
    ax[3][1].text(.5, .75, "Energy: Python\nPhase: Python", 
                  transform=ax[3][1].transAxes, **default_text)


    fig.text(.95/2, .05, "Energy (keV)", fontsize=30, ha='center', va='center')
    fig.text(.05, .95/2, "Counts", rotation='vertical', ha='center', va='center', fontsize=30)
    fig.savefig("Filtercompare.pdf")

def xspecplots(datafile_on, datafile_off, datafile_sub):

    df_on = pd.read_csv(datafile_on, skiprows=3, 
                        delimiter=" ", header=None)
    df_off = pd.read_csv(datafile_off, skiprows=3, 
                         delimiter=" ", header=None)
    df_sub = pd.read_csv(datafile_sub, skiprows=3, 
                         delimiter=" ", header=None)
    df_on.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
    df_off.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
    df_sub.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']

    fig, ax = plt.subplots(3, 1, figsize=(15, 20))
    
    errorbarparams = dict(ls=' ', marker='.', color='#58508d')
    modelplotparams = dict(ls='--', color='#ffa600', lw=4)
    ax[0].errorbar(df_on['energy'], df_on['counts'], 
                   yerr=df_on['counts_err'],
                   **errorbarparams)
    ax[1].errorbar(df_off['energy'], df_off['counts'], 
                   yerr=df_off['counts_err'],
                   **errorbarparams)
    ax[2].errorbar(df_sub['energy'], df_sub['counts'],
                   yerr=df_sub['counts_err'],
                   **errorbarparams)

    ax[0].plot(df_on['energy'], df_on['model'],  **modelplotparams)
    ax[1].plot(df_off['energy'], df_off['model'], **modelplotparams)
    ax[2].plot(df_sub['energy'], df_sub['model'], **modelplotparams)

    ax[0] = plotparams(ax[0])
    ax[1] = plotparams(ax[1])
    ax[2] = plotparams(ax[2])

    ax[0].set_ylim(ymin=df_on['counts'].min(), ymax=df_on['counts'].max())
    ax[1].set_ylim(ymin=df_off['counts'].min(), ymax=df_off['counts'].max())
    ax[2].set_ylim(ymin=df_sub['counts'].min(), ymax=df_sub['counts'].max())

    fig.text(.05, 0.5, 'Normalized Counts', fontsize=30, ha='center', 
             va='center', rotation='vertical')
    fig.text(.5, .05, 'Channel', fontsize=30, ha='center', va='center', 
             rotation='horizontal')

    ax[0].text(.5, .75, 
               "On-peak, " + r'$0\leq\phi\leq0.06$' +"\nNo background subtraction",
               transform=ax[0].transAxes, fontsize=25)
    ax[1].text(.5, .75, 
               "Off-peak " + r'$0.2\leq\phi\leq0.5$' + "\n'Background' spectra'",
               transform=ax[1].transAxes, fontsize=25)
    ax[2].text(.5, .75, "Background subtracted",
               transform=ax[2].transAxes, fontsize=25)
    fig.savefig("XspecModelPlots.pdf")

def xspec_residuals(datafile):
    df = pd.read_csv(datafile, skiprows=3, delimiter=" ", header=None)
    df.columns = ['energy', 'energy_err', 'counts', 'counts_err', 'model']
    
    errorbarparams = dict(ls=' ', marker='.', color='#58508d')
    modelplotparams = dict(ls='--', color='#ffa600', lw=4)

    fig = plt.figure(figsize=(16,8))
    gs1 = gridspec.GridSpec(2, 1, height_ratios=[2,1])
    ax0 = plt.subplot(gs1[0])
    ax1 = plt.subplot(gs1[1])

    ax0.errorbar(df['energy'], df['counts'], 
                   yerr=df['counts_err'], **errorbarparams)
    ax0.plot(df['energy'], df['model'], **modelplotparams)

    ax0 = plotparams(ax0)

    residuals = np.subtract(df['counts'], df['model'])
    ax1.errorbar(df['energy'], residuals, yerr=df['counts_err'],
                   **errorbarparams)
    ax1.axhline(0, ls='--', color='0.8', lw=4)
    ax1 = plotparams(ax1)

    fig.savefig("ResidualPlot.pdf")


def xspec_allresiduals(datafile_on, datafile_off, datafile_sub):
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

    """
    comparsionplot('spec_onpeak.pha', 'phase_onpeak.wi',
                   'spec_offpeak.pha', 'phase_offpeak.wi',
                   '1821_combined.evt')
    """
    #filter_compare()
    #xselectplots('onpeak_data.txt', 'offpeak_data.txt', 'subtracted_data.txt')
    #xselect_residuals("subtracted_data.txt")
    #xselect_allresiduals('onpeak_data.txt', 'offpeak_data.txt', 'subtracted_data.txt')
    xspec_allresiduals('onpeak_data.txt', 'offpeak_data.txt', '1821_subtracted.txt')
