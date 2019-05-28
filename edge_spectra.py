#!/usr/bin/env python

import argparse
from astropy import log
import collections
from fuzzywuzzy import process
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
import pexpect
import time
from tqdm import tqdm

import genspectra
from LCClass import LightCurve
import multispectra
import profile_utils
from spectraplots import plotparams

#Dom Rowan 2019

desc="""
Similar to multispectra, produce plots of unfolded spectra
fit with power law model for leading and trailing edges
of primary peak component
"""

def gen_edge_spectra(evt, off1, off2,
                     leading_ranges, trailing_ranges,
                     leading_mincounts, trailing_mincounts,
                     lower_energy=1.0):
    
    assert(os.path.isfile("autofitting.xcm"))
    assert(os.path.isfile("runsetup.xcm"))
    
    clobber=True

    log.info("Generating Off-Peak Spectra")
    genspectra.gen_spectra(evt, off1, off2, 0, 1200,
                           leading_mincounts[0], save_pha="offpeak.pha", run_xspec=False)

    log.info("Generating Leading Edge Spectra")
    for i, tup in enumerate(tqdm(leading_ranges)):
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, leading_mincounts[i],
                               save_pha=f"leading_{i}.pha", run_xspec=False,
                               verbose=False)
        xspec = pexpect.spawn("xspec")
        #This is opening xspec in python and doing basic fitting 
        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 leading_{i}.pha")
        xspec.expect("XSPEC12>")
        #This is an xspec script that loads data and fits model
        xspec.sendline("@runsetup.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline(f"ig **-{lower_energy}, 10.-**")
        xspec.expect("XSPEC12>")
        
        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")

        #Save the model params with this command
        xspec.sendline(f"save model leading_model_{i}")
        if os.path.isfile(f"leading_model_{i}.xcm") and clobber:
            xspec.sendline("y")
        xspec.expect("XSPEC12>")
        #Set the xaxis to be keV for our plots (and txt files)
        xspec.sendline("setplot energy")
        xspec.expect("XSPEC12>")
        #Plot both at same time then save txt file
        xspec.sendline(f"plot ufspec delchi")
        time.sleep(3)
        xspec.expect("XSPEC12>")
        xspec.sendline("ipl")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"wdata data_leading_{i}.txt")
        if os.path.isfile(f"data_leading_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")
        xspec.sendline("exit")

    log.info("Generating Trailing Edge Spectra")
    for i, tup in enumerate(tqdm(trailing_ranges)):
        genspectra.gen_spectra(evt, tup[0], tup[1], 
                               0, 1200, trailing_mincounts[i],
                               save_pha=f"trailing_{i}.pha", run_xspec=False,
                               verbose=False)

        xspec = pexpect.spawn("xspec")
        #This is opening xspec in python and doing basic fitting 
        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 trailing_{i}.pha")
        xspec.expect("XSPEC12>")
        #This is an xspec script that loads data and fits model
        xspec.sendline("@runsetup.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline(f"ig **-{lower_energy}, 10.-**")
        xspec.expect("XSPEC12>")
        
        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")

        #Save the model params with this command
        xspec.sendline(f"save model trailing_model_{i}")
        if os.path.isfile(f"trailing_model_{i}.xcm") and clobber:
            xspec.sendline("y")
        xspec.expect("XSPEC12>")
        #Set the xaxis to be keV for our plots (and txt files)
        xspec.sendline("setplot energy")
        xspec.expect("XSPEC12>")
        #Plot both at same time then save txt file
        xspec.sendline(f"plot ufspec delchi")
        time.sleep(3)
        xspec.expect("XSPEC12>")
        xspec.sendline("ipl")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"wdata data_trailing_{i}.txt")
        if os.path.isfile(f"data_trailing_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")
        xspec.sendline("exit")

def edge_phase_ranges(evt, off1, off2, nsigma=2, nranges=4):
    edge_tup = profile_utils.find_edge(evt, off1, off2, nsigma=2)
    
    leading_ranges = [ (edge_tup.min, edge_tup.peak) ]
    while len(leading_ranges) < nranges:
        if leading_ranges[-1][0] + 0.01 == edge_tup.peak:
            break
        else:
            leading_ranges.append( (leading_ranges[-1][0]+0.01, edge_tup.peak) )

    trailing_ranges = [ (edge_tup.peak, edge_tup.max) ]
    while len(trailing_ranges) < nranges:
        if trailing_ranges[-1][0] - 0.01 == edge_tup.peak:
            break
        else:
            trailing_ranges.append( (edge_tup.peak, trailing_ranges[-1][1]-0.01) )

    return (leading_ranges, trailing_ranges)

#Plotting routine 
def plot_multi_edge_spec(sourcename, leadingtxts, trailingtxts, 
                         leading_ranges, trailing_ranges, 
                         leading_mincounts, trailing_mincounts,
                         output="edgespectra.pdf",
                         leading_models=[],
                         trailing_models=[]):

    log.info("Plotting Spectra")
    #Init figure
    fig = plt.figure(figsize=(10, 11))
    plt.subplots_adjust(top=.98, right=.98, hspace=.15, left=.15)
    #Outer gridspec of size 2
    outer = gridspec.GridSpec(2, 1, height_ratios=[1,1])

    #Each spec of outer contains ufspec and delchi
    inner_l = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                               hspace=0, height_ratios=[3, 1])
    inner_t = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],
                                               hspace=0, height_ratios=[3,1])

    #Make axes we can plot onto
    axl1 = plt.Subplot(fig, inner_l[1])
    axl0 = plt.Subplot(fig, inner_l[0], sharex=axl1)
    axt1 = plt.Subplot(fig, inner_t[1])
    axt0 = plt.Subplot(fig, inner_t[0], sharex=axt1)
    
    #Fill lists of xspecdata objects
    leading_data = []
    trailing_data = []

    use_counts_label=False #Setting this manually for now 
    for i in range(len(leadingtxts)):
        xd = multispectra.xspecdata(leadingtxts[i])
        xd.set_phaserange(round(leading_ranges[i][0], 2), 
                          round(leading_ranges[i][1], 2))
        if use_counts_label:
            xd.set_counts(leading_mincounts[i])
        if len(leading_models) != 0:
            xd.phot_index_from_file(leading_models[i])
        leading_data.append(xd)
    for i in range(len(trailingtxts)):
        xd = multispectra.xspecdata(trailingtxts[i])
        xd.set_phaserange(round(trailing_ranges[i][0],2), 
                          round(trailing_ranges[i][1],2))
        if use_counts_label:
            xd.set_counts(trailing_mincounts[i])
        if len(trailing_models) != 0:
            xd.phot_index_from_file(trailing_models[i])
        trailing_data.append(xd)

    #Make one list with both to easily iterate through
    alldata = [leading_data, trailing_data]

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR_B1821-24', 'PSR_B1937+21'],
                                 limit=1)[0][0]
    assert(sourcename in ['PSR_B1821-24', 'PSR_B1937+21'])

    #Labels for each plot
    labels=["Leading Edge", "Trailing Edge"]

    #Plot data
    colors = ["#d5483a",
            "#70c84c",
            "#853bce",
            "#d4ae2f",
            "#625cce",
            "#c24ebe"]

    #Iterate through xspecdata and axes
    for i, ax in enumerate([axl0, axt0]):
        for j in range(len(alldata[i])):
            ax.errorbar(alldata[i][j].data['energy'], 
                        alldata[i][j].data['counts'],
                        xerr = alldata[i][j].data['energy_err'],
                        yerr = alldata[i][j].data['counts_err'],
                        ls=' ', marker='o', color=colors[j],
                        label=alldata[i][j].get_label(),
                        zorder=i)
            ax.plot(alldata[i][j].data['energy'],
                    alldata[i][j].data['model'],
                    ls='-', lw=3, color=colors[j],
                    zorder=len(alldata[i])+i, 
                    label='_nolegend_')
            #Testing isolation of errorbars
            #ax = isolate_errorbars.flag_energies(alldata[i][j], ax, 2, colors[j]) 
        #Set plot parameters
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        if sourcename == 'PSR_B1837+21':
            ax.text(.05, .95, labels[i], transform=ax.transAxes, 
                    fontsize=20, ha='left', va='top')
        else:
            ax.text(.95, .95, labels[i], transform=ax.transAxes, 
                    fontsize=20, ha='right', va='top')
        ax.legend(loc=(.30, 0.05), fontsize=13, edgecolor='black')
        ax.set_xlim(right=10)
        fig.add_subplot(ax)

    #Plot residuals
    for i, ax in enumerate([axl1, axt1]):
        for j in range(len(alldata[i])):
            ax.errorbar(
                    alldata[i][j].residuals['energy'].astype(float), 
                    alldata[i][j].residuals['delchi'].astype(float),
                    xerr=alldata[i][j].residuals['energy_err'].astype(float), 
                    yerr=alldata[i][j].residuals['delchi_err'].astype(float),
                    ls=' ', marker='.', color=colors[j], alpha=0.8,
                    zorder=i)
        ax = plotparams(ax)
        ax.axhline(0, ls=':', lw=1.5, color='gray')
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        ax.set_ylabel(r'Residuals ($\sigma$)', fontsize=15)
        fig.add_subplot(ax)

    #Dont want to show xtick labels for ufspec
    plt.setp(axl0.get_xticklabels(), visible=False)
    plt.setp(axt0.get_xticklabels(), visible=False)

    #Add axes labels
    fig.text(.03, .55, "Normalized Cts/S", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axt1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig(output, dpi=300)

def main(evt, off1, off2, lower_energy, mcl, mct):
    log.info("Finding leading and trailing phase ranges")
    leading_ranges, trailing_ranges = edge_phase_ranges(evt, off1, off2)
    
    log.info("Using variable mincounts")
    mincounts_tup = multispectra.variable_mincounts(
            leading_ranges, trailing_ranges,
            mcl, mct, scalefactor=1.25)

    mincounts_leading = mincounts_tup.primary
    mincounts_trailing = mincounts_tup.interpulse
    log.info("Generating Spectra")
    #gen_edge_spectra(evt, off1, off2, 
    #                 leading_ranges, trailing_ranges, 
    #                 mincounts_leading, mincounts_trailing,
    #                 lower_energy=lower_energy)

    leadingtxts = [f"data_leading_{i}.txt" for i in range(len(leading_ranges))]
    trailingtxts = [f"data_trailing_{i}.txt" for i in range(len(trailing_ranges))]
    leading_models = [f"leading_model_{i}.xcm" for i in range(len(leading_ranges))]
    trailing_models = [f"trailing_model_{i}.xcm" for i in range(len(trailing_ranges))]

    plot_multi_edge_spec(evt, leadingtxts, trailingtxts, 
                         leading_ranges, trailing_ranges, 
                         mincounts_leading, mincounts_trailing,
                         output="edgespectra.pdf",
                         leading_models=leading_models,
                         trailing_models=trailing_models)


if __name__ == '__main__':
    main("../PSR_B1937+21_combined.evt", .2, .4, .8, 200, 800)
