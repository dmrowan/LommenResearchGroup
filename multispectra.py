#!/usr/bin/env python
import spectraplots
import collections
import genspectra
import pexpect
import time
import os
import pandas as pd
from LCClass import LightCurve
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from fuzzywuzzy import process
from spectraplots import plotparams
import numpy as np
from tqdm import tqdm
from astropy import log

#Dom Rowan 2019

desc="""
Generate multiple spectra with different phase selections

Wrapper function (equiv to main): 
    (1) Finds phase regions 
    (2) Generates spectra 
    (3) Makes a plot

    Arguments:
        param evt : input event file 
        type evt: str

        param lower : lower phase boundary of OFF_PULSE
        type lower: float

        param upper: upper phase boundary of OFF_PULSE
        type upper: float

        param mincounts: default mininum counts used for grppha for primary
        type mincounts: int

        param mincounts_interpulse: default minimum counts used for grppha
                                    for interpulse
        type mincounts_interpulse: int

        param use_variable_mincounts: change the grppha mincounts dependending                                      on width of phase selection
        type use_variable_mincounts: bool (default False)

        param scalefactor: How much to adjust grppha mincounts by if using 
                           variable mincounts
        type scalefactor: float (default 1.00)

        param output: output plot filename
        type output: str (default "multispectra.pdf")
"""

#Fix phasetups that overlap zero
def phase_correction(phasetup):
    if phasetup[1] <= phasetup[0]:
        phasetup = (phasetup[0], phasetup[1]+1)
    return phasetup

def gen_multispectra(evt, onpeak_ranges, interpulse_ranges, 
                     offpeak_range, mincounts, mincounts_interpulse):
    assert(os.path.isfile("autofitting.xcm"))
    clobber=True
    print("*****Generating Off-Peak Spectra*****")
    genspectra.gen_spectra(evt, 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           mincounts[0], 
                           save_pha="offpeak.pha", run_xspec=False)

    print("*****Generating On-peak spectra*****")
    #Define phasebins kinda manually, iterate through them
    for i, tup in enumerate(tqdm(onpeak_ranges)):
        #This generates the pha file we normally open in xspec manually
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts[i], 
                               save_pha=f"onpeak_{i}.pha", run_xspec=False,
                               verbose=False)

        #This is opening xspec in python and doing basic fitting 
        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 onpeak_{i}.pha")
        xspec.expect("XSPEC12>")
        #This is an xspec script that loads data and fits model
        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")
        #Save the model params with this command
        xspec.sendline(f"save model onpeak_model_{i}")
        if os.path.isfile(f"onpeak_model_{i}.xcm") and clobber:
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
        xspec.sendline(f"wdata data_onpeak_{i}.txt")
        if os.path.isfile(f"data_onpeak_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")
        xspec.sendline("exit")


    print("*****Generating interpulse spectra*****")
    for i, tup in enumerate(tqdm(interpulse_ranges)):
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts_interpulse[i], 
                               save_pha=f"interpulse_{i}.pha", 
                               run_xspec=False, verbose=False)

        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 interpulse_{i}.pha")
        xspec.expect("XSPEC12>")
        #This is an xspec script that loads data and fits model
        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")
        #Save the model params with this command
        xspec.sendline(f"save model interpulse_model_{i}")
        if os.path.isfile(f"interpulse_model_{i}.xcm") and clobber:
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
        xspec.sendline(f"wdata data_interpulse_{i}.txt")
        if os.path.isfile(f"data_interpulse_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")
        xspec.sendline("exit")


#This class reads in the txt file output
class xspecdata:
    def __init__(self, filename):
        #Read in file and find where table breaks
        assert(os.path.isfile(filename))
        with open(filename) as h:
            lines = h.readlines()
        breakidx = np.where(np.array(lines) == 'NO NO NO NO NO\n')[0][0]
        #First table is the spectra
        df0 = pd.read_csv(filename, skiprows=3, delimiter=" ", header=None,
                          nrows=breakidx-3)
        #Second table gives delchi
        df1 = pd.read_csv(filename, skiprows=breakidx+1, 
                          delimiter=" ", header=None)
        df0.columns = ['energy', 'energy_err', 
                       'counts', 'counts_err', 'model']
        df1.columns = ['energy', 'energy_err', 
                       'delchi', 'delchi_err', 'model']
        self.data = df0
        self.residuals = df1

        self.width = None
        self.lower = None
        self.upper = None
        self.component = None
    def set_component(self, component):
        self.component = process.extract(component, ['primary', 'interpulse'],
                                         limit=1)[0][0]
    def set_phaserange(self, p1, p2):
        self.width=p2-p1
        self.lower = p1
        self.upper = p2
    def get_label(self):
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            return f"Phase: {self.lower} -- {self.upper}"

#Plotting routine (still needs comments)
def plot_multi_ufspec(sourcename, primarytxts, interpulsetxts, 
                      p_ranges, i_ranges, output="multispectra.pdf"):

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
    

    primary_data = []
    interpulse_data = []
    for i in range(len(primarytxts)):
        xd = xspecdata(primarytxts[i])
        xd.set_phaserange(p_ranges[i][0], p_ranges[i][1])
        primary_data.append(xd)
    for i in range(len(interpulsetxts)):
        xd = xspecdata(interpulsetxts[i])
        xd.set_phaserange(i_ranges[i][0], i_ranges[i][1])
        interpulse_data.append(xd)
    alldata = [primary_data, interpulse_data]

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR_B1821-24', 'PSR_B1937+21'],
                                 limit=1)[0][0]
    assert(sourcename in ['PSR_B1821-24', 'PSR_B1937+21'])


    labels=["Primary Pulse", "Interpulse"]


    #Plot data
    colors = ["#0062b2",
            "#ff1e33",
            "#00d5d3",
            "#d131e7"]
    labels=["Primary Pulse", "Interpulse"]
    for i, ax in enumerate([axp0, axi0]):
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
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.text(.95, .95, labels[i], transform=ax.transAxes, 
                fontsize=20, ha='right', va='top')
        ax.legend(loc=(.40, 0.05), fontsize=13, edgecolor='black')
        ax.set_xlim(right=10)
        fig.add_subplot(ax)

    #Plot residuals
    for i, ax in enumerate([axp1, axi1]):
        for j in range(len(alldata[i])):
            ax.errorbar(
                    alldata[i][j].residuals['energy'].astype(float), 
                    alldata[i][j].residuals['delchi'].astype(float),
                    xerr=alldata[i][j].residuals['energy_err'].astype(float), 
                    yerr=alldata[i][j].residuals['delchi_err'].astype(float),
                    ls=' ', marker='.', color=colors[j], alpha=0.8,
                    zorder=i)
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        fig.add_subplot(ax)


    plt.setp(axi0.get_xticklabels(), visible=False)
    plt.setp(axp0.get_xticklabels(), visible=False)

    fig.text(.03, .55, "Normalized Cts/S", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axi1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig(output, dpi=300)

#Find unique phase regions corresponding to different values of sigma
def find_phase_ranges(evt, lower, upper, nranges=4):
    lc = LightCurve(evt)
    lc.generate()
    tup = lc.peak_cutoff(lower, upper, 2)
    primary_ranges = [phase_correction((tup.min_phase, tup.max_phase))]
    interpulse_ranges = [phase_correction((tup.min_phase_ip, tup.max_phase_ip))]
    sigma = [2]
    for n in range(3, 10):
        if len(sigma)==nranges:
            break
        tup = lc.peak_cutoff(lower, upper, nsigma=n)
        primary = phase_correction((tup.min_phase, tup.max_phase))
        interpulse = phase_correction((tup.min_phase_ip, tup.max_phase_ip))
        if (primary_ranges[-1] == primary or 
                interpulse_ranges[-1] == interpulse):
            continue
        else:
            primary_ranges.append(primary)
            interpulse_ranges.append(interpulse)
            sigma.append(n)

    rangetup = collections.namedtuple('rangetup', 
            ['primary', 'interpulse', 'sigma'])
    tup = rangetup(primary_ranges, interpulse_ranges, sigma)
    return tup

#Adjust the minimum counts depending on the width of selection
def variable_mincounts(primary_ranges, interpulse_ranges, 
        mincounts, mincounts_interpulse, scalefactor=1.00):

    mincounts_both = [np.zeros(len(primary_ranges), dtype=int), 
                      np.zeros(len(interpulse_ranges), dtype=int)]
    initial_values = [mincounts, mincounts_interpulse]
    for i, ranges in enumerate([primary_ranges, interpulse_ranges]):
        width = np.array([ r[1] - r[0] for r in ranges ])
        max_idx = np.where(width == max(width))[0][0]
        init_width = width[max_idx]
        for j in range(len(width)):
            ratio = width[j] / init_width
            if ratio != 1:
                new_value = int(round(ratio*initial_values[i]*scalefactor))
                if new_value > initial_values[i]:
                    new_value = initial_values[i]
                mincounts_both[i][j] = new_value
            else:
                mincounts_both[i][j] = int(round(initial_values[i]))
    
    mincounts_tup = collections.namedtuple('mincounts_tup',
            ['primary', 'interpulse', 'scalefactor'])
    tup = mincounts_tup(mincounts_both[0], mincounts_both[1], scalefactor)
 
    return tup


#Finds phase regions, generates spectra, and plots
#lower and upper are for offpeak region
def wrapper(evt, lower, upper, mincounts, mincounts_interpulse, 
            use_variable_mincounts=False, scalefactor=1.00, 
            output="multispectra.pdf"):

    log.info(f"Finding phase ranges for {evt}")
    rangetup = find_phase_ranges(evt, lower, upper)
    if use_variable_mincounts:
        log.info(f"Using variable mincounts with initial {mincounts} "\
                  f"for primary and {mincounts_interpulse} for interpulse")
        mincounts_tup = variable_mincounts(
                rangetup.primary, rangetup.interpulse,
                mincounts, mincounts_interpulse)
        mincounts = mincounts_tup.primary
        mincounts_interpulse = mincounts_tup.interpulse
    else:
        mincounts = [mincounts]*len(rangetup.primary)
        mincounts_interpulse = [mincounts_interpulse]*len(
                rangetup.interpulse)

    gen_multispectra(evt, rangetup.primary, rangetup.interpulse, 
                     (lower, upper), mincounts, mincounts_interpulse)

    primarytxts = [f"data_onpeak_{i}.txt" for i in range(4)]
    interpulsetxts = [f"data_interpulse_{i}.txt" for i in range(4)]
    plot_multi_ufspec(evt, primarytxts, interpulsetxts,
                 rangetup.primary, rangetup.interpulse, output=output)



if __name__ == '__main__':
    wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200, 800, 
            use_variable_mincounts=True, scalefactor=1.25, 
            output="1821multispectra.pdf")
    #wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 400)

