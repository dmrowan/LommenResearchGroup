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

#Dom Rowan 2019

desc="""
Generate multiple spectra with different phase selections
"""

#Fix phasetups that overlap zero
def phase_correction(phasetup):
    if phasetup[1] <= phasetup[0]:
        phasetup = (phasetup[0], phasetup[1]+1)
    return phasetup

#Default fitting procedure
def autofitting(fname_head, clobber=True):
    assert(os.path.isfile(f"{fname_head}.pha"))
    #This file must exist in current directoy
    assert(os.path.isfile("autofitting.xcm"))
    xspec = pexpect.spawn("xspec")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"data 1:1 {fname_head}.pha")
    xspec.expect("XSPEC12>")
    #This is an xspec script that loads data and fits model
    xspec.sendline("@autofitting.xcm")
    xspec.expect("XSPEC12>")
    #Save the model params with this command
    xspec.sendline(f"save model model_{fname_head}")
    if os.path.isfile(f"model_{fname_head}.xcm") and clobber:
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
    xspec.sendline(f"wdata data_{fname_head}.txt")
    if os.path.isfile(f"data_{fname_head}.txt") and clobber:
        xspec.sendline("yes")
    xspec.expect("XSPEC12>")
    xspec.sendline("exit")


def gen_multispectra(evt, onpeak_ranges, interpulse_ranges, 
                     offpeak_range, mincounts, mincounts_interpulse):
    assert(os.path.isfile("autofitting.xcm"))
    print("*****Generating Off-Peak Spectra*****")
    genspectra.gen_spectra(evt, 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           mincounts, 
                           save_pha="offpeak.pha")

    print("*****Generating On-peak spectra*****")
    #Define phasebins kinda manually, iterate through them
    for i, tup in enumerate(onpeak_ranges):
        #This generates the pha file we normally open in xspec manually
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts, 
                               save_pha=f"onpeak_{i}.pha")

        #This is opening xspec in python and doing basic fitting 
        #autofitting(f"onpeak_{i}")

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
    for i, tup in enumerate(interpulse_ranges):
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts_interpulse, 
                               save_pha=f"interpulse_{i}.pha")

        #autofitting(f"interpulse_{i}")
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
        assert(os.path.isfile(filename))
        with open(filename) as h:
            lines = h.readlines()
        breakidx = np.where(np.array(lines) == 'NO NO NO NO NO\n')[0][0]
        df0 = pd.read_csv(filename, skiprows=3, delimiter=" ", header=None,
                          nrows=breakidx-3)
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
def multi_ufspec(sourcename, primarytxts, interpulsetxts, p_ranges, i_ranges):

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

    fig.text(.03, .55, "Normalized Flux", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axi1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig("plotunfolded.pdf", dpi=300)

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
        if (primary_ranges[-1] == primary and 
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

#Finds phase regions, generates spectra, and plots
#lower and upper are for offpeak region
def wrapper(evt, lower, upper, mincounts, mincounts_interpulse):
    rangetup = find_phase_ranges(evt, lower, upper)
    gen_multispectra(evt, rangetup.primary, rangetup.interpulse, 
                     (lower, upper), mincounts, mincounts_interpulse)

    primarytxts = [f"data_onpeak_{i}.txt" for i in range(4)]
    interpulsetxts = [f"data_interpulse_{i}.txt" for i in range(4)]
    multi_ufspec(evt, primarytxts, interpulsetxts,
                 rangetup.primary, rangetup.interpulse)


if __name__ == '__main__':
    #wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200)
    wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 400)

