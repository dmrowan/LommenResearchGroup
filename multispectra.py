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
import subprocess

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
    log.info("Generating Off-Peak Spectra")
    genspectra.gen_spectra(evt, 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           mincounts[0], 
                           save_pha="offpeak.pha", run_xspec=False)

    log.info("Generating On-Peak Spectra")
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


    log.info("Generating Interpulse spectra")
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
        self.mincounts = None
        self.phot_index = None
        self.filename = filename
    def set_component(self, component):
        self.component = process.extract(component, ['primary', 'interpulse'],
                                         limit=1)[0][0]
    def set_phaserange(self, p1, p2):
        self.width=p2-p1
        self.lower = p1
        self.upper = p2

    def set_counts(self, n):
        self.mincounts = n


    def phot_index_from_file(self, fname):
        assert(os.path.isfile(fname))
        self.phot_index = parse_model(fname)

    def get_label(self):
        ask_for_error = True #temp variable
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            #label = f"Phase: {self.lower} -- {self.upper}"
            label = f"Phase: {self.lower}" + r'$-$' + str(self.upper)

        if self.mincounts is not None:
            label = label + f", Mincounts: {self.mincounts}"

        if self.phot_index is not None:
            if ask_for_error:
                error = input(f"Enter error for {self.filename}:")
                error = float(error)
            label = label + f", PhotIndex: {round(self.phot_index, 2)}" + r'$\pm$' + str(error)
            
        return label


#Plotting routine (still needs comments)
def plot_multi_ufspec(sourcename, primarytxts, interpulsetxts, 
                      p_ranges, i_ranges, mincounts, mincounts_interpulse,
                      output="multispectra.pdf",
                      primary_models=[],
                      interpulse_models=[]):

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
    use_counts_label=False #Setting this manually for now
    for i in range(len(primarytxts)):
        xd = xspecdata(primarytxts[i])
        xd.set_phaserange(p_ranges[i][0], p_ranges[i][1])
        if use_counts_label:
            xd.set_counts(mincounts[i])
        if len(primary_models) != 0:
            xd.phot_index_from_file(primary_models[i])
        primary_data.append(xd)
    for i in range(len(interpulsetxts)):
        xd = xspecdata(interpulsetxts[i])
        xd.set_phaserange(i_ranges[i][0], i_ranges[i][1])
        if use_counts_label:
            xd.set_counts(mincounts_interpulse[i])
        if len(interpulse_models) != 0:
            xd.phot_index_from_file(interpulse_models[i])
        interpulse_data.append(xd)
    alldata = [primary_data, interpulse_data]

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR_B1821-24', 'PSR_B1937+21'],
                                 limit=1)[0][0]
    assert(sourcename in ['PSR_B1821-24', 'PSR_B1937+21'])


    labels=["Primary Pulse", "Interpulse"]


    #Plot data
    colors = ["#d5483a",
            "#70c84c",
            "#853bce",
            "#d4ae2f",
            "#625cce",
            "#c24ebe"]
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
        ax.legend(loc=(.30, 0.05), fontsize=13, edgecolor='black')
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
        ax.axhline(0, ls=':', lw=1.5, color='gray')
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        ax.set_ylabel(r'Residuals ($\sigma$)', fontsize=15)
        fig.add_subplot(ax)


    plt.setp(axi0.get_xticklabels(), visible=False)
    plt.setp(axp0.get_xticklabels(), visible=False)

    fig.text(.03, .55, "Normalized Cts/S", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axi1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig(output, dpi=300)

#Find unique phase regions corresponding to different values of sigma
def find_phase_ranges(evt, lower, upper, nranges=4, sigma0=2):
    assert(0 < sigma0 < 10)
    lc = LightCurve(evt)
    lc.generate()
    tup = lc.peak_cutoff(lower, upper, sigma0)
    primary_ranges = [phase_correction((tup.min_phase, tup.max_phase))]
    interpulse_ranges = [
            phase_correction((tup.min_phase_ip, tup.max_phase_ip))]
    sigma = [sigma0]
    for n in np.arange(sigma0+.25, 10, .25):
        if len(sigma)==nranges:
            break
        tup = lc.peak_cutoff(lower, upper, nsigma=n)
        if tup.min_phase_ip is None:
            continue

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

def find_phase_ranges_v2(evt, lower, upper, nranges, sigma0=2):
    assert(0 < sigma0 < 10)
    lc = LightCurve(evt)
    lc.generate()
    tup = lc.peak_cutoff(lower, upper, sigma0)
    primary_ranges = [phase_correction((tup.min_phase, tup.max_phase))]
    interpulse_ranges = [phase_correction((tup.min_phase_ip, tup.max_phase_ip))]
    while (len(primary_ranges) < nranges):
        if primary_ranges[-1][1] - primary_ranges[-1][0] > 0.02:
            newtup = (primary_ranges[-1][0]+.01, primary_ranges[-1][1]-.01)
            primary_ranges.append(newtup)
        else:
            break

    while (len(interpulse_ranges) < nranges):
        if interpulse_ranges[-1][1] - interpulse_ranges[-1][0] > 0.02:
            newtup = (interpulse_ranges[-1][0]+.01, interpulse_ranges[-1][1]-.01)
            interpulse_ranges.append(newtup)
        else:
            break
    
    rangetup = collections.namedtupe('rangetup', ['primary', 'interpulse', 'nfound'])
    tup = rangetup(primary_ranges, interpulse_ranges, 
                   min(len(primary_ranges), len(interpulse_ranges)))
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

def manual_mincounts(pulse):
    s = input(f"Enter minimum counts for {pulse} --- ")
    output = []
    while (s != ''):
        try:
            output.append(int(s))
        except:
            print("Invalid input")
        s = input(f"Enter minimum counts for {pulse} --- ")
    print("\n")
    return output

def parse_model(xcm):
    assert(os.path.isfile(xcm))
    with open(xcm) as f:
        lines = np.array(f.readlines())
    prefix_end = [ i for i in range(len(lines)) if lines[i].find("model") >= 0 ][0]
    phot_index = float(lines[prefix_end+1].split()[0])
    return phot_index
    

#Finds phase regions, generates spectra, and plots
#lower and upper are for offpeak region
def wrapper(evt, lower, upper, mincounts, mincounts_interpulse, 
            use_variable_mincounts=False, scalefactor=1.00, 
            output="multispectra.pdf", use_variable_ranges=True):

    if use_variable_mincounts and use_variable_ranges:
        log.info(f"Finding phase ranges for {evt}")
        rangetup = find_phase_ranges(evt, lower, upper)
        primary_ranges = rangetup.primary
        interpulse_ranges = rangetup.interpulse
        nfiles = len(rangetup.primary)

        log.info(f"Using variable mincounts with initial {mincounts} "\
                  f"for primary and {mincounts_interpulse} for interpulse")
        mincounts_tup = variable_mincounts(
                rangetup.primary, rangetup.interpulse,
                mincounts, mincounts_interpulse)
        mincounts = mincounts_tup.primary
        mincounts_interpulse = mincounts_tup.interpulse

    #In this case we cant use the variable_mincounts function because
    # we only have one phase range
    elif ((use_variable_mincounts) and (not use_variable_ranges)):
        log.info(f"Finding single phase range for {evt}")
        rangetup = find_phase_ranges(evt, lower, upper, 
                                     nranges=1, sigma0=3)

        log.info("Using variable mincounts. Enter manually:")
        mincounts = manual_mincounts("primary")
        mincounts_interpulse = manual_mincounts("interpulse")
        assert(len(mincounts) == len(mincounts_interpulse))
        nfiles = len(mincounts)
        primary_ranges = rangetup.primary*nfiles
        interpulse_ranges = rangetup.interpulse*nfiles

    elif not use_variable_mincounts and use_variable_ranges:
        log.info(f"Finding phase ranges for {evt}")
        rangetup = find_phase_ranges(evt, lower, upper)
        primary_ranges = rangetup.primary
        interpulse_ranges = rangetup.interpulse
        #primary_ranges = rangetup.primary[:rangetup.nfound]
        #interpulse_ranges = rangetup.interpulse[:rangetup.nfound]
        nfiles = len(rangetup.primary)
        #nfiles = rangetup.nfound

        log.info(f"Using static mincounts {mincounts}"\
                 f"and {mincounts_interpulse}")
        mincounts = [mincounts]*nfiles
        mincounts_interpulse = [mincounts_interpulse]*nfiles

    else:
        nfiles = 1
        log.info(f"Finding single phase range for {evt}")
        rangetup = find_phase_ranges(evt, lower, upper, 
                                     nranges=1, sigma0=3)
        primary_ranges = rangetup.primary
        interpulse_ranges = rangetup.interpulse
        mincounts = [mincounts]
        mincounts_interpulse = [mincounts_interpulse]
    
    print(primary_ranges, interpulse_ranges)
    #gen_multispectra(evt, primary_ranges, interpulse_ranges,
    #                 (lower, upper), mincounts, mincounts_interpulse)

    primarytxts = [f"data_onpeak_{i}.txt" for i in range(nfiles)]
    interpulsetxts = [f"data_interpulse_{i}.txt" for i in range(nfiles)]
    primary_models = [f"onpeak_model_{i}.xcm" for i in range(nfiles)]
    interpulse_models = [f"interpulse_model_{i}.xcm" for i in range(nfiles)]
    log.info("Plotting spectra")
    plot_multi_ufspec(evt, primarytxts, interpulsetxts,
                      primary_ranges, interpulse_ranges, 
                      mincounts, mincounts_interpulse, output=output,
                      primary_models=primary_models, 
                      interpulse_models=interpulse_models)



if __name__ == '__main__':
    #These 1821 calls all work fine
    wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200, 800, 
            use_variable_mincounts=True, scalefactor=1.25, 
            output="spectra_CT_RT.pdf", use_variable_ranges=True)
    """
    wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200, 800, 
            use_variable_mincounts=False, output="spectra_CF_RT.pdf", 
            use_variable_ranges=True)

    wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200, 800, 
            use_variable_mincounts=True, output="spectra_CT_RF.pdf", 
            use_variable_ranges=False)

    wrapper("../PSR_B1821-24_combined.evt", .1, .4, 1200, 800, 
            use_variable_mincounts=False, output="spectra_CF_RF.pdf",
            use_variable_ranges=False)

    #Some 1937 Calls
    wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 1000, 
            use_variable_mincounts=True, scalefactor=1.25, 
            output="spectra_CT_RT.pdf", use_variable_ranges=True)

    wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 1000, 
            use_variable_mincounts=False, output="spectra_CF_RT.pdf", 
            use_variable_ranges=True)

    wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 1000,
            use_variable_mincounts=True, output="spectra_CT_RF.pdf", 
            use_variable_ranges=False)

    wrapper("../PSR_B1937+21_combined.evt", .2, .4, 2000, 1000,
            use_variable_mincounts=False, output="spectra_CF_RF.pdf",
            use_variable_ranges=False)

    """ 
