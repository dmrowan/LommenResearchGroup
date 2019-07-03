#!/usr/bin/env python
import argparse
import ast
from astropy import log
from astropy.table import Table
import collections
from fuzzywuzzy import process
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pexpect
import subprocess
import time
from tqdm import tqdm
from LCClass import LightCurve
from spectraplots import plotparams
import isolate_errorbars
import profile_utils
import genspectra
import _pickle as pickle

#Dom Rowan 2019

desc="""
Generate multiple spectra with different phase selections and plot
"""

#Function to return correct number of significant figures given an error
def round_sigfigs(val, err):
    idx = len(str(err))
    skip1=True
    for i in range(len(str(err))):
        if str(err)[i] in ['-', '0', '.']:
            continue
        elif str(err)[i] == '1' and skip1:
            skip1=False
            continue
        else:
            idx=i
            break
    err_rounded = round(err, idx-1)
    val_rounded = round(val, idx-1)
    return val_rounded, err_rounded

#Fix phasetups that overlap zero
def phase_correction(phasetup):
    if type(phasetup) != tuple:
        raise TypeError
    if phasetup[1] <= phasetup[0]:
        phasetup = (phasetup[0], phasetup[1]+1)
    return phasetup

#Function to read through output xcm file for phot index
def parse_model(xcm):
    assert(os.path.isfile(xcm))
    with open(xcm) as f:
        lines = np.array(f.readlines())
    prefix_end = [ i for i in range(len(lines)) if lines[i].find("model") >= 0 ][0]
    phot_index = float(lines[prefix_end+1].split()[0])
    return phot_index

#Parse through pexpect xspec line to find the phot index
def parse_log_old(txt):
    """
    Information needed from log:
    Photon index
    Photon index error
    Absorbing column
    absorbing column error
    Chi2 
    DOF
    """
    with open(txt, 'r') as f:
        lines = np.array(f.readlines())

    idx_params = np.where( (np.array(lines) == ' par  comp\n') | 
                           (np.array(lines) == '# par  comp\n') )[0][-1]
    phot_index = float(lines[idx_params+1].split()[4])
    nH = float(lines[idx_params+1].split()[4])

    idx_conf = np.where( (np.array(lines) == 'XSPEC12>error 1 3\n') | 
                         (np.array(lines) == '!XSPEC12> error 1 3\n') )[0][-1]
    if "Cannot do error calc" in lines[idx_conf+1]:
        phot_err = float('NaN')
        nH_err = float('NaN')
    else:
        phot_tuple = ast.literal_eval(lines[idx_conf+2].split()[3])
        nH_tuple = ast.literal_eval(lines[idx_conf+3].split()[3])
        phot_err = np.mean( [abs(val) for val in phot_tuple] )
        nH_err = np.mean( [abs(val) for val in nH_tuple] )

    for i in np.arange(idx_params, idx_conf):
        if "Reduced" in lines[i]:
            chi2 = float(lines[i].split()[3])
            dof = int(lines[i].split()[5])

    phot_index, phot_err = round_sigfigs(phot_index, phot_err)
    nH, nH_err = round_sigfigs(nH, nH_err)
    log_tuple = collections.namedtuple(
            'log_tuple', ['phot', 'phot_err',
                          'nH', 'nH_err', 'chi2', 'dof'])

    out = log_tuple(phot_index, phot_err, nH, nH_err, chi2, dof)
    return out

#A new version of this function to work with 
# different log formats
def parse_log(txt):
    """
    Information needed from log:
    Photon index
    Photon index error
    Absorbing column
    absorbing column error
    Chi2 
    DOF
    """
    #Read txt file
    with open(txt, 'r') as f:
        lines = np.array(f.readlines())

    #Find index of params info
    idx_params = [ i for i in range(len(lines)) if 'par  comp' in lines[i] ][-1]

    phot_index = float([ lines[idx_params+1].split()[i+1] 
                         for i in range(len(lines[idx_params+1].split())) 
                         if 'PhoIndex' in lines[idx_params+1].split()[i] ][0])

    nH = float([ lines[idx_params+3].split()[i+1]
                 for i in range(len(lines[idx_params+3].split()))
                 if '10^22' in lines[idx_params+3].split()[i] ][0])


    #Find confidence interval calculation
    idx_conf = [ i for i in range(len(lines)) if 'error 1' in lines[i] ][-1]

    if "Cannot do error calc" in lines[idx_conf+1]:
        phot_err = float('NaN')
        nH_err = float('NaN')
    else:
        #Read in tuple for each component
        phot_tuple = ast.literal_eval(lines[idx_conf+2].split()[-1])

        #Have to check for warnings
        if "Warning" not in lines[idx_conf+3]:
            nH_tuple = ast.literal_eval(lines[idx_conf+3].split()[-1])
        else:
            nH_tuple = ast.literal_eval(lines[idx_conf+4].split()[-1])
        phot_err = np.mean( [abs(val) for val in phot_tuple] )
        nH_err = np.mean( [abs(val) for val in nH_tuple] )

    #Find chi2 and dof in file
    idx_chisq = [ i for i in np.arange(idx_params, idx_conf) if 'Reduced' in lines[i] ][0]
    chi2 = float([ lines[idx_chisq].split()[i+1]
                    for i in range(len(lines[idx_chisq].split()))
                    if '=' in lines[idx_chisq].split()[i] ][0])
    dof = float([ lines[idx_chisq].split()[i+1]
                  for i in range(len(lines[idx_chisq].split()))
                  if 'for' in lines[idx_chisq].split()[i] ][0])


    #Round both phot index and nH
    phot_index, phot_err = round_sigfigs(phot_index, phot_err)
    nH, nH_err = round_sigfigs(nH, nH_err)

    #Return named tuple
    log_tuple = collections.namedtuple(
            'log_tuple', ['phot', 'phot_err',
                          'nH', 'nH_err', 'chi2', 'dof'])

    out = log_tuple(phot_index, phot_err, nH, nH_err, chi2, dof)
    return out
                                        
#Make csv and latex tables
def make_table(evt, output, first_ranges, second_ranges, firstlogs, secondlogs):

    log.info("Making output table from logs")

    #Calculate number of counts for each range
    phase_ranges = first_ranges + second_ranges
    print(phase_ranges)
    evt = Table.read(evt, hdu=1)
    ncounts = []
    for r in phase_ranges:
        if not r[0] <= 1 <= r[1]:
            n = len(np.where( ( evt['PULSE_PHASE'] >= r[0] ) &
                                         ( evt['PULSE_PHASE'] <= r[1] ) )[0])

        else:
            n = len(np.where( (evt['PULSE_PHASE'] >= r[0] ) |
                                         (evt['PULSE_PHASE'] <= r[1]-1) )[0])


        ncounts.append(n)


    logs = firstlogs + secondlogs

    phot_index = []
    nH = []
    chi2 = np.zeros(len(logs))
    dof = np.zeros(len(logs))

    for i in range(len(logs)):
        tup = parse_log(logs[i])
        phot_index.append(f"{tup.phot}" + r'$\pm$' + f"{tup.phot_err}")
        nH.append(f"{tup.nH}" + r'$\pm$' + f"{tup.nH_err}")
        chi2[i] = tup.chi2
        dof[i] = tup.dof

    phase_ranges = [ f"{round(tup[0], 2)} -- {round(tup[1],2)}" for tup in phase_ranges ]

    df = pd.DataFrame({"Phase Range": phase_ranges,
                       r'$\Gamma$':phot_index, 
                       r'$N_{\rm{H}}$':nH,
                       r'$\chi^2$':chi2,
                       "Degrees of":dof,
                       "Number of":ncounts,
                       })

    df.to_csv(f"{output}.csv", index=False)
    df.to_latex(f"{output}.tex", index=False, escape=False, column_format='lccccr')
    
    with open(f"{output}.tex", 'r') as f:
        table_lines = f.readlines()
        f.close()

    table_lines.insert(3, "& &" + r'($10^{22}$ cm$^2$)' + " & & Freedom & Photons\\\ \n")
    with open(f"{output}.tex", 'w') as f:
        contents = "".join(table_lines)
        f.write(contents)
        f.close()

#Varation of V1. Starts with nsigma and increments smaller
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
    
    rangetup = collections.namedtuple('rangetup', ['primary', 'interpulse', 'nfound'])
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


def edge_phase_ranges(evt, off1, off2, nsigma=2, nranges=4):
    edge_tup = profile_utils.find_edge(evt, off1, off2, nsigma=2)
    
    leading_ranges = [ (round(edge_tup.min, 2), round(edge_tup.peak, 2)) ]
    while len(leading_ranges) < nranges:
        if leading_ranges[-1][0] + 0.01 == edge_tup.peak:
            break
        else:
            leading_ranges.append( (round(leading_ranges[-1][0]+0.01,2), round(edge_tup.peak,2)) )

    trailing_ranges = [ (round(edge_tup.peak,2), round(edge_tup.max, 2)) ]
    while len(trailing_ranges) < nranges:
        if trailing_ranges[-1][0] - 0.01 == edge_tup.peak:
            break
        else:
            trailing_ranges.append( (round(edge_tup.peak, 2), 
                                    round(trailing_ranges[-1][1]-0.01, 2)) )

    for i in range(len(trailing_ranges)):
        pair = trailing_ranges[i]
        if pair[0] >= 1.0 and pair[1] >= 1.0:
            pair = (pair[0] - 1.0, pair[1] - 1.0)
        trailing_ranges[i] = pair

    for i in range(len(leading_ranges)):
        pair = leading_ranges[i]
        if pair[0] >= 1.0 and pair[1] >= 1.0:
            pair = (pair[0] - 1.0, pair[1] - 1.0)
        leading_ranges[i] = pair

    return (leading_ranges, trailing_ranges)


def gen_multispectra(evt, first_ranges, second_ranges, offpeak_range, 
                     first_mincounts, second_mincounts, 
                     lower_energies_first, 
                     lower_energies_second):

    assert(os.path.isfile("autofitting.xcm"))
    assert(os.path.isfile("runsetup.xcm"))

    clobber=True

    log.info("Generating Off-Peak Background Spectra")

    genspectra.gen_spectra(evt, 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           first_mincounts[0], 
                           save_pha="offpeak.pha", run_xspec=False)

    log.info("Generating Top Panel Spectra")

    for i, tup in enumerate(tqdm(first_ranges)):
        genspectra.gen_spectra(evt, tup[0], tup[1], 0, 1200, 
                               first_mincounts[i],
                               save_pha=f"first_{i}.pha", run_xspec=False, 
                               verbose=False)

        #This is opening xspec in python and doing basic fitting 
        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 first_{i}.pha")
        xspec.expect("XSPEC12>")
        #This is an xspec script that loads data and fits model
        xspec.sendline("@runsetup.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline("ig bad")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"ig **-{lower_energies_first[i]}, 10.-**")
        xspec.expect("XSPEC12>")

        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline("error 1 3")
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
        xspec.sendline(f"wdata data_first_{i}.txt")

        #Save a log file
        if os.path.isfile(f"data_first_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")
        lines = []
        xspec.timeout=2
        for j in range(300):
            try:
                lines.append(xspec.readline(j).decode("utf-8"))
            except:
                break
        with open(f"log_first_{i}.txt", 'w') as handle:
            for l in lines:
                handle.write(l)

        xspec.sendline("exit")

    log.info("Generating second panel spectra")
    for i, tup in enumerate(tqdm(second_ranges)):
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, second_mincounts[i],
                               save_pha=f"second_{i}.pha", 
                               run_xspec=False, verbose=False)

        xspec = pexpect.spawn("xspec")
        xspec.expect("XSPEC12>")
        xspec.sendline(f"data 1:1 second_{i}.pha")
        xspec.expect("XSPEC12>")

        xspec.sendline("@runsetup.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline(f"ig **-{lower_energies_second[i]}, 9.-**")
        xspec.expect("XSPEC12>")

        xspec.sendline("@autofitting.xcm")
        xspec.expect("XSPEC12>")

        xspec.sendline("error 1 3")
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
        xspec.sendline(f"wdata data_second_{i}.txt")
        if os.path.isfile(f"data_second_{i}.txt") and clobber:
            xspec.sendline("yes")
        xspec.expect("XSPEC12>")

        lines = []
        xspec.timeout=2
        for j in range(300):
            try:
                lines.append(xspec.readline(j).decode("utf-8"))
            except:
                break
        with open(f"log_second_{i}.txt", 'w') as handle:
            for l in lines:
                handle.write(l)

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
        self.lower = round(p1, 2)
        self.upper = round(p2, 2)

    def set_counts(self, n):
        self.mincounts = n


    def phot_index_from_log(self, fname):
        assert(os.path.isfile(fname))
        tup = parse_log(fname)
        self.phot_index = tup.phot
        self.phot_index_err = tup.phot_err


    def get_label(self):
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            #label = f"Phase: {self.lower} -- {self.upper}"
            label = f"Phase: {self.lower}" + r'$-$' + str(self.upper)

        if self.mincounts is not None:
            label = label + f", Mincounts: {self.mincounts}"

        if self.phot_index is not None:
            label = label + r', $\Gamma=$'+'{:.2f}'.format(self.phot_index)
            if self.phot_index_err is not None:
                label += r'$\pm$' + '{:.2f}'.format(self.phot_index_err)

        return label

#Plotting routine 
def plot_multi_ufspec(sourcename, firsttxts, secondtxts, 
                      first_ranges, second_ranges, first_mincounts, second_mincounts,
                      first_logs,
                      second_logs,
                      first_label,
                      second_label,
                      output="multispectra.pdf",
                      vertical=True):

    #Init figure
    if vertical:
        fig = plt.figure(figsize=(10,11))
        plt.subplots_adjust(top=.98, right=.98, hspace=.15, left=.15)
        outer = gridspec.GridSpec(2, 1, height_ratios=[1,1])
    else:
        fig = plt.figure(figsize=(20, 6.5))
        plt.subplots_adjust(top=.98, right=.98, wspace=.08, left=.05, bottom=.17)
        outer = gridspec.GridSpec(1, 2, width_ratios=[1,1])

    #Each spec of outer contains ufspec and delchi
    inner_f = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                               hspace=0, height_ratios=[3, 1])
    inner_s= gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],
                                               hspace=0, height_ratios=[3,1])

    #Make axes we can plot onto
    axf1 = plt.Subplot(fig, inner_f[1])
    axf0 = plt.Subplot(fig, inner_f[0], sharex=axf1)
    axs1 = plt.Subplot(fig, inner_s[1])
    axs0 = plt.Subplot(fig, inner_s[0], sharex=axs1)
    
    #Fill lists of xspecdata objects
    first_data = []
    second_data = []

    use_counts_label=False #Setting this manually for now 
    for i in range(len(firsttxts)):
        xd = xspecdata(firsttxts[i])
        xd.set_phaserange(first_ranges[i][0], first_ranges[i][1])
        if use_counts_label:
            xd.set_counts(first_mincounts[i])
        if len(first_logs) != 0:
            xd.phot_index_from_log(first_logs[i])
        first_data.append(xd)
    for i in range(len(secondtxts)):
        xd = xspecdata(secondtxts[i])
        xd.set_phaserange(second_ranges[i][0], second_ranges[i][1])
        if use_counts_label:
            xd.set_counts(second_mincounts[i])
        if len(second_logs) != 0:
            xd.phot_index_from_log(second_logs[i])
        second_data.append(xd)

    #Make one list with both to easily iterate through
    alldata = [first_data, second_data]

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR B1821-24', 'PSR B1937+21'],
                                 limit=1)[0][0]
    assert(sourcename in ['PSR B1821-24', 'PSR B1937+21'])

    #Labels for each plot
    labels=[first_label, second_label]


    #Plot data
    old_colors = ["#d5483a",
                  "#70c84c",
                  "#853bce",
                  #"#d4ae2f",
                  #"#625cce",
                  #"#c24ebe", 
                  "xkcd:azure"]

    colors = [['xkcd:crimson', 'xkcd:orangered', 'xkcd:azure', 'xkcd:darkblue'],
            ['xkcd:green', 'xkcd:darkgreen', 'xkcd:violet', 'xkcd:indigo']]

    #Iterate through xspecdata and axes
    for i, ax in enumerate([axf0, axs0]):
        for j in range(len(alldata[i])):
            ax.errorbar(alldata[i][j].data['energy'], 
                        alldata[i][j].data['counts'],
                        xerr = alldata[i][j].data['energy_err'],
                        yerr = alldata[i][j].data['counts_err'],
                        ls=' ', marker='o', color=colors[i][j],
                        label=alldata[i][j].get_label(),
                        zorder=i)
            ax.plot(alldata[i][j].data['energy'],
                    alldata[i][j].data['model'],
                    ls='-', lw=3, color=colors[i][j],
                    zorder=len(alldata[i])+i, 
                    label='_nolegend_')
            #Testing isolation of errorbars
            #ax = isolate_errorbars.flag_energies(alldata[i][j], ax, 2, colors[j]) 
        #Set plot parameters
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_yscale('log')

        if sourcename == 'PSR B1937+21':
            if i==0:
                ax.set_ylim(top=ax.get_ylim()[1]*2)

                ax.set_ylim(bottom=ax.get_ylim()[0]*.5)

        ax.text(.95, .95, sourcename, transform=ax.transAxes, 
                ha='right', va='top', fontsize=20)
        ax.text(.95, .85, labels[i], transform=ax.transAxes, 
                fontsize=15, ha='right', va='top')

        ax.legend(loc=(.20, 0.05), fontsize=13, edgecolor='black', framealpha=.9)
        ax.set_xlim(right=11)
        fig.add_subplot(ax)

    if sourcename == 'PSR B1937+21':
        axs0.set_ylim(top=axs0.get_ylim()[1]*2)
        axf0.set_ylim(top=axf0.get_ylim()[1]*2)

    #Plot residuals
    for i, ax in enumerate([axf1, axs1]):
        for j in range(len(alldata[i])):
            ax.errorbar(
                    alldata[i][j].residuals['energy'].astype(float), 
                    alldata[i][j].residuals['delchi'].astype(float),
                    xerr=alldata[i][j].residuals['energy_err'].astype(float), 
                    yerr=alldata[i][j].residuals['delchi_err'].astype(float),
                    ls=' ', marker='.', color=colors[i][j], alpha=0.8,
                    zorder=i)
        ax = plotparams(ax)
        ax.axhline(0, ls=':', lw=1.5, color='gray')
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        if vertical:
            ax.set_ylabel(r'Residuals ($\chi$)', fontsize=15)
        else:
            ax.set_xlabel("Energy (keV)", fontsize=30)
        fig.add_subplot(ax)

    #Dont want to show xtick labels for ufspec
    plt.setp(axs0.get_xticklabels(), visible=False)
    plt.setp(axf0.get_xticklabels(), visible=False)

    #Add axes labels
    if vertical:
        axs1.set_xlabel("Energy (keV)", fontsize=30)
        fig.text(.03, .55, "Counts/Sec", ha='center', va='center', 
                 rotation='vertical', fontsize=30)
    else:
        axf0.set_ylabel("Counts/sec", fontsize=30)
        axf1.set_ylabel(r'Residuals ($\sigma$)', fontsize=18)
    fig.savefig(output, dpi=2000)


def wrapper(evt, lower_back, upper_back, 
            first_label, second_label,
            output,
            mincounts_scalefactor=1.0,
            vertical=True, generate_new=True):

        source = process.extract(evt, ['PSR B1821-24', 'PSR B1937+21'],
                                 limit=1)[0][0]
        if source == 'PSR B1821-24':
            lower_energy_primary = 0.8
            lower_energy_interpulse = 0.8
            lower_energy_leading = 0.9
            lower_energy_trailing = 0.7
            mincounts_primary_init = 1200
            mincounts_interpulse_init = 800
            mincounts_leading_init = 200
            mincounts_trailing_init = 600
        else:
            lower_energy_primary = 1.0
            lower_energy_interpulse = 1.0
            lower_energy_leading = 0.7
            lower_energy_trailing = 0.8
            mincounts_primary_init = 1500
            mincounts_interpulse_init = 800
            mincounts_leading_init = 200
            mincounts_trailing_init = 800



        if os.path.isfile("ranges.pickle"):
            log.info(f"Loading phase ranges from file")
            allranges = pickle.load( open("ranges.pickle", 'rb') )
            primary_ranges = allranges[0]
            interpulse_ranges = allranges[1]
            leading_ranges = allranges[2]
            trailing_ranges = allranges[3]
        else:

            log.info(f"Finding phase ranges for {evt}")


            pi_rangetup = find_phase_ranges_v2(evt, lower_back, upper_back, nranges=2)
            primary_ranges = pi_rangetup.primary
            interpulse_ranges = pi_rangetup.interpulse

            leading_ranges, trailing_ranges = edge_phase_ranges(
                    evt, lower_back, upper_back, nranges=2)

            pickle.dump([primary_ranges, interpulse_ranges, 
                        leading_ranges, trailing_ranges], 
                        open("ranges.pickle", 'wb'))

        pi_mincounts_tup = variable_mincounts(primary_ranges, interpulse_ranges, 
                                              mincounts_primary_init, mincounts_interpulse_init,
                                              scalefactor=1.25)


        lt_mincounts_tup = variable_mincounts(leading_ranges, trailing_ranges,
                                              mincounts_leading_init, mincounts_trailing_init,
                                              scalefactor=1.25)

        first_ranges = primary_ranges + interpulse_ranges
        lower_energies_first = [lower_energy_primary]*len(primary_ranges) + [lower_energy_interpulse]*len(interpulse_ranges)
        first_mincounts = np.concatenate((pi_mincounts_tup.primary, pi_mincounts_tup.interpulse))
        print(first_mincounts)

        second_ranges = leading_ranges + trailing_ranges
        lower_energies_second = [lower_energy_leading]*len(leading_ranges) + [lower_energy_trailing]*len(trailing_ranges)
        second_mincounts = np.concatenate((lt_mincounts_tup.primary, lt_mincounts_tup.interpulse))
        print(second_mincounts)
        
        log.info(f"Using variable mincounts with initial {first_mincounts} "\
                 f"and {second_mincounts}")

        if generate_new:
            log.info("Generating Spectra")

            gen_multispectra(evt, first_ranges, second_ranges, (lower_back, upper_back),
                             first_mincounts, second_mincounts,
                             lower_energies_first, lower_energies_second)

        firsttxts = [f"data_first_{i}.txt" for i in range(len(first_ranges))]
        secondtxts = [f"data_second_{i}.txt" for i in range(len(second_ranges))]
        firstlogs = [f"log_first_{i}.txt" for i in range(len(first_ranges))]
        secondlogs = [f"log_second_{i}.txt" for i in range(len(second_ranges))]

        log.info("Plotting Spectra")

        plot_multi_ufspec(evt, firsttxts, secondtxts,
                          first_ranges, second_ranges, first_mincounts, second_mincounts,
                          firstlogs, secondlogs, first_label, second_label, output, 
                          vertical=vertical)


        make_table(evt, f"Model{source.replace('PSR ', '')[:-3]}", 
                   first_ranges, second_ranges, 
                   firstlogs, secondlogs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("evt", help='Event file', type=str)
    parser.add_argument("output", help="Output image file", type=str, nargs='?', 
                        default="plot_spectra.pdf")
    parser.add_argument("--h", help="Make plot horizontal instead of vertical",
                        action='store_true', default=False)
    parser.add_argument("--p", help="Plot with existing files in dir",
                        action='store_true', default=False)
    args= parser.parse_args()


    wrapper(args.evt, (.2, .4), (.7, .9), "Primary & Interpulse", "Leading & Trailing", 
            args.output,
            vertical=(not args.h), generate_new=(not args.p))





