#!/usr/bin/env python
import spectraplots
import argparse
import ast
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
import isolate_errorbars
import multispectra

#Dom Rowan 2019

#Plotting routine 
def plot_multi_telescope(txtfiles_1821, txtfiles_1937, 
                         labels_1821, labels_1937,
                         output):

    #Init figure
    fig = plt.figure(figsize=(10, 9))
    plt.subplots_adjust(top=.98, right=.98, hspace=.15, left=.15)
    #Outer gridspec of size 2
    outer = gridspec.GridSpec(2, 1, height_ratios=[1,1])

    #Each spec of outer contains ufspec and delchi
    inner_1821 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                               hspace=0, height_ratios=[3, 1])
    inner_1937 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],
                                               hspace=0, height_ratios=[3,1])

    #Make axes we can plot onto
    ax_top1 = plt.Subplot(fig, inner_1821[1])
    ax_top0 = plt.Subplot(fig, inner_1821[0], sharex=ax_top1)
    ax_bot1 = plt.Subplot(fig, inner_1937[1])
    ax_bot0 = plt.Subplot(fig, inner_1937[0], sharex=ax_bot1)
    
    #Fill lists of xspecdata objects
    data_1821 = []
    data_1937 = []

    for i in range(len(txtfiles_1821)):
        xd = multispectra.xspecdata(txtfiles_1821[i])
        xd.set_label(labels_1821[i])
        data_1821.append(xd)

    for i in range(len(txtfiles_1937)):
        xd = multispectra.xspecdata(txtfiles_1937[i])
        xd.set_label(labels_1937[i])
        data_1937.append(xd)

    alldata = [data_1821, data_1937]
    sourcenames = ['PSR B1821-24', 'PSR B1937+21']

    #Plot data
    colors = ["#d5483a",
            "#70c84c",
            "#853bce",
            #"#d4ae2f",
            #"#625cce",
            #"#c24ebe", 
            "xkcd:azure"]

    #Iterate through xspecdata and axes
    for i, ax in enumerate([ax_top0, ax_bot0]):
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

        ax.text(.95, .95, sourcenames[i], transform=ax.transAxes, 
                ha='right', va='top', fontsize=20)

        ax.legend(loc=(.05, 0.05), fontsize=16, edgecolor='black', framealpha=.9)
        fig.add_subplot(ax)

    #Plot residuals
    for i, ax in enumerate([ax_top1, ax_bot1]):
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
        ax.set_ylabel(r'Residuals ($\sigma$)', fontsize=15)
        fig.add_subplot(ax)

    #Dont want to show xtick labels for ufspec
    plt.setp(ax_top0.get_xticklabels(), visible=False)
    plt.setp(ax_bot0.get_xticklabels(), visible=False)

    ax_top0.set_ylim(ax_top0.get_ylim()[0] * .4)
    #Add axes labels
    fig.text(.03, .55, "Normalized Cts/S", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    ax_bot1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig(output, dpi=2000)

def main():
    txtfiles_1821 = [ 'PSR_B1821-24/OtherSpectra/' + s 
                       for s in [ 'data_onpeak_1.txt', 
                                  'nustar_data.txt', 
                                  'xte_data.txt' ] ]

    txtfiles_1937 = [ 'PSR_B1937+21/OtherSpectra/' + s
                      for s in [ 'data_onpeak_1.txt', 
                                 'nustar_data.txt',
                                 'xmm_data.txt' ] ]

    labels_1821 = [ r'$\it{NICER}$', r'$\it{NuSTAR}$', r'$\it{RXTE}$' ]
    labels_1937 = [ r'$\it{NICER}$', r'$\it{NuSTAR}$', r'$\it{XMM}$' ] 

    photInd_1937 = [ (1.03, .06 ), (1.2, .1), (1.1, .1) ]
    photInd_1821 = [ (1.16, .07 ), (1.42, .07), (1.2, .2) ]

    for i in range(len(labels_1821)):
        labels_1821[i] += r': $\Gamma=$'+str(photInd_1821[i][0])+r'$\pm$'+str(photInd_1821[i][1])

    for i in range(len(labels_1937)):
        labels_1937[i] += r': $\Gamma=$'+str(photInd_1937[i][0])+r'$\pm$'+str(photInd_1937[i][1])

    plot_multi_telescope(txtfiles_1821, txtfiles_1937, 
                         labels_1821, labels_1937, 
                         "TelescopeCompare.jpeg")


if __name__ == '__main__':
    main()
