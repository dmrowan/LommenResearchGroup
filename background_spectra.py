#!/usr/bin/env python

import argparse
import glob
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import ConnectionPatch, Rectangle
import mpl_toolkits.axes_grid1.inset_locator as il
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd
import pexpect
import subprocess
import time

import genspectra
import niutils
import xspeclog

#Dom Rowan 2020

desc="""
XSPEC analysis of BKGD_RXTE data from NICER
"""

rc('text', usetex=True)

def xspec_routine(pha):
    assert(os.path.isfile(pha))

    datafile = pha.replace('.pha', '_data.txt')

    xspec = pexpect.spawn("xspec")
    xspec.expect('XSPEC12>')

    xspec.sendline(f'data 1 {pha}')
    xspec.expect('XSPEC12>')

    xspec.sendline('cpd /xs')
    xspec.expect('XSPEC12>')

    xspec.sendline('ig **-.2, 10.-**')
    xspec.expect('XSPEC12>')

    xspec.sendline('setplot energy')
    xspec.expect('XSPEC12>')

    xspec.sendline('plot data')
    xspec.expect('XSPEC12>')

    time.sleep(3)

    xspec.sendline('ipl')
    xspec.expect('PLT>')

    xspec.sendline(f"wdata {datafile}")
    if os.path.isfile(datafile):
        xspec.expect('exists, reuse it?')
        xspec.sendline("yes")
    
    xspec.expect('PLT>')
    xspec.sendline('exit')


def xspec_environ_spectra(bkg):
    assert(os.path.isfile(bkg))

    datafile = bkg.replace('_bkg.pha', '_bkg_data.txt')

    xspec = pexpect.spawn("xspec")
    xspec.expect('XSPEC12>')

    xspec.sendline(f'data 1 {bkg}')
    xspec.expect('XSPEC12>')

    xspec.sendline('@environ_script.xcm')
    xspec.expect('XSPEC12>')

    xspec.sendline('ipl')
    xspec.expect('PLT>')

    xspec.sendline(f"wdata {datafile}")
    if os.path.isfile(datafile):
        xspec.expect('exists, reuse it?')
        xspec.sendline("yes")
    
    xspec.expect('PLT>')
    xspec.sendline('exit')


def plot_bkgd_with_environmental(data, bkg):
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax = niutils.plotparams(ax)
    colors = ["#9749b9", "#669c55", "#b54958", "#757cad", "#bb8140"]
    colors_bkg = ['gray']*len(colors)

    df = pd.read_csv(data, skiprows=3, delimiter=" ", header=None)
    df.columns = ['energy', 'energy_err', 'counts', 'counts_err']

    df_bkg = pd.read_csv(bkg, skiprows=3, delimiter=" ", header=None)
    df_bkg.columns = ['energy', 'energy_err', 'counts', 'counts_err']

    ax.errorbar(df['energy'], df['counts'], 
                xerr=df['energy_err'], yerr=df['counts_err'],
                ls=' ', marker='.', color=colors[2], 
                label='60--80')

    ax.errorbar(df_bkg['energy'], df_bkg['counts'], 
                xerr=df_bkg['energy_err'], yerr=df_bkg['counts_err'],
                ls=' ', marker='.', color=colors_bkg[2], 
                label='BKG Model')

    ax.legend(edgecolor='black', fontsize=20)
    ax.set_ylabel(r'Normalized Counts (s$^{-1}$ keV$^{-1}$)', fontsize=20)
    ax.set_xlabel('Energy (keV)', fontsize=20)

    fig.savefig("first_environmental_attempt.pdf")

def plot_all_environmental(bkg_list, labels):
    fig, ax = plt.subplots(1, 1, figsize=(12,6))
    ax = niutils.plotparams(ax)

    for bkg, l in zip(bkg_list, labels):
        df_bkg = pd.read_csv(bkg, skiprows=3, delimiter=" ", header=None)
        df_bkg.columns = ['energy', 'energy_err', 'counts', 'counts_err']

        ax.errorbar(df_bkg['energy'], df_bkg['counts'], 
                    xerr=df_bkg['energy_err'], yerr=df_bkg['counts_err'],
                    ls=' ', marker='.',
                    label=l)

    ax.legend(edgecolor='black', fontsize=20)
    ax.set_ylabel(r'Normalized Counts (s$^{-1}$ keV$^{-1}$)', fontsize=20)
    ax.set_xlabel('Energy (keV)', fontsize=20)

    plt.show()

def plot_bkgd_spectra(data_list, ranges, environ=None, bkg3c50=None):
    fig, (ax, ax1) = plt.subplots(2, 1, figsize=(12, 12))
    plt.subplots_adjust(hspace=.1)
    ax = niutils.plotparams(ax)
    ax1 = niutils.plotparams(ax1)

    colors = niutils.get_colors()['br_colors']

    if environ is not None:
        data_list.append(environ)
        colors.append(niutils.get_colors()['environ_color'])
        ranges.append("Environmental Model")

    if bkg3c50 is not None:
        data_list.append(bkg3c50)
        colors.append(niutils.get_colors()['3c50_color'])
        ranges.append("3C50 Model")


    for i in range(len(data_list)):
        df = pd.read_csv(data_list[i], skiprows=3, delimiter=" ", header=None)
        df.columns = ['energy', 'energy_err', 'counts', 'counts_err']

        ax.errorbar(df['energy'], df['counts'], 
                    xerr=df['energy_err'], yerr=df['counts_err'],
                    ls=' ', marker='.', color=colors[i], 
                    label=ranges[i])

        ax1.errorbar(df['energy'], df['counts'], 
                    xerr=df['energy_err'], yerr=df['counts_err'],
                    ls=' ', marker='.', color=colors[i], 
                    label=ranges[i])


    ax.legend(edgecolor='black', fontsize=20)

    for a in [ax,ax1]:
        a.set_ylabel(r'Normalized Counts (s$^{-1}$ keV$^{-1}$)', fontsize=20)
        a.set_xscale('log')
        a.set_yscale('log')
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    ax1.xaxis.set_minor_formatter(mticker.FormatStrFormatter('%.01f'))

    ax1.set_xlabel('Energy (keV)', fontsize=20)

    ax1.set_xlim(left=.28, right=.82)
    ax1.set_ylim(top=.92, bottom=.3)

    con1 = ConnectionPatch(xyA=(ax1.get_xlim()[0], ax1.get_ylim()[1]),
                           xyB=(ax1.get_xlim()[0], ax1.get_ylim()[0]),
                           coordsA='data', coordsB='data',
                           axesA=ax1, axesB=ax, 
                           color='black', ls='-')

    con2 = ConnectionPatch(xyA=(ax1.get_xlim()[1], ax1.get_ylim()[1]),
                           xyB=(ax1.get_xlim()[1], ax1.get_ylim()[0]),
                           coordsA='data', coordsB='data',
                           axesA=ax1, axesB=ax,
                           color='black', ls='-')


    rect = Rectangle((ax1.get_xlim()[0], ax1.get_ylim()[0]), 
                     (ax1.get_xlim()[1]-ax1.get_xlim()[0]),
                     (ax1.get_ylim()[1]-ax1.get_ylim()[0]), 
                     facecolor='none', edgecolor='black', 
                     alpha=.5)
    ax.add_patch(rect)
    ax1.add_patch(con1)
    ax1.add_patch(con2)

    fig.savefig("BKGD_spectra.pdf")


def make_spectra(): 
    evtlist = glob.glob('*.evt') 
    evtlist = ['br_earth_--40.evt', 'br_earth_40--60.evt', 
               'br_earth_60--80.evt',
               'br_earth_80--180.evt', 'br_earth_180--.evt']
    nchan = [-999, 1000, 1000, 1000, 1000] 
    for i in range(len(evtlist)): 
        genspectra.gen_bkgd_spectra(evtlist[i], nchan[i], 
                                    save_pha=evtlist[i].replace('evt', 'pha')) 

if __name__ == '__main__':
    pha_list = ['br_earth_--40.pha', 'br_earth_40--60.pha', 
                'br_earth_60--80.pha',
                'br_earth_80--180.pha', 'br_earth_180--.pha']
    ranges = ['0--40', '40--60', '60--80', '80--180', '180--200']
    """
    make_spectra()

    pha_list = ['br_earth_--40.pha', 'br_earth_40--60.pha', 
                'br_earth_60--80.pha',
                'br_earth_80--180.pha', 'br_earth_180--.pha']

    ranges = ['0--40', '40--60', '60--80', '80--180', '180--200']
    for pha in pha_list:
        xspec_routine(pha)
    """
    plot_bkgd_spectra([pha.replace('.pha', '_data.txt') for pha in pha_list ],
                      ranges, 
                      environ='environ_all/data_environ_all.txt', 
                      bkg3c50='3c50/data_3c50.txt')

    #plot_bkgd_with_environmental('br_earth_60--80_data.txt', 'br_earth_60--80_bkg_data.txt')

