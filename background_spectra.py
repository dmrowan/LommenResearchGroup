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

#Generate a text data file from spectra using xspec
def xspec_routine(pha):
    assert(os.path.isfile(pha))

    datafile = pha.replace('.pha', '_data.txt')

    #Start xspec child
    xspec = pexpect.spawn("xspec")
    xspec.expect('XSPEC12>')

    #Load data
    xspec.sendline(f'data 1 {pha}')
    xspec.expect('XSPEC12>')

    #Set plotting device
    xspec.sendline('cpd /xs')
    xspec.expect('XSPEC12>')

    #Choose energy range
    xspec.sendline('ig **-.2, 10.-**')
    xspec.expect('XSPEC12>')

    #Change axis
    xspec.sendline('setplot energy')
    xspec.expect('XSPEC12>')

    #Plot in xs
    xspec.sendline('plot data')
    xspec.expect('XSPEC12>')

    time.sleep(3)

    xspec.sendline('ipl')
    xspec.expect('PLT>')

    #Save txt file
    xspec.sendline(f"wdata {datafile}")
    if os.path.isfile(datafile):
        xspec.expect('exists, reuse it?')
        xspec.sendline("yes")
    
    xspec.expect('PLT>')
    xspec.sendline('exit')


def xspec_environ_spectra(bkg):
    assert(os.path.isfile(bkg))

    datafile = bkg.replace('.pha', '_data.txt')

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


def reformat_input(l):
    for i in range(len(l)):
        if type(l[i]) == str:
            l[i] = [l[i]]
    return l

def check_length(l):
    l = reformat_input(l)
    return len(l), [ len(elem) for elem in l ]

def plot_bkgd_spectra(br_data_list, br_ranges,
                      sa_ranges=None,
                      environ=None, bkg3c50=None,
                      two_panel=False, savefig=None,
                      zoom_top=0.3, zoom_bottom=0.9,
                      zoom_left=0.3, zoom_right=0.9):

    assert(len(br_data_list) == len(br_ranges))

    br_data_list = reformat_input(br_data_list)

    if any( [ l > 1 for l in check_length(br_data_list)[1] ] ):
        assert(sa_ranges is not None)

        sa_ranges = reformat_input(sa_ranges)
        assert(check_length(br_data_list) == check_length(sa_ranges))



    linestyles = ['-', '--', '-.', ':']

    data_flat = []
    br_ranges_flat = []
    sa_ranges_flat = []
    ls_flat = []

    for i in range(len(br_data_list)):
        for j in range(len(br_data_list[i])):
            data_flat.append(br_data_list[i][j])
            br_ranges_flat.append(br_ranges[i])
            if sa_ranges is not None:
                sa_ranges_flat.append(sa_ranges[i][j])
            ls_flat.append(linestyles[j])

    colors_flat = [ niutils.map_colors()[r] for r in br_ranges_flat ]
    if sa_ranges is None:
        labels_flat = [ niutils.br_label(r) for r in br_ranges_flat ]
    else:
        labels_flat = [ f'{niutils.br_label(br_ranges_flat[i])}; {niutils.sa_label(sa_ranges_flat[i])}' for i in range(len(br_ranges_flat)) ]

    if environ is not None:
        data_flat.append(environ)
        colors_flat.append(niutils.get_colors()['environ_color'])
        labels_flat.append("Environmental Model")
        ls_flat.append('-')

    if bkg3c50 is not None:
        data_flat.append(bkg3c50)
        colors_flat.append(niutils.get_colors()['3c50_color'])
        labels_flat.append('3C50 Model')
        ls_flat.append('-')

    if two_panel:
        fig, (ax, ax1) = plt.subplots(2, 1, figsize=(16, 16))
        plt.subplots_adjust(hspace=.1)
        ax_list = [ax, ax1]

    else:
        fig, ax = plt.subplots(1, 1, figsize=(15, 8))
        ax_list = [ax]

    for a in ax_list: a = niutils.plotparams(a)

    for i in range(len(data_flat)):
        df = pd.read_csv(data_flat[i], skiprows=3, delimiter=" ", header=None)
        df.columns = ['energy', 'energy_err', 'counts', 'counts_err']

        for a in ax_list:
            a.errorbar(df['energy'], df['counts'],
                       xerr=df['energy_err'], yerr=df['counts_err'],
                       ls=ls_flat[i], marker='.', color=colors_flat[i],
                       label=labels_flat[i])

    ax.legend(edgecolor='black', fontsize=20)

    for a in ax_list:
        a.set_ylabel(r'Normalized Counts (s$^{-1}$ keV$^{-1}$)', fontsize=20)
        a.set_xscale('log')
        a.set_yscale('log')
    ax_list[-1].set_xlabel('Energy (keV)', fontsize=20)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    if two_panel:
        ax1.xaxis.set_minor_formatter(mticker.FormatStrFormatter('%.01f'))
        ax1.set_xlim(zoom_bottom, zoom_top)
        ax1.set_ylim(zoom_left, zoom_right)

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

    if savefig is not None:
        fig.savefig(savefig)
        if savefig.endswith('pdf'):
            fig.savefig(savefig.replace('pdf', 'png'))

    else:
        plt.show()


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

    plot_bkgd_spectra(['br_earth_--40_data.txt',
                       'br_earth_40--60_data.txt', 
                       'br_earth_60--80_data.txt', 
                       'br_earth_80--180_data.txt', 
                       'br_earth_180--_data.txt'], 
                       ['--40', '40--60', '60--80', '80--180', '180--'],
                       environ='environ_all/data_environ_all.txt',
                       bkg3c50='3c50/data_3c50.txt',
                       savefig='BKGD_spectra.pdf')

