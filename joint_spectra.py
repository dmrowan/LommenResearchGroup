#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib import rc
import numpy as np
import os
import pandas as pd
from fuzzywuzzy import process

from niutils import plotparams

#Dom Rowan 2019

desc="""
Plot a simulatneously modeled spectra for multiple telescopes
from XSPEC output

Example command line calls:

    $ joint_spectra.py data_all_three.txt -s "PSR B1821-24" -o OutputSpec.pdf

    $ joint_spectra.py datafile.txt # Creates a spec with no sourcelabel
                                         # And filename datafile_spec.pdf

Example function calls:

    >>> plot_joint_telescope("datafile.txt", "PSR B1821-24", "myspec.pdf")

No function return -- output spec is saved in the current directory. 

"""

#Setting style parameters for plots
rc('text', usetex=True)

#Class for parsing the joint data file from simultaneous fitting
class joint_data:

    def __init__(self, fname):
    
        if not os.path.isfile(fname):
            raise FileNotFoundError

        self.fname = fname
        self.source = fname

        self.parse_file()
   
    #Option to set sourcename unique from fname
    def set_source(self, source):
        self.source = source

    #Method to set the telescope names for the joint table
    def set_tel(self, i, name):
        try:
            type(self.tel_names) == list
        except:
            self.parse_file()

        if i >= len(self.tel_names):
            raise ValueError("Invalid telescope index")

        self.tel_names[i] = name


    #File parser method
    def parse_file(self):
        with open(self.fname) as f:
            lines = f.readlines()

        #Find all break indicies
        length_break = len(lines[4].split())
        break_string = ('NO '*length_break).rstrip()+"\n"
        breakidx = np.where(np.array(lines) == break_string)[0]

        #read in fist df seperately with skiprows=3
        df0 = pd.read_csv(self.fname, skiprows=3, delimiter=" ", header=None,
                          nrows=breakidx[0]-3)
        #reset column names
        try:
            df0.columns = ['energy', 'energy_err', 'counts', 
                           'counts_err', 'model']
        except:
            df0.columns = ['energy', 'energy_err', 'counts', 
                           'counts_err', 'model', '5', '6']

        self.df_list = [df0]

        #iterate through other indicies to create new data frames
        for i in range(len(breakidx)):
            if i != len(breakidx)-1:
                df = pd.read_csv(self.fname, skiprows=breakidx[i]+1, 
                                 delimiter=" ", header=None, 
                                 nrows=breakidx[i+1]-breakidx[i]-1)

            else:
                df = pd.read_csv(self.fname, skiprows=breakidx[i]+1, 
                                 delimiter=" ", header=None)

            try:
                df.columns = ['energy', 'energy_err', 
                               'counts', 'counts_err', 'model']
            except:
                df.columns = ['energy', 'energy_err', 
                               'counts', 'counts_err', 'model', '5', '6']

            self.df_list.append(df)

        #Check that we have the correct number of data frames & lengths
        self.check_parse()

        #Attributes for data and residuals
        self.data = [ self.df_list[i] 
                      for i in range(int(len(self.df_list)/2)) ]
        self.residuals = [ self.df_list[i] 
                           for i in range(
                           int(len(self.df_list)/2), len(self.df_list)) ]

        #Reset column names for residuals
        for df in self.residuals:
            try:
                df.columns = ['energy', 'energy_err',
                              'delchi', 'delchi_err', 'model']
            except:
                df.columns = ['energy', 'energy_err',
                              'delchi', 'delchi_err', 'model', '5', '6']



        #Setup empty telescope names list
        self.tel_names = [ "" for i in range(len(self.data)) ]

        return 0

    #Method to verify parse
    def check_parse(self):
        #file parse must be attempted before verification
        try:
            type(self.df_list) == list
        except:
            self.parse_file()

        #Use lengths to verify. Should be two dfs with each length
        df_lengths = [ len(df) for df in self.df_list ]
        if len(df_lengths) != 2*len(set(df_lengths)):
            raise ValueError("Error in file parse inconsistent lengths")

        #Order of lengths should go data dfs then residual dfs
        for i in range(int(len(self.df_list)/2)):
            if (len(self.df_list[i]) != 
                    len(self.df_list[i+int(len(self.df_list)/2)])):
                raise ValueError("Error in file parse inconsistent lengths")

        return 0

#Plotting routine for joint fit
def plot_joint_telescope(fname, source, output):

    #Init figure
    fig = plt.figure(figsize=(8.5, 5.5))
    plt.subplots_adjust(top=.98, right=.98, wspace=.1, 
                        left=.15, bottom=.18)
    #Outer gridspec of size 2

    #Each spec of outer contains ufspec and delchi
    inner = gridspec.GridSpec(2, 1, height_ratios=[3,1.2], hspace=0)

    #Make axes we can plot onto
    ax1 = plt.Subplot(fig, inner[1])
    ax0 = plt.Subplot(fig, inner[0], sharex=ax1)

    #Use joint data class to parse text file
    xd = joint_data(fname)
    xd.set_source(source)

    source = process.extract(source, [r'PSR B1937$+$21', r'PSR B1821$-$24',
                                      r'PSR J0218$+$4232'], limit=1)[0][0]
    if source == r'PSR B1937$+$21':
        labels = [ r'$NICER$', r'$NuSTAR$',r'$XMM-Newton$' ]

        colors = ["#d5483a",
                "#70c84c",
                "#853bce",
                #"#d4ae2f",
                #"#625cce",
                #"#c24ebe", 
                "xkcd:azure"]
    elif source == r'PSR B1821$-$24':
        labels = [ r'$NICER$', r'$NuSTAR$', r'$RXTE$']

        colors = ["#d5483a",
                "#70c84c",
                #"#853bce",
                #"#d4ae2f",
                "#625cce",
                #"#c24ebe", 
                "xkcd:azure"]

    elif source == r'PSR J0218$+$4232':
        labels = [ r'$NICER$', r'$NuSTAR$',r'$XMM-Newton$' ]

        colors = ["#d5483a",
                "#70c84c",
                "#853bce",
                #"#d4ae2f",
                #"#625cce",
                #"#c24ebe", 
                "xkcd:azure"]

    else:
        return -1

    zorders = [2, 1.5, 1]
    #Iterate through each data and residual pair
    for i in range(len(xd.data)):
        ax0.errorbar(xd.data[i]['energy'],
                     xd.data[i]['counts'],
                     xerr=xd.data[i]['energy_err'],
                     yerr=xd.data[i]['counts_err'],
                     ls=' ', marker='.', 
                     color=colors[i],
                     label=labels[i],
                     zorder=zorders[i])
        ax0.plot(xd.data[i]['energy'],
                 xd.data[i]['model'],
                 ls='-', lw=3, color=colors[i],
                 zorder=zorders[i],
                 label='_nolegend_')

        ax1.errorbar(xd.residuals[i]['energy'].astype(float),
                     xd.residuals[i]['delchi'].astype(float),
                     xerr=xd.residuals[i]['energy_err'].astype(float),
                     yerr=xd.residuals[i]['delchi_err'].astype(float),
                     ls=' ', marker='.', color=colors[i], alpha=0.8,
                     zorder=zorders[i])


    #Misc plot params
    ax0 = plotparams(ax0)
    ax1 = plotparams(ax1)
    #ax0.set_ylim(bottom=5e-7)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.text(.95, .95, source, transform=ax0.transAxes, 
             ha='right', va='top', fontsize=20)
    ax0.legend(loc=(.05, 0.05), fontsize=16, 
               edgecolor='black', framealpha=.9)
    ax1.axhline(0, ls=':', lw=1.5, color='gray')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    ax1.set_xlabel('Energy (keV)', fontsize=20)

    #Add subplots to figure
    fig.add_subplot(ax0)
    fig.add_subplot(ax1)

    #Dont want to show xtick labels for ufspec
    plt.setp(ax0.get_xticklabels(), visible=False)

    #ax_top0.set_ylim(ax_top0.get_ylim()[0] * .4)

    #Add axes labels
    ax0.set_ylabel(r'Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$',
                   ha='center', va='center', fontsize=20)
    ax1.set_ylabel(r'Residuals ($\chi$)', fontsize=15)
    ax0.yaxis.set_label_coords(
            -0.115, ax0.yaxis.get_label().get_position()[1])

    #Save figure
    fig.savefig(output, dpi=2000)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fname", help='Xspec output file', type=str)
    parser.add_argument("-s", help="PSR source name", dest='source',
                        type=str, required=False, default=None)
    parser.add_argument("-o", help="Output file", dest='output',
                        type=str, required=False, default=None)
    args = parser.parse_args()

    if args.output is None:
        base = os.path.splitext(os.path.basename(args.fname))[0]
        args.output = base + "_spec.pdf"


    args = parser.parse_args()
    
    plot_joint_telescope(args.fname, args.source, args.output)
