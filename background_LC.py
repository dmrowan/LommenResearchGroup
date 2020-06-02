#!/usr/bin/env python

import argparse
from astropy import log
import datetime
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
import pandas as pd
import random
from tqdm import tqdm

import niutils

rc('text', usetex=True)

desc="""
Quick tests of LCs for BKGD RXTE
"""

def convert_time(time):
    timezero = datetime.datetime(year=2014, month=1, 
                                 day=1, hour=0, minute=0, second=0)
    new_time = timezero+datetime.timedelta(seconds=time)

    return new_time

def lc_from_evt(fname, savefig=None):

    tab = Table.read(fname, hdu=1)
    mintime = tab['TIME'].min()
    maxtime = tab['TIME'].max()
    nbins = len(tab)//1000
    bins = np.linspace(mintime, maxtime, nbins)

    counts = np.zeros(len(bins))

    for i in tqdm(range(len(tab))):
        idx = np.where(tab['TIME'][i] >= bins)[0][-1]
        counts[idx] += 1 

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax = niutils.plotparams(ax)

    timezero = datetime.datetime(year=2014, month=1, 
                                 day=1, hour=0, minute=0, second=0)
    bins = [ timezero + datetime.timedelta(seconds=b) for b in bins ]

    ax.scatter(bins, counts, color='xkcd:azure')
    ax.set_xlabel("Time (s)", fontsize=20)
    ax.set_ylabel("Number of Events", fontsize=20)
    ax.xaxis.set_minor_locator(MonthLocator())


    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)

def events_per_second_plot(flist, savefig=None, ax=None):
    log.info("Generating EXP per s plot")
    d = {}
    for f in tqdm(flist):
        tab = Table.read(f, hdu=1)
        date = tab.meta['DATE-OBS']
        exposure = tab.meta['EXPOSURE']
        table_length = len(tab)

        if exposure == 0:
            continue
        if date in d.keys():
            d[date] += table_length/exposure
        else:
            d[date] = table_length/exposure

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        plt.subplots_adjust(top=.98, right=.98)
        ax = niutils.plotparams(ax)
        created_fig=True
    else:
        created_fig=False

    dates = list(d.keys())
    format_specifier = '%Y-%m-%dT%H:%M:%S'
    dates_formatted = [ datetime.datetime.strptime(s, format_specifier)
                        for s in dates ]

    temp = []
    for i in range(len(dates_formatted)):
        dtemp = datetime.datetime(year=2017, month=11, day=1)
        if dates_formatted[i] < dtemp:
            temp.append(d[dates[i]])

    print(len(temp))
    print(np.median(temp))
    print(np.mean(temp))

    #Need to fix this


    df = pd.DataFrame({'dates':dates_formatted,
                       'vals':list(d.values())})
    df.sort_values('dates', inplace=True)
    print(df)


    color = niutils.get_colors()['br_colors'][-1]
    ax.plot_date(df['dates'], df['vals'], color=color, ls='-')

    ax.xaxis.set_minor_locator(MonthLocator())

    ax.set_xlabel("Observation Start Date", fontsize=20)
    ax.set_ylabel(r'Count Rate (s$^{-1}$)', fontsize=20)

    if (savefig is None) and (created_fig):
        plt.show()
    elif savefig is not None:
        fig.savefig(savefig)
    else:
        return ax


def cumulative_exposure_plot(flist, savefig=None, ax=None):
    log.info("Generating cumulative exposure plot")

    d = {}
    format_specifier = '%Y-%m-%dT%H:%M:%S'

    for f in tqdm(flist):
        tab = Table.read(f, hdu=1)
        date = tab.meta['DATE-OBS']
        exposure = tab.meta['EXPOSURE']
        if exposure == 0:
            continue

        date_formatted = datetime.datetime.strptime(date, format_specifier)
        date_formatted = date_formatted.replace(hour=0, minute=0, second=0)
        if date_formatted in d.keys():
            d[date_formatted] += exposure
        else:
            d[date_formatted] = exposure

    sorted_dates = list(d.keys())
    sorted_dates.sort()

    mindate = min(list(d.keys()))
    maxdate = max(list(d.keys()))

    iterdate = mindate - datetime.timedelta(1)
    date_list = []
    exposure_cumulative = []
    while iterdate <= maxdate + datetime.timedelta(1):

        date_list.append(iterdate)

        if len(exposure_cumulative) == 0:
            exposure_next = 0
        elif iterdate in d.keys():
            exposure_next = exposure_cumulative[-1] + d[iterdate]
        else:
            exposure_next = exposure_cumulative[-1]

        exposure_cumulative.append(exposure_next)

        iterdate = iterdate + datetime.timedelta(1)


    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        plt.subplots_adjust(top=.98, right=.98, bottom=.16)
        ax = niutils.plotparams(ax)
        created_fig=True
    else:
        created_fig = False

    exposure_cumulative = [ e/(10e10) for e in exposure_cumulative ]
    ax.plot(date_list, exposure_cumulative, 
            color=niutils.get_colors()['3c50_color'], lw=3)
    
    ax.xaxis.set_minor_locator(MonthLocator())

    ax.set_xlabel("Observation Start Date", fontsize=20)
    ax.set_ylabel(r'Cumulative Expsoure ($\times10^{10}$ s)', fontsize=20)

    if (savefig is None) and (created_fig):
        plt.show()
    elif savefig is not None:
        fig.savefig(savefig)
    else:
        return ax


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--evt", help="evt for LC", default=None, type=str)
    parser.add_argument("--cexp", help='Generate plot of cumulative exposure',
                        default=False, action='store_true')
    parser.add_argument("--eps", help="Generate plot of events per second",
                        default=False, action='store_true')
    parser.add_argument("--path", help="evt list for plot", 
                        default=None, type=str)
    parser.add_argument("--abs_path", help="Append absoulte path for cexp",
                        default=None, type=str)
    parser.add_argument("--savefig", help="output name",
                        default=None, type=str, nargs='+')
    args = parser.parse_args()

    if args.evt is not None:
        lc_from_evt(args.evt, savefig=args.savefig)

    if (args.cexp) or (args.eps):
        with open(args.path, 'r') as f:
            paths = f.readlines()
        paths = [ p.strip('\n') for p in paths ]
        if args.abs_path is not None:
            paths = [ f'{args.abs_path}{p}' for p in paths ]
        
    if args.cexp:
        cumulative_exposure_plot(paths, savefig=args.savefig[0])

    if args.eps:
        events_per_second_plot(paths, savefig=args.savefig[-1])

