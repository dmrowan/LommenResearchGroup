#!/usr/bin/env python

import argparse
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import NullLocator

import niutils

rc('text', usetex=True)

desc="""
Parse and Plot ACE data
"""

#Make the column name formatted for y axis label
def format_ace_column(cname, table=None):
    #capitalize and remove underscore
    formatted= ' '.join([s.capitalize() for s in cname.split('_')])
    #If we are given a table, find the units as well
    if table is not None:
        unit_string = get_ace_units(table, cname)
        formatted = f'{formatted} {unit_string}'

    return formatted

#Get the units from the table 'header'
def get_ace_units(table, column):

    #Read in the file line by line
    with open(table, 'r') as f:
        lines = f.readlines()

    #Search for the part that lists column description
    a = np.array([l.strip('\n').split()[0] if l!='\n' else '' for l in lines])
    idx = np.where(a == 'fp_year')[0][0]

    #Iterate through and find the units at the end of the lines
    for i in range(idx, idx+13):
        if column in lines[i].split()[0]:
            if '(' and ')' in lines[i].split()[-1]:
                unitstring = lines[i].split()[-1]
                break

    #We want just the parenthesis
    while unitstring[-1] != ')':
        unitstring = unitstring[:-1]
    while unitstring[0] != '(':
        unitstring = unitstring[1:]

    #We need to change to latex math if we have an exponent
    if '^' in unitstring: 
        #Find where the carrot is
        idx_carrot = np.where(np.array(list(unitstring)) == '^')[0][0]
        #Get everything from before the exp
        front = unitstring[:idx_carrot]
        #We need to also have curly brackets if there are multiple numbers
        # ~~ just realized I need to make this work for powers>9 
        if unitstring[idx_carrot+1] == '-':
            back= '$^{'+unitstring[idx_carrot+1]+unitstring[idx_carrot+2]+'}$)'
        else:
            back= '$^{'+unitstring[idx_carrot+1]+'}$)'
        #Turn into r string with formatting
        unitstring = r'{0}'.format(front+back)

    return unitstring

def read_table(fname):

    with open(fname, 'r') as f:
        lines = f.readlines()

    a = np.array([l.strip('\n').split()[0] if l!='\n' else '' for l in lines])
    idx = np.where(a == 'year')[0][0]

    df = pd.read_csv(fname, header=None, skiprows=idx+2, 
                     delim_whitespace=True,
                     names=lines[idx].strip('\n').split())

    return df

def merge_tables(fnames):

    df_list = [read_table(f) for f in fnames]

    df_merged = df_list[0]

    if len(df_list) > 1:
        for i in range(1, len(df_list)):
            df_merged = df_merged.append(df_list[i])

    df_merged = df_merged.reset_index(drop=True)

    return df_merged

def add_datetime_col(df):

    dt_list = []

    for i in range(len(df)):
        
        dt = datetime.datetime(year=df['year'][i], month=1, day=1)
        dt = dt + datetime.timedelta(days=int(df['day'][i]), 
                                     hours=int(df['hr'][i]),
                                     minutes=int(df['min'][i]), 
                                     seconds=int(df['sec'][i]))

        dt_list.append(dt)

    df['datetime'] = dt_list

    return df


def plot_data(fnames, column, savefig=None):

    if niutils.check_iter(fnames):
        df = merge_tables(fnames)
        tab_for_format = fnames[0]
        ntab = len(fnames)
    else:
        df = read_table(fnames)
        tab_for_format = fname
        ntab = 1

    if column is None:
        raise ValueError("Column can not be None")

    df = add_datetime_col(df)


    fig, ax = plt.subplots(1, 1, figsize=(12*ntab, 6))
    ax = niutils.plotparams(ax)


    if not niutils.check_iter(column):
        column = [column]

    mask = [ np.where(df[column[i]] != -9999.9)[0] 
              for i in range(len(column)) ]

    colors = niutils.get_colors()['br_colors']
    colors = [colors[0], colors[2]]

    ax.plot(df['datetime'][mask[0]], df[column[0]][mask[0]], 
            color=colors[0], label=column[0])

    ax.set_xlabel("Date", fontsize=20)
    ax.set_ylabel(format_ace_column(column[0], table=tab_for_format), 
                  fontsize=20)

    if len(column) > 1:
        ax, ax1 = niutils.DoubleY(ax, colors=(colors[0], colors[1]))

        ax1.plot(df['datetime'][mask[1]], df[column[1]][mask[1]], 
                 color=colors[1], label=column[1])

        ax1.set_ylabel(format_ace_column(column[1], table=tab_for_format), 
                       fontsize=20)

    if max(df['day']) > 180:
        ax.xaxis.set_minor_locator(MonthLocator())
    else:
        ax.xaxis.set_major_locator(MonthLocator())
        ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
        myfmt = DateFormatter('%Y-%m')
        ax.xaxis.set_major_formatter(myfmt)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)

def plot_data_split(fnames, column, savefig=None):

    if niutils.check_iter(fnames):
        df_list = [read_table(f) for f in fnames]
    else:
        raise ValueError("Must input multiple filenames")

    if column is None:
        raise ValueError("Column can not be None")

    df_list = [add_datetime_col(df) for df in df_list]

    fig, ax = plt.subplots(len(df_list), 1, figsize=(12, 5.2*len(df_list)))
    plt.subplots_adjust(top=.98, bottom=.07)

    for a in ax: a = niutils.plotparams(a)

    if not niutils.check_iter(column):
        column = [column]

    colors = niutils.get_colors()['br_colors']
    colors = [colors[0], colors[2]]

    for i in range(len(df_list)):
        df = df_list[i]
        mask = [ np.where(df[column[i]] != -9999.9)[0] 
                  for i in range(len(column)) ]

        ax[i].plot(df['datetime'][mask[0]], df[column[0]][mask[0]], 
                      color=colors[0], label=column[0])

        ax[i].set_xlabel("Date", fontsize=20)
        ax[i].set_ylabel(format_ace_column(column[0], table=fnames[0]), 
                         fontsize=20)
        if len(column) > 1:
            ax[i], ax1 = niutils.DoubleY(ax[i], colors=(colors[0], colors[1]))

            ax1.plot(df['datetime'][mask[1]], df[column[1]][mask[1]], 
                     color=colors[1], label=column[1])

            ax1.set_ylabel(format_ace_column(column[1], 
                           table=fnames[0]), fontsize=20)
        if max(df['day']) > 180:
            ax[i].xaxis.set_minor_locator(MonthLocator())
        else:
            ax[i].xaxis.set_major_locator(MonthLocator())
            ax[i].xaxis.set_minor_locator(MonthLocator(bymonthday=15))
            myfmt = DateFormatter('%Y-%m')
            ax[i].xaxis.set_major_formatter(myfmt)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--table", help="table fname", type=str,
                        default=None, nargs='+')
    parser.add_argument("--column", help="column to plot", type=str,
                        default=None, nargs='+')
    parser.add_argument("--savefig", help="output plot name", type=str,
                        default=None)
    parser.add_argument("--split", default=False, help="Use multiple panels",
                        action='store_true')

    args = parser.parse_args()

    if args.split:
        plot_data_split(args.table, args.column, savefig=args.savefig)
    else:
        plot_data(args.table, args.column, savefig=args.savefig)

