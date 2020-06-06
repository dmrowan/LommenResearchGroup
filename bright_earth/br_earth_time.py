#!/usr/bin/env python

from astropy.table import Table
import niutils
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tqdm import tqdm


def dfsplit(tab, tbreak):
    breaks = []
    for i in range(len(tab['TIME'])):
        if i!=0:
            if tab['TIME'][i] - tab['TIME'][i-1] >= tbreak:
                breaks.append(i)

    data = np.split(tab, breaks)
    return data


def test(mkf_list):

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    found=False
    for mkf in tqdm(mkf_list):
        tab = Table.read(mkf, hdu=1)
        tab_list = dfsplit(tab, 2)

        for t in tab_list:
            if len(t) <=1500:
                continue

            if touches_200(t, 'BR_EARTH') and (not all(t['BR_EARTH'] == 100)):
                ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
                ax.plot(t['TIME']-t['TIME'].min(), t['ELV'])
                found=True
                break
        if found:break


    plt.show()

def paper_plot(mkf_list):

    """
    Go through all mkfs
    find out which split tables are longest in each of the four categories
    Plot the earth elevation and bright earth angles
    """


    fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    d = {'increasing':['', 0],
         'decreasing':['', 0],
         'touches200':['', 0],
         'localextrema':['', 0]}

    for mkf in tqdm(mkf_list):
        tab = Table.read(mkf, hdu=1)
        tab_list = dfsplit(tab, 2)

        for t in tab_list:
            if len(t) <=5:
                continue

            if  touches_200(t, 'BR_EARTH'):
                if len(t) > d['touches200'][1]:
                    d['touches200'][0] = t
                    d['touches200'][1] = len(t)
            elif always_increasing(t, 'BR_EARTH'):
                if len(t) > d['increasing'][1]:
                    d['increasing'][0] = t
                    d['increasing'][1] = len(t)
            elif always_decreasing(t, 'BR_EARTH'):
                if len(t) > d['decreasing'][1]:
                    d['decreasing'][0] = t
                    d['decreasing'][1] = len(t)
            elif local_extrema(t, 'BR_EARTH'):
                if len(t) > d['localextrema'][1]:
                    d['localextrema'][0] = t
                    d['localextrema'][1] = len(t)
            else:
                if len(groups) ==4:
                    print(f'Problem with groups for {mkf}')

    ax = niutils.plotparams(ax)

    colors = niutils.get_colors()['br_colors']
    darkcolors = niutils.get_colors()['br_colors_dark']

    for i in range(len(d.keys())):
        t = d[list(d.keys())[i]][0]
        ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'],
                color=colors[i], lw=2)
        ax.plot(t['TIME']-t['TIME'].min(), t['ELV'],
                color=darkcolors[i], lw=2)

    plt.show()

def main(mkf_list, groups=['increasing', 'decreasing', 
                           'extrema', 'touches']):


    fig, (ax, ax1) = plt.subplots(2, 1, figsize=(12, 12))
    plt.subplots_adjust(hspace=0)
    ax = niutils.plotparams(ax)
    ax1 = niutils.plotparams(ax1)

    for mkf in tqdm(mkf_list):
        tab = Table.read(mkf, hdu=1)
        tab_list = dfsplit(tab, 2)

        for t in tab_list:
            if len(t) <= 5:
                continue

            if  touches_200(t, 'BR_EARTH') and 'touches' in groups:
                ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            elif always_increasing(t, 'BR_EARTH') and 'increasing' in groups:
                ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            elif always_decreasing(t, 'BR_EARTH') and 'decreasing' in groups:
                ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            elif local_extrema(t, 'BR_EARTH') and 'extrema' in groups:
                ax.plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            else:
                if len(groups) ==4:
                    print(f'Problem with groups for {mkf}')
            #ax1.plot(t['TIME']-t['TIME'].min(), t['ELV'])

    ax1.set_xlabel("Time from start (s)", fontsize=20)
    ax.set_ylabel("Bright Earth Angle", fontsize=20)
    ax1.set_ylabel("Earth Elevation Angle", fontsize=20)

    plt.show()

def always_increasing(tab, key):
    for i in range(len(tab)):
        if i != len(tab)-1:
            if tab[key][i+1] <= tab[key][i]:
                return False
    return True

def always_decreasing(tab, key):
    for i in range(len(tab)):
        if i != len(tab)-1:
            if tab[key][i+1] >= tab[key][i]:
                return False
    return True

def local_extrema(tab, key):
    return (not always_increasing(tab, key)) and (not always_decreasing(tab, key))

def touches_200(tab, key):
    if any(tab[key] == 200):
        return True
    else:
        return False
    
def multi_window_plot(mkf_list, savefig=None):

    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    plt.subplots_adjust(left=.1, bottom=.08, right=.99, top=.97, wspace=.11, hspace=.13)
    for a in ax.reshape(-1): a = niutils.plotparams(a)

    for mkf in tqdm(mkf_list):
        tab = Table.read(mkf, hdu=1)
        tab_list = dfsplit(tab, 2)

        for t in tab_list:
            if len(t) <= 5:
                continue

            if  touches_200(t, 'BR_EARTH'):
                ax.reshape(-1)[0].plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            elif always_increasing(t, 'BR_EARTH'):
                ax.reshape(-1)[1].plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            elif always_decreasing(t, 'BR_EARTH'):
                ax.reshape(-1)[2].plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])
            else: 
                assert(local_extrema(t, 'BR_EARTH'))
                ax.reshape(-1)[3].plot(t['TIME']-t['TIME'].min(), t['BR_EARTH'])

    for a, l in zip(ax.reshape(-1), ['Touches 200', 'Increasing', 'Decreasing', 'Extrema']):
        a.text(.95, .95, l, transform=a.transAxes, ha='right', va='top', fontsize=20)

    fig.text(.05, .5, 'Bright Earth Angle', ha='center', 
             va='center', rotation='vertical', fontsize=25)
    fig.text(.55, .03, 'Time (s)', ha='center', 
             va='center', fontsize=25)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)
    
if __name__ == '__main__':

    with open('bkgd_all_mkf_list', 'r') as f:
        all_paths = f.readlines()

    all_paths = [a.strip('\n') for a in all_paths ]
    paper_plot(all_paths)

    with open('mini_mkf_list', 'r') as f:
        paths = f.readlines()

    paths = [p.strip('\n') for p in paths ]
