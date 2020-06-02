#!/usr/bin/env python
import os
import argparse
from astropy.io import fits
from astropy import log
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import subprocess
import sys
import pint
import pandas as pd
import numpy as np
import pickle

import niutils

#Dom Rowan 2019

rc('text', usetex=True)

def mkf_hist(table, key, hdu=1, bins=50, pickle_file=None, save=None):

	if pickle_file is None:
		tab = Table.read(table, hdu=hdu)
	else:
		tab = pickle.load(open(pickle_file, 'rb'))
	fig, ax = plt.subplots(1, 1, figsize=(12, 6))
	ax = niutils.plotparams(ax)

	assert(key in tab.colnames)

	ax.hist(tab[key], bins=bins, edgecolor='black', 
            color=niutils.get_colors()['br_colors'][-1])
	ax.set_ylabel("N Rows", fontsize=20)
	ax.set_xlabel(key, fontsize=20)

	if save is not None:
		fig.savefig(save)
	else:
		plt.show()

def br_earth_hist(table, bins=50, save=None):

    tab = Table.read(table, hdu=1)
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax = niutils.plotparams(ax)
    
    plt.subplots_adjust(top=.98, right=.98)

    br_earth = [ b for b in tab['BR_EARTH'] if b < 180 ]
    ax.hist(br_earth, bins=bins, edgecolor='black', 
            color=niutils.get_colors()['br_colors'][-1])
    ax.set_ylabel("N Rows", fontsize=20)
    ax.set_xlabel('Bright Earth Angle', fontsize=20)
    plt.setp(ax.get_yticklabels()[0], visible=False)

    if save is not None:
        fig.savefig(save)
    else:
        plt.show()

def br_hist(table, bins=50, save=None):
    
    tab = Table.read(table, hdu=1)

    fig, ax = plt.subplots(2, 1, figsize=(12, 12))
    for a in ax: 
        a = niutils.plotparams(a)
        a.set_xlabel("BR\_EARTH", fontsize=25)
        a.set_ylabel("NRows", fontsize=25)

    ax[0].hist(tab[np.where(tab['SUNSHINE']==1)[0]]['BR_EARTH'],
               bins=bins, color='#FF5903', alpha=.7, edgecolor='black')
    ax[1].hist(tab[np.where(tab['SUNSHINE']==0)[0]]['BR_EARTH'],
               bins=bins, color='#6902CB', edgecolor='black')

    ax[0].text(.95, .95, r'SUNSHINE$=1$', ha='right', va='top',
               transform=ax[0].transAxes, fontsize=20)
    ax[1].text(.95, .95, r'SUNSHINE$=0$', ha='right', va='top',
               transform=ax[1].transAxes, fontsize=20)

    #Save or show
    if save is not None:
        fig.savefig(save)
    else:
        plt.show()
    

    
def reg_earth_hists(table, bins=50, save=None):
    tab = Table.read(table, hdu=1)

    fig, ax = plt.subplots(2, 1, figsize=(12, 12))
    for a in ax: 
        a = niutils.plotparams(a)
        a.set_xlabel("BR\_EARTH", fontsize=25)
        a.set_ylabel("NRows", fontsize=25)

    ax[0].hist(tab[np.where(tab['SUNSHINE']==1)[0]]['ELV'],
               bins=bins, color='#FF5903', alpha=.7, edgecolor='black')
    ax[1].hist(tab[np.where(tab['SUNSHINE']==0)[0]]['ELV'],
               bins=bins, color='#6902CB', edgecolor='black')

    ax[0].text(.95, .95, r'SUNSHINE$=1$', ha='right', va='top',
               transform=ax[0].transAxes, fontsize=20)
    ax[1].text(.95, .95, r'SUNSHINE$=0$', ha='right', va='top',
               transform=ax[1].transAxes, fontsize=20)

    #Save or show
    if save is not None:
        fig.savefig(save)
    else:
        plt.show()

if __name__ == '__main__':
    br_earth_hist('/students/pipeline/heasoft6.26/bkgd_merged.mkf',
                  save='BR_EARTH_hist.pdf')

