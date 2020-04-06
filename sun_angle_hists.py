#!/usr/bin/env python

from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import colors as mpc
from matplotlib import rc
import numpy as np
from tqdm import tqdm

#Dom Rowan 2019

import niutils

rc('text', usetex=True)

def sun_angle_hist(mkf_list, ranges, savefig=None, multiple_panels=False, 
                   density=True):

    if multiple_panels:
        fig, ax = plt.subplots(3, 2, figsize=(20, 16))
        for a in ax.reshape(-1): a = niutils.plotparams(a)
        plt.subplots_adjust(top=.98, right=.98, hspace=0, wspace=.1)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax = niutils.plotparams(ax)

    for i in tqdm(range(len(mkf_list))):
        tab = Table.read(mkf_list[i], hdu=1)
        if multiple_panels:
            a = ax.reshape(-1)[i]
            a.text(.95, .95, niutils.br_label(ranges[i]), 
                   fontsize=20, ha='right', va='top',
                   transform=ax.reshape(-1)[i].transAxes, 
                   bbox=dict(facecolor='white', edgecolor='black',
                             boxstyle='round'))
            a.axvline(np.median(tab['SUN_ANGLE']), color='gray', ls=':',
                      lw=3)
            if i in [3,4]:
                a.set_xlabel('Sun Angle', fontsize=30)
            else:
                a.set_xticklabels([])
        else:
            a = ax

        a.hist(tab['SUN_ANGLE'], bins=50, facecolor='none', 
               edgecolor=niutils.map_colors()[ranges[i]], 
               label=niutils.br_label(ranges[i]), density=density, lw=4)
        a.set_xlim(50, 190)
        plt.setp(a.get_yticklabels()[0], visible=False)

    if density:
        ylabel='Normalized Number of Rows'
    else:
        ylabel='Number of Rows'

    if multiple_panels:
        fig.text(.05, .5, ylabel, fontsize=40,
                 rotation='vertical', ha='center', va='center')
        ax.reshape(-1)[-1].axis('off')
    else:
        ax.legend(edgecolor='black', fontsize=20)
        ax.set_xlabel("Sun Angle", fontsize=20)
        ax.set_ylabel(ylabel, fontsize=20)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)
        if savefig.endswith('.pdf'):
            fig.savefig(savefig.replace('.pdf', '.png'), dpi=300)



if __name__ == '__main__':
    ranges = ['--40', '40--60', '60--80', '80--180', '180--']
    mkf_list = [f'br_earth_{r}.mkf3' for r in ranges]
    sun_angle_hist(mkf_list, ranges, savefig='sun_angle_onepanel.pdf')
    sun_angle_hist(mkf_list, ranges, multiple_panels=True, density=False,
                   savefig='sun_angle_multipanel.pdf')
