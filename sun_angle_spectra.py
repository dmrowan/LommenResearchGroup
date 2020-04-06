#!/usr/bin/env python

from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd
import subprocess

import genspectra
import background_spectra
import niutils

#Dom Rowan 2020

rc('text', usetex=True)

def format_sa_range(r):
    split = r.split('--')
    assert(len(split)==2)
    if split[0] == '':
        return f'(SUN_ANGLE < {split[1]})'
    else:
        return f'(SUN_ANGLE > {split[0]})'

def make_sun_angle_gtis(mkf, br_range, sa_ranges):
    assert(os.path.isfile(mkf))

    formatted_ranges = [format_sa_range(r) for r in sa_ranges]

    gti_names = [f'gti_br_earth_{br_range}_SA_{s}.gti' for s in sa_ranges]

    for i in range(len(formatted_ranges)):
        cmd = ['nimaketime', 'infile={}'.format(mkf), 
                'outfile={}'.format(gti_names[i]),'nicersaafilt=YES',
                'saafilt=NO','trackfilt=YES', 'ang_dist=0.015', 'elv=20', 
                'br_earth=0', 'min_fpm=38', 'underonly_range=0-200', 
                'expr={0}'.format(formatted_ranges[i]), 'clobber=YES']
        subprocess.run(cmd)

def extract_sun_angle_events(gti, evt, mkf, sa_range):

    new_evt = gti.replace('gti', '')[1:]+'evt'
    new_mkf = new_evt.replace('evt', 'mkf3')

    cmd = ['niextract-events', f'filename={evt}', f'eventsout={new_evt}',
           f'timefile={gti}', 'gti=GTI', 'clobber=yes']

    subprocess.run(cmd)

    sa_expr = format_sa_range(sa_range)
    cmd = ['ftcopy', f'{mkf}[{sa_expr}]', new_mkf, 
           'clobber=yes', 'history=yes']

    subprocess.run(cmd)


def plot_sun_angle_fill(br_data_list, br_ranges, sa_ranges, 
                        savefig=None, text=None):

    fig, ax = plt.subplots(3, 2, figsize=(20, 15))
    plt.subplots_adjust(top=.98, right=.98, hspace=0, wspace=.1, bottom=.05)
    for i in range(len(br_data_list)):
        df_bottom = pd.read_csv(br_data_list[i][0], skiprows=3,
                                delimiter=" ", header=None)
        df_bottom.columns = ['energy', 'energy_err', 'counts', 'counts_err']
        df_top = pd.read_csv(br_data_list[i][1], skiprows=3,
                             delimiter=" ", header=None)
        df_top.columns = ['energy', 'energy_err', 'counts', 'counts_err']
        color = niutils.map_colors()[br_ranges[i]]
        dark_color = niutils.map_colors(dark=True)[br_ranges[i]]

        a = ax.reshape(-1)[i]
        a.fill_between(df_bottom['energy'],
                       df_bottom['counts'],
                       df_top['counts'], color='gray',zorder=1, alpha=.6)
        a.plot(df_bottom['energy'], df_bottom['counts'], color=color,
               label=niutils.sa_label(sa_ranges[i][0]), lw=3)
        a.plot(df_top['energy'], df_top['counts'], color=dark_color,
               label=niutils.sa_label(sa_ranges[i][1]), lw=3)
        a = niutils.plotparams(a)
        a.legend(edgecolor='black', fontsize=20, loc='upper right', 
                 bbox_to_anchor=(.98, .88), bbox_transform=a.transAxes)

        if i in [3,4]:
            a.set_xlabel('Energy (keV)', fontsize=20)
        else:
            a.set_xticklabels([])
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_ylim(bottom=.009, top=130)
        a.set_xlim(left=.15, right=11)
        a.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))

        a.text(.95, .95, niutils.br_label(br_ranges[i]), fontsize=20,
                ha='right', va='top', transform=a.transAxes)


    fig.text(.05, .55, r'Normalized Counts (s$^{-1}$ keV$^{-1}$)',
             fontsize=35, rotation='vertical', ha='center', va='center')
    ax.reshape(-1)[-1].axis('off')

    axt = ax.reshape(-1)[-1]
    axt.axis('off')
    axt.text(.05, .8, text, ha='left', va='top', fontsize=30)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)
        if savefig.endswith('pdf'):
            fig.savefig(savefig.replace('pdf', 'png'), dpi=300)


def wrapper():

    br_ranges = ['--40', '40--60', '60--80', '80--180', '180--']
    mkf_list = [f'br_earth_{r}.mkf3' for r in br_ranges]
    evt_list = [f'br_earth_{r}.evt' for r in br_ranges]


    for i in range(len(br_ranges)):
        make_sun_angle_gtis(mkf_list[i], br_ranges[i], ['--100', '100--'])

        extract_sun_angle_events(f'gti_br_earth_{br_ranges[i]}_SA_--100.gti',
                                 evt_list[i], mkf_list[i], 
                                 '--100')
        extract_sun_angle_events(f'gti_br_earth_{br_ranges[i]}_SA_100--.gti',
                                 evt_list[i], mkf_list[i], 
                                 '100--')
        genspectra.gen_bkgd_spectra(
                f'br_earth_{br_ranges[i]}_SA_--100.evt', 
                -999,
                save_pha=f'br_earth_{br_ranges[i]}_SA_--100.pha')
        genspectra.gen_bkgd_spectra(
                f'br_earth_{br_ranges[i]}_SA_100--.evt', 
                -999,
                save_pha=f'br_earth_{br_ranges[i]}_SA_100--.pha')


        background_spectra.xspec_routine(
                f'br_earth_{br_ranges[i]}_SA_--100.pha')
        background_spectra.xspec_routine(
                f'br_earth_{br_ranges[i]}_SA_100--.pha')

    
def wrapper2(gen=True):

    br_ranges = ['--40', '40--60', '60--80', '80--180', '180--']
    mkf_list = [f'br_earth_{r}.mkf3' for r in br_ranges]
    evt_list = [f'br_earth_{r}.evt' for r in br_ranges]

    sa_medians = []
    for mkf in mkf_list:
        tab = Table.read(mkf, hdu=1)
        sa_medians.append(int(round(np.median(tab['SUN_ANGLE']))))

    if gen:
        for i in range(len(br_ranges)):
            make_sun_angle_gtis(mkf_list[i], br_ranges[i], 
                                [f'--{sa_medians[i]}', f'{sa_medians[i]}--'])
            extract_sun_angle_events(
                f'gti_br_earth_{br_ranges[i]}_SA_--{sa_medians[i]}.gti',
                evt_list[i], mkf_list[i], f'--{sa_medians[i]}')
            extract_sun_angle_events(
                f'gti_br_earth_{br_ranges[i]}_SA_{sa_medians[i]}--.gti',
                evt_list[i], mkf_list[i], f'{sa_medians[i]}--')

            genspectra.gen_bkgd_spectra(
                    f'br_earth_{br_ranges[i]}_SA_--{sa_medians[i]}.evt',
                    -999,
                    save_pha=f'br_earth_{br_ranges[i]}_SA_--{sa_medians[i]}.pha')
            genspectra.gen_bkgd_spectra(
                    f'br_earth_{br_ranges[i]}_SA_{sa_medians[i]}--.evt',
                    -999,
                    save_pha=f'br_earth_{br_ranges[i]}_SA_{sa_medians[i]}--.pha')

            background_spectra.xspec_routine(
                    f'br_earth_{br_ranges[i]}_SA_--{sa_medians[i]}.pha')
            background_spectra.xspec_routine(
                    f'br_earth_{br_ranges[i]}_SA_{sa_medians[i]}--.pha')

    

    data_list = []
    sa_ranges = []
    for i in range(len(br_ranges)):
        data0 = f'br_earth_{br_ranges[i]}_SA_--{sa_medians[i]}_data.txt'
        data1 = f'br_earth_{br_ranges[i]}_SA_{sa_medians[i]}--_data.txt'
        data_list.append([data0, data1])
        sar = [f'--{sa_medians[i]}', f'{sa_medians[i]}--']
        sa_ranges.append(sar)

    text="For each bright earth angle range we extract\ntwo spectra: one with sun angle below the \nmedian, and one above. For each bright earth\nangle range we plot the two \nspectra with the area between shaded in"

    plot_sun_angle_fill(data_list, br_ranges, sa_ranges, text=text, 
                        savefig='sun_angle_spectra.pdf')
