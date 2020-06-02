#!/usr/bin/env python

import pexpect
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc
import matplotlib.ticker as mticker
import os
import niutils

rc('text', usetex=True)

def xspec_bsbs(rxte_bkg, bkg_3c50, bkg_environ):

    assert(os.path.isfile(rxte_bkg))
    assert(os.path.isfile(bkg_3c50))
    assert(os.path.isfile(bkg_environ))

    datafile_3c50 = rxte_bkg.replace('.pha', '_bsbs_3c50_data.txt')
    datafile_environ = rxte_bkg.replace('.pha', '_bsbs_environ_data.txt')

    for b, d in zip([bkg_3c50, bkg_environ], [datafile_3c50, datafile_environ]):
        xspec = pexpect.spawn('xspec')
        xspec.expect('XSPEC12>')

        xspec.sendline(f'data 1 {rxte_bkg}')
        xspec.expect('XSPEC12>')

        xspec.sendline(f'back 1 {b}')
        xspec.expect('XSPEC12>')
        xspec.sendline('cpd /xs')
        xspec.expect('XSPEC12>')
        xspec.sendline('ig **-.2')
        xspec.expect('XSPEC12>')
        xspec.sendline('setplot energy')
        xspec.expect('XSPEC12>')
        xspec.sendline('plot')
        xspec.expect('XSPEC12>')

        xspec.sendline('ipl')
        xspec.expect('PLT>')

        xspec.sendline(f"wdata {d}")
        if os.path.isfile(d):
            xspec.expect('exists, reuse it?')
            xspec.sendline("yes")

        xspec.expect('PLT>')
        xspec.sendline('exit')


def plot_bsbs(br_earth_range, datafile_3c50, datafile_environ):
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax = niutils.plotparams(ax)

    labels=['3c50', 'Environ']

    for d, l in zip([datafile_3c50, datafile_environ], labels):
        assert(os.path.isfile(d))

        df_bkg = pd.read_csv(d, skiprows=3, delimiter=" ", header=None)
        df_bkg.columns = ['energy', 'energy_err', 'counts', 'counts_err']


        ax.errorbar(df_bkg['energy'], df_bkg['counts'],
                    xerr=df_bkg['energy_err'], yerr=df_bkg['counts_err'],
                    ls=' ', marker='.',
                    label=l)

    ax.legend(edgecolor='black', fontsize=20)
    ax.set_ylabel(r'Normalized Counts (s$^{-1}$ keV$^{-1}$)', fontsize=20)
    ax.set_xlabel('Energy (keV)', fontsize=20)

    ax.text(.95, .95, br_earth_range, ha='right', va='top', fontsize=20, transform=ax.transAxes, 
            bbox=dict(facecolor='white', edgecolor='none', alpha=.6))

    ax.set_xscale('log')
    ax.set_xlim(.3, 4)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))

    plt.show()

if __name__ == '__main__':
    for br_earth_pha in ['br_earth_--40.pha', 'br_earth_40--60.pha', 'br_earth_60--80.pha',
                         'br_earth_80--180.pha', 'br_earth_180--.pha']:

        xspec_bsbs(br_earth_pha, '3c50/merged_3C50_bkg.pha', 'merge_environ/merged_environ.pha')

    #plot_bsbs(r, 'br_earth_60--80_bsbs_3c50_data.txt', 
    #          'br_earth_60--80_bsbs_environ_data.txt')
