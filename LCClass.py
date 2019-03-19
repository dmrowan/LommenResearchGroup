#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import argparse
import collections
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
import os
import pexpect
import sys
import subprocess
from astropy.table import Table
import spectraplots
#Dom Rowan 2019

class LightCurve:
    def __init__(self, evtfile):
        assert(type(evtfile) == str)
        assert(os.path.isfile(evtfile))

        self.tab = Table.read(evtfile, hdu=1)
        self.pi = self.tab['PI']
        self.ph = self.tab['PULSE_PHASE']
        self.counts = None

    def mask(self, lower_pi=0, upper_pi=10, lower_ph=0, upper_ph=1):
        en_mask = (self.pi > lower_pi) & (self.pi < upper_pi)
        ph_mask = (self.ph > lower_ph) & (self.ph < upper_ph)
        full_mask = en_mask & ph_mask
        self.pi = self.pi[full_mask]
        self.ph = self.ph[full_mask]
    def generate(self, n_phase=2, bs=.01):
        self.bs = bs 
        self.phasebins = np.array([round(b, 4) for b in np.arange(0, 1, bs)])
        self.counts = np.zeros(len(self.phasebins))
        for phase in self.ph:
            idx_phasebins = np.where(self.phasebins <= phase)[0].max()
            self.counts[idx_phasebins] +=1

        self.counts_extended = np.array([])
        for i in range(n_phase):
            self.counts_extended = np.append(
                    self.counts_extended, self.counts)
        self.phasebins_extended = np.array(
                [round(b, 2) for b in np.arange(0, n_phase, bs)])

    def plot(self, n_phase=2, bs=.01, output=None, 
             show=False, extension='pdf', l1=None, l2=None):
        if output is not None:
            assert(type(output) == str)
            assert(type(extension) == str)

        if self.counts is None:
            self.generate(n_phase, bs=bs)

        fig, ax = plt.subplots(1, 1, figsize=(8,4))
        plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
        ax.plot(self.phasebins_extended, self.counts_extended, marker='.',
                ls='-', color='xkcd:violet')
        ax = spectraplots.plotparams(ax)
    
        if l1 is not None and l2 is not None:
            cutofftup = self.peak_cutoff(l1, l2)
            ax.axhline(cutofftup.median, ls='--', color='gray')
            ax.axhline(cutofftup.threesigma, ls=':', color='gray')
            ax.axvline(cutofftup.min_phase, ls='--', color='black')
            ax.axvline(cutofftup.max_phase, ls='--', color='black')
        ax.set_xlabel('Phase', fontsize=25)
        ax.set_ylabel('Counts', fontsize=25)
        

        if output is not None:
            fig.savefig(f"{output}.{extension}")
        elif show:
            plt.show()
        else:
            return ax
    def peak_cutoff(self, l1, l2):
        if self.counts is None:
            self.generate()
        off_pc = [ self.counts[i] for i in 
                            range(len(self.phasebins)) 
                            if l1 <= self.phasebins[i] <= l2 ]
        on_peak_phases = self.phasebins[np.where(
            self.counts >= (np.median(off_pc) + 3*np.std(off_pc)))[0]]

        wp = []
        for i in range(len(on_peak_phases)):
            if on_peak_phases[i] > .9:
                wp.append(on_peak_phases[i])
        for i in range(len(on_peak_phases)):
            if i != len(on_peak_phases) - 1:
                if on_peak_phases[i+1] - on_peak_phases[i] > .02:
                    max_phase = on_peak_phases[i]
                    break
        if len(wp) == 0:
            min_phase = on_peak_phases[0]
        else:
            if max(wp) >= .98:
                min_phase = wp[0]
                wp.reverse()
                for i in range(len(wp)):
                    if i != len(wp) -1:
                        if abs(wp[i+1] - wp[i]) > .02:
                            min_phase = wp[i]
        CutoffTup = collections.namedtuple('CutoffTup', 
                ['min_phase', 'max_phase', 'median', 'threesigma'])
        tup = CutoffTup(min_phase, max_phase, np.median(off_pc), 
                        np.median(off_pc)+3*np.std(off_pc))

        return tup

if __name__ == '__main__':
    lc = LightCurve("PSR_B1821-24_combined.evt")
    lc.plot(output="1821LC")
