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

#Andrea adding her name for testing purposes.

#This is a comment that Dom made today!

#LC Class for pulsar profile
class LightCurve:
    def __init__(self, evtfile):
        assert(type(evtfile) == str)
        assert(os.path.isfile(evtfile))

        self.tab = Table.read(evtfile, hdu=1)
        self.pi = self.tab['PI']
        self.ph = self.tab['PULSE_PHASE']
        self.counts = None # Initialize counts to none
        self.name = None

    # Apply energy mask
    def mask(self, lower_pi=0, upper_pi=10, lower_ph=0, upper_ph=1):
        en_mask = (self.pi > lower_pi) & (self.pi < upper_pi)
        ph_mask = (self.ph > lower_ph) & (self.ph < upper_ph)
        full_mask = en_mask & ph_mask
        self.pi = self.pi[full_mask]
        self.ph = self.ph[full_mask]

    # Give a name to include in plots
    def set_name(self, name):
        self.name = name

    #Generate count/phasebin information
    def generate(self, n_phase=2, bs=.01):
        assert(n_phase > 0)
        self.bs = bs 
        self.n_phase = n_phase
        #Array of phase bins
        self.phasebins = np.array([round(b, 4) for b in np.arange(0, 1, bs)])
        #Initialize array of counts
        self.counts = np.zeros(len(self.phasebins))
        for phase in self.ph:
            idx_phasebins = np.where(self.phasebins <= phase)[0].max()
            self.counts[idx_phasebins] +=1

        #If n_phase is greater than 1 need to extend both axes
        if n_phase > 1:
            self.counts_extended = np.array([])
            for i in range(n_phase):
                self.counts_extended = np.append(
                        self.counts_extended, self.counts)
            self.phasebins_extended = np.array(
                    [round(b, 2) for b in np.arange(0, n_phase, bs)])
        else:
            self.counts_extended = None
            self.phasebins_extended = None

    # Produce plot of puslar profile
    def plot(self, n_phase=2, bs=.01, output=None, 
             show=False, extension='pdf', 
             l1=None, l2=None, nsigma=3):

        if output is not None:
            assert(type(output) == str)
            assert(type(extension) == str)

        if self.counts is None:
            self.generate(n_phase, bs=bs)

        #Initialize matplotlib figure
        fig, ax = plt.subplots(1, 1, figsize=(8,4))
        plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
        ax.plot(self.phasebins_extended, self.counts_extended, 
                marker='.', ls='-', color='xkcd:violet')

        #Default plot paramaters
        ax = spectraplots.plotparams(ax)
        ax.set_xlabel('Phase', fontsize=25)
        ax.set_ylabel('Counts', fontsize=25)
        ax.set_xlim(left=-.025, right=n_phase+.025)
    
        #Use l1 and l2 to define an offpeak region & determine onpeak region
        if l1 is not None and l2 is not None:
            cutofftup = self.peak_cutoff(l1, l2, nsigma=nsigma)
            #ax.axhline(cutofftup.median, ls='--', color='gray')
            ax.axhline(cutofftup.two_sigma, ls=':', color='darkblue', 
                       label=r'$2\sigma$')
            ax.axhline(cutofftup.three_sigma, ls=':', color='gray',
                       label=r'$3\sigma$')
            default_line = dict(ls='--', color='black')
            default_span = dict(alpha=.2, color='gray')
            for i in range(n_phase):
                ax.axvspan(cutofftup.min_phase_ip+i, cutofftup.max_phase_ip+i,
                           **default_span)

            if cutofftup.min_phase > cutofftup.max_phase:
                for i in range(n_phase):
                    ax.axvspan(i, cutofftup.max_phase+i, **default_span)
                    ax.axvspan(cutofftup.min_phase+i, i+1, **default_span)
            else:
                for i in range(n_phase):
                    ax.axvspan(cutofftup.min_phase+i, cutofftup.max_phase+i, 
                               **default_span)
        
        #ax.legend(loc=(.85, .85), fontsize=20, edgecolor='black')
        if self.name is not None:
            ax.text(.95, .95, self.name, ha='right', va='top', 
                    transform=ax.transAxes, fontsize=20,
                    bbox=dict(facecolor='white', edgecolor='black',
                              alpha=.5))
        #Save/display/return plot
        if output is not None:
            fig.savefig(f"{output}.{extension}", dpi=300)
        elif show:
            plt.show()
        else:
            return ax

    #Use offpeak region to determine pulse phase limits
    def peak_cutoff(self, l1, l2, nsigma=3, interpulse_lower=0.3):
        assert(all([ 0 <= l <= 1 for l in [l1, l2]]))
        assert(l1 < l2)
        if self.counts is None:
            self.generate()

        #Collect all counts in selected off peak region
        off_pc = [ self.counts[i] for i in 
                   range(len(self.phasebins)) 
                   if l1 <= self.phasebins[i] <= l2 ]

        #collect phases where counts are greater than nsigma away
        on_peak_phases = self.phasebins[np.where(
            self.counts >= (np.median(off_pc) + nsigma*np.std(off_pc)))[0]]
        print(on_peak_phases)
        #If phases wrap around 0 need to take extra care
        wp = []
        for i in range(len(on_peak_phases)):
            if on_peak_phases[i] > .9:
                wp.append(on_peak_phases[i])

        #Look for gap in selected phases to determine pulse region
        for i in range(len(on_peak_phases)):
            if i != len(on_peak_phases) - 1:
                if on_peak_phases[i+1] - on_peak_phases[i] > .02:
                    max_phase = on_peak_phases[i]
                    break
        #If we have no wrapped phases the min of pulse is first index
        if len(wp) == 0:
            min_phase = on_peak_phases[0]
        #Else repeat process for max in other direction
        else:
            if max(wp) >= .98:
                min_phase = wp[0]
                wp.reverse()
                for i in range(len(wp)):
                    if i != len(wp) -1:
                        if abs(wp[i+1] - wp[i]) > .02:
                            min_phase = wp[i]

        #Search for interpuse
        ip_found = False
        for i in range(len(on_peak_phases)):
            if ip_found:
                if i != len(on_peak_phases) - 1:
                    if on_peak_phases[i+1] - on_peak_phases[i] > .02:
                        max_phase_interpulse = on_peak_phases[i]
                        break
                else:
                    max_phase_interpulse = on_peak_phases[i]
            else:
                if ((on_peak_phases[i] >= interpulse_lower) and 
                   (on_peak_phases[i+1] - on_peak_phases[i] < .02)):
                    min_phase_interpulse = on_peak_phases[i]
                    ip_found = True

        #Define return tuple
        CutoffTup = collections.namedtuple('CutoffTup', 
                ['min_phase', 'max_phase', 
                 'min_phase_ip', 'max_phase_ip', 
                 'median', 'two_sigma', 'three_sigma'])
        tup = CutoffTup(min_phase, max_phase, 
                        min_phase_interpulse,
                        max_phase_interpulse,
                        np.median(off_pc), 
                        np.median(off_pc)+2*np.std(off_pc),
                        np.median(off_pc)+3*np.std(off_pc))

        #Doms new edit

        return tup

# My new edit
