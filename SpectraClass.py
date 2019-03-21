#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import Table
import spectraplots
#Dom Rowan 2019

#Class for NICER pulsar spectra generated from evt or pha files
class Spectra:
    def __init__(self, fname, filters={}):
        assert(os.path.isfile(fname))
        assert(type(fname) == str)
        assert(type(filters) == dict)

        #Init counts to none
        self.counts = None

        #Different init for each filetype
        if fname.endswith('evt'):
            self.ftype = 'evt'
            self.init_evt(fname)
        else:
            self.ftype = 'pha'
            self.init_pha(fname)

        #Filters is a directory to store filter info
        self.filters = filters
        self.exposure = self.get_meta()['EXPOSURE']
        
        self.fname = fname

    #Initialization for event file
    def init_evt(self, evtfile):
        assert(os.path.isfile(evtfile))
        assert(evtfile.endswith('.evt'))

        self.evtfile = evtfile
        self.tab = Table.read(evtfile, hdu=1)

    #Initialization for spectra pha file
    def init_pha(self, phafile):
        assert(os.path.isfile(phafile))
        assert(phafile.endswith('.pha') or phafile.endswith('.pi'))

        self.phafile = phafile
        self.tab = Table.read(self.phafile, hdu=1)
        
        #Already contains energy and counts information
        self.keV = np.true_divide(np.array(self.tab['CHANNEL']), 100)
        self.counts = np.array(self.tab['COUNTS'])

    #Apply energy mask
    def filter_energy(self, energy_lower, energy_upper):
        assert(all([type(val) in [int, float] for 
               val in [energy_lower, energy_upper]]))
        #Different procedures for each filetype
        if self.ftype == 'evt':
            self.filter_energy_evt(energy_lower, energy_upper)
        else:
            self.filter_energy_pha(energy_lower, energy_upper)
    
    #Energy filter for event files
    def filter_energy_evt(self, pi_lower, pi_upper):
        self.tab = self.tab[np.where( (self.tab['PI'] >= pi_lower) &
                                      (self.tab['PI'] <= pi_upper) )[0]]
        #Add applied filter to dictionary
        self.filters['Energy'] = (pi_lower, pi_upper)
    
    #Energy filter for pha files
    def filter_energy_pha(self, channel_lower, channel_upper):
        self.tab = self.tab[np.where( 
            (self.tab['CHANNEL'] >= channel_lower) & 
            (self.tab['CHANNEL'] <= channel_upper) )[0]]
        self.keV = np.true_divide(np.array(self.tab['CHANNEL']), 100)
        self.counts = np.array(self.tab['COUNTS'])

        #Add applied filter to dictionary
        self.filters['Energy'] = (channel_lower, channel_upper)
    
    
    #Apply phase filter (only available for event files)
    def filter_phase(self, phase_lower, phase_upper):
        assert(all([type(val) in [int, float] for 
               val in [phase_lower, phase_upper]]))

        if self.ftype == 'pha':
            print("No phase filtering available for pha file input")
        else:

            assert(self.ftype == 'evt')
            #Different indexing depending on phase wrapping
            if (0 <= phase_lower <= phase_upper <= 1):
                self.tab = self.tab[np.where( 
                        (self.tab['PULSE_PHASE'] >= phase_lower) &
                        (self.tab['PULSE_PHASE'] <= phase_upper) )[0]]
                #Add applied filter to dictionary
                self.filters['Phase'] = (phase_lower, phase_upper)
            elif (0 <= phase_lower <= 1 <= phase_upper):
                self.tab = self.tab[np.where(
                        (self.tab['PULSE_PHASE'] <= (phase_upper - 1)) |
                        (self.tab['PULSE_PHASE'] >= phase_lower))[0]]
                #Add applied filter to dictionary
                self.filters['Phase'] = (phase_lower, phase_upper)
            else:
                print("Invalid phase selection")
    
    #Bin photons by energy (i.e. produce spectra)
    def count_photons(self):
        #If file type is pha this step is unnecessary
        if self.ftype == 'pha':
            return self.keV, self.counts
        else:
            assert(self.ftype == 'evt')
            #Init energy axis
            energy = np.arange(self.tab['PI'].min(), 
                               self.tab['PI'].max()+1, 1)
            #Init counts
            counts = np.zeros(len(energy))
            #Fill counts iteratively
            for pi in self.tab['PI']:
                idx_energy = np.where(pi >= energy)[0].max()
                counts[idx_energy] += 1

            self.counts = counts
            #Use keV instead of channel PI
            self.keV = np.true_divide(energy, 100)
            return self.keV, self.counts
    
    #Return meta information
    def get_meta(self):
        d = {}
        #Clean up dictionary a bit
        for key in self.tab.meta:
            if key != '':
                d[key] = self.tab.meta[key]
        return d
    
    #Produce plot of spectra
    def plot(self, output=None, show=False, 
             extension='pdf', color='xkcd:violet'):

        if output is not None:
            assert(type(output) == str)
        assert(type(extension) == str)

        #Init matplotlib plot
        fig, ax = plt.subplots(1,1, figsize=(8,4))
        self.count_photons()
        ax.plot(self.keV, self.counts, marker='o', 
                ls='-', color=color)

        #Use default plot params
        ax = spectraplots.plotparams(ax)
        ax.set_xlabel("Energy (keV)", fontsize=25)
        ax.set_ylabel("Counts", fontsize=25)

        #Save show or return
        if output is not None:
            fig.saavefig(f"{output}.{extension}")
        elif show:
            plt.show()
        else:
            return ax
    
    def __repr__(self):
        if self.counts is None:
            self.count_photons()

        return(f"""
        Spectra for filename: {self.fname}
        Uses filters: {self.filters}
        """)
        

