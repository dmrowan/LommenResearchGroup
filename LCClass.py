#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import collections
from math import log10, floor
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
import os
from astropy.table import Table
from fuzzywuzzy import process
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy import exp
import spectraplots

#Dom Rowan and Lauren Lugo 2019

def gaus(x, a, x0, sigma, b):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + b

def round_1sigfig(x):
    return round(x, -int(floor(log10(abs(x)))))

def two_gaus(x, 
             a_0, x0_0, sigma_0, 
             a_1, x0_1, sigma_1, b,
             ):
               
    return a_0*exp(-(x-x0_0)**2/(2*sigma_0**2)) + a_1*exp(-(x-x0_1)**2/(2*sigma_1**2)) + b

#LC Class for pulsar profile
class LightCurve:
    def __init__(self, evtfile):
        assert(type(evtfile) == str)
        assert(os.path.isfile(evtfile))

        self.tab = Table.read(evtfile, hdu=1)
        self.pi = self.tab['PI']
        self.ph = self.tab['PULSE_PHASE']
        self.piratio = self.tab['PI_RATIO']
        self.counts = None # Initialize counts to none
        self.name = process.extract(evtfile, 
                                    ['PSR B1821-24', 'PSR B1937+21', 
                                     'PSR J0218+4232'],
                                    limit=1)[0][0]
                                               

    # Apply energy mask
    def mask(self, lower_pi=0, upper_pi=10, lower_ph=0, upper_ph=1):
        en_mask = (self.pi > lower_pi) & (self.pi < upper_pi)
        ph_mask = (self.ph > lower_ph) & (self.ph < upper_ph)
        full_mask = en_mask & ph_mask
        self.pi = self.pi[full_mask]
        self.ph = self.ph[full_mask]
    
    # Apply Trumpet Cut using this mask and changing the threshold	
    def TrumpMask(self,fastconst = 1.1):
        #fastconst = 1.1 is the overall ratio threshold (normal ratio is 1.0 with a tolerance of 0.1=10%)
        # I recommend changing the threshold number by small amounts to see the difference in info
        
        t = np.arange(0,1201,1)
        newLine = []
        for pi in t:
            newLine.append(fastconst +((1200/10)/pi))
        # Creating temporary storage for data that will be lost once the mask is        # applied

        oldData =[]
        oldEnergy =[]
        #Creating a temp array to make a mask
        t_mask = []
        erase = []

        # loop goes through the rows one by one
        for py in range(len(self.piratio)):

            #Creating a mask by deciding on boolean
            t_mask.append(self.piratio[py] < newLine[self.pi[py]]) 
            #Saving all "false" information so it can be plotted before
            if (self.piratio[py] > newLine[self.pi[py]]):
                #print (self.pi[py])
                erase.append(py)
                oldData.append(self.piratio[py])
                oldEnergy.append(self.pi[py])

        #Applying the mask
        self.piratio = self.piratio[t_mask]
        self.pi = self.pi[t_mask]
        
	#removing the rows that ly outside the trumpet cut        
	self.tab.remove_rows(erase)
        print(len(self.tab['PI']))
                  # num = num-1

        self.newFile = self.tab.write(self.newFile, format = 'fits', overwrite = True)
        return self.newFile
        # run loop to create differnet cuts then return a list of evt files
        #call them in newCreateSpecs.py
        plt.scatter(oldEnergy,oldData, s=1)

        #colormag = np.vstack([self.pi,self.piratio])
        #z = gaussian_kde(colormag)(colormag)
        plt.ylim(bottom = 0)
        plt.scatter(self.pi,self.piratio, s= 1, label = 'newTrumpetCute')
        plt.show()  
	
    #unfinished
    def multiTrumpetCut(self):
        fastConst = [1.5,1.75]
        evtFiles = []
        for i in range(len(fastConst)):
           evtFiles.append(TrumpMask(fastConst))

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
    def plot(self, n_phase=2, bs=.01, 
             output=None, extension='pdf', 
             l1=None, l2=None, nsigma=3, ax=None,
             label=True):

        if output is not None:
            assert(type(output) == str)
            assert(type(extension) == str)

        if self.counts is None:
            self.generate(n_phase, bs=bs)

        if ax is None:
            #Initialize matplotlib figure
            fig, ax = plt.subplots(1, 1, figsize=(8,4))
            plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
            created_fig=True
        else:
            created_fig=False

        #Default plot paramaters
        ax = spectraplots.plotparams(ax)
        if label:
            ax.set_xlabel('Phase', fontsize=25)
            ax.set_ylabel('Counts', fontsize=25)
        ax.set_xlim(left=-.025, right=n_phase+.025)

        ax.plot(self.phasebins_extended, self.counts_extended, 
                marker='.', ls='-', color='xkcd:violet')

        #Use l1 and l2 to define an offpeak region & determine onpeak region
        if l1 is not None and l2 is not None:
            cutofftup = self.peak_cutoff(l1, l2, nsigma=nsigma)
            #ax.axhline(cutofftup.median, ls='--', color='gray')
            ax.axhline(cutofftup.nsigma, ls=':', color='darkblue', 
                       label=str(cutofftup.n)+r'$\sigma$')
            default_span = dict(alpha=.2, color='gray')
            if cutofftup.min_phase_ip is not None:
                for i in range(n_phase):
                    ax.axvspan(cutofftup.min_phase_ip+i, 
                               cutofftup.max_phase_ip+i,
                               **default_span)

            if cutofftup.min_phase > cutofftup.max_phase:
                for i in range(n_phase):
                    ax.axvspan(i, cutofftup.max_phase+i, **default_span)
                    ax.axvspan(cutofftup.min_phase+i, i+1, **default_span)
            else:
                for i in range(n_phase):
                    ax.axvspan(cutofftup.min_phase+i, cutofftup.max_phase+i, 
                               **default_span)
                    #ax.legend()
        
        #ax.legend(loc=(.85, .85), fontsize=20, edgecolor='black')
        if self.name is not None and label:
            ax.text(.95, .95, self.name, ha='right', va='top', 
                    transform=ax.transAxes, fontsize=20)
        #Save/display/return plot
        if output is not None:
            fig.savefig(f"{output}.{extension}", dpi=500)

        #If the figure is within class, show
        if created_fig:
            plt.show()
            return 0
        #If appending to input axis, return modifying axis
        else:
            return ax

    #Use offpeak region to determine pulse phase limits
    def peak_cutoff(self, l1, l2, nsigma=3, interpulse_lower=0.3):

        if self.counts is None:
            self.generate()

        #If using one bkgd range
        if all([ type(l) in [int, float] for l in [l1, l2] ]):
            if not all([ 0 <= l <= 1 for l in [l1, l2] ]):
                raise ValueError("Invalid phase val")

            if l1 >= l2:
                raise ValueError("Invalid phase range")

            #Collect all counts in selected off peak region
            off_pc = [ self.counts[i] for i in 
                       range(len(self.phasebins)) 
                       if l1 <= self.phasebins[i] <= l2 ]

        #If using two background ranges
        else:
            if not all([ type(l) in [list, tuple, np.ndarray] for l in [l1, l2] ]):
                raise TypeError("Phase limits must be int/float or list, tuple array")

            if not all( [ all([ 0 <= ll <= 1 for ll in l]) for l in [l1, l2]  ]):
                raise ValueError("Invalid phase val")

            if any( [l[0] >= l[1] for l in [l1, l2] ]):
                raise ValueError("invalid phase range")

            #Collect all counts in selected off peak region
            off_pc = [ self.counts[i] for i in range(len(self.phasebins)) 
                       if ( (l1[0] <= self.phasebins[i] <= l1[1]) 
                          or (l2[0] <= self.phasebins[i] <= l2[1]) ) ]

                        

        #collect phases where counts are greater than nsigma away
        on_peak_phases = self.phasebins[np.where(
            self.counts >= (np.median(off_pc) + nsigma*np.std(off_pc)))[0]]
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
            else:
                max_phase = on_peak_phases[i]
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
                if i != len(on_peak_phases) - 1:
                    if ((on_peak_phases[i] >= interpulse_lower) and 
                       (on_peak_phases[i+1] - on_peak_phases[i] < .02)):
                        min_phase_interpulse = on_peak_phases[i]
                        ip_found = True
        if not ip_found:
            min_phase_interpulse = None
            max_phase_interpulse = None


        #Define return tuple
        CutoffTup = collections.namedtuple('CutoffTup', 
                ['min_phase', 'max_phase', 
                 'min_phase_ip', 'max_phase_ip', 
                 'median', 'nsigma', 'n'])
        tup = CutoffTup(min_phase, max_phase, 
                        min_phase_interpulse,
                        max_phase_interpulse,
                        np.median(off_pc), 
                        np.median(off_pc)+nsigma*np.std(off_pc),
                        nsigma)

        return tup

    def peak_center(self):
        if self.counts is None:
            self.generate()

        idx = np.where(np.array(self.counts_extended) == max(self.counts_extended))[0]
        return self.phasebins_extended[idx]
       
    def interpulse_center(self):
        if self.counts is None:
            self.generate()
        
        interpulse_counts = [ self.counts[i] for i in range(len(self.counts))
                              if 0.2 <= self.phasebins[i] <= 0.8 ]

        idx = np.where(np.array(self.counts) == max(interpulse_counts))[0]
        return self.phasebins[idx]


    def fit_gauss(self, component, include_phases=None, ax=None):
        if self.counts is None:
            self.generate()
        
        #Parse user entered component type
        component = process.extract(component, ['primary', 'interpulse'], 
                                    limit=1)[0][0]
    
        if include_phases is None:
            if component == 'primary':
                phase_min = .75
                phase_max = 1.25
            elif component == 'interpulse':
                phase_min = .4 
                phase_max = .7
            else:
                print("Invalid component")
                return -1
        else:
            phase_min = include_phases[0]
            phase_max = include_phases[1]

        phasebins_fitting = np.array([ p for p in self.phasebins_extended 
                                       if phase_min <= p <= phase_max ])
        counts_fitting = np.array([ self.counts_extended[i] 
                                    for i in range(len(self.counts_extended)) 
                                    if phase_min <= self.phasebins_extended[i] <= phase_max ])

        

        p0_a = max(counts_fitting)
        if component == 'primary':
            p0_x0 = self.peak_center()[0] + 1.0
        elif component == 'interpulse':
            p0_x0 = self.interpulse_center()[0]
        else:
            print("Invalid component")
            return -1

        p0_sigma = 0.1
        p0_b = min(counts_fitting)
        popt, pcov = curve_fit(gaus, phasebins_fitting, 
                               counts_fitting, 
                               p0=[p0_a, p0_x0, p0_sigma, p0_b])

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 4))
            plt.subplots_adjust(bottom=.2, top=.98, right=.98, left=.15)
            created_fig = True
            ax.set_xlabel('Phase', fontsize=25)
            ax.set_ylabel('Counts', fontsize=25)
        else:
            created_fig=False

        ax = spectraplots.plotparams(ax)
        ax.set_xlim(left=phase_min-.025, right=phase_max+.025)

        ax.plot(phasebins_fitting, counts_fitting,
                marker='.', ls='-', color='xkcd:violet', zorder=2)
        ax.plot(phasebins_fitting, 
                gaus(phasebins_fitting, *popt),
                color='xkcd:azure', ls='-', zorder=1, alpha=.4, lw=6)

        chisq, p = chisquare(counts_fitting, gaus(phasebins_fitting, *popt))
        chisq = chisq / (len(counts_fitting)-len(popt))
        ratio = popt[0] / np.std(counts_fitting)
        if any(np.sqrt(np.diag(pcov)) == np.inf):
            ratio=-1.0

        ax.text(.95, .95, f"{self.name}", ha='right', va='top', 
                fontsize=20, transform=ax.transAxes)
        ax.text(.95, .85, f"Center: {round(popt[1], 3)}",
                ha='right', va='top', fontsize=15, transform=ax.transAxes)

        if created_fig:
            plt.show()
            return popt, ratio
        else:
            return ax, popt, ratio


        
    def fit_two_gauss(self, include_phases=None, ax=None, annotate=True):
        if self.counts is None:
            self.generate()
        
    
        if include_phases is None:
            phase_min = .75
            phase_max = 1.75
        else:
            phase_min = include_phases[0]
            phase_max = include_phases[1]

        phasebins_fitting = np.array([ p for p in self.phasebins_extended 
                                       if phase_min <= p < phase_max ])
        counts_fitting = np.array([ self.counts_extended[i] 
                                    for i in range(len(self.counts_extended)) 
                                    if phase_min <= self.phasebins_extended[i] < phase_max ])

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 4))
            created_fig = True

        else:
            created_fig = False

        ax.set_xlim(left=phase_min-.025, right=phase_max+.025)
        ax = spectraplots.plotparams(ax)

        if self.name == 'PSR B1821-24':
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.5
            p0_sigma_0 = .01
            p0_sigma_1 = .01
            p0_x0_0 = 1.0
            p0_x0_1 = 1.55

            bounds =([0,      0.9, 0, 0,      1.4, 0, 0],
                     [np.inf, 1.1, 1, np.inf, 1.6, 1, max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, p0_a_1, p0_x0_1, p0_sigma_1, p0_b]

        #Initial values and bounds for 1937
        elif self.name == 'PSR B1937+21':
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.1
            p0_sigma_0 = .01
            p0_sigma_1 = .01
            p0_x0_0 = self.peak_center()[0]+1.0
            p0_x0_1 = self.interpulse_center()[0]+1.0

            bounds =([0,      0.9, 0, 0,      1.4, 0, 0],
                     [np.inf, 1.1, 1, np.inf, 1.6, 1, max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, p0_a_1, p0_x0_1, p0_sigma_1, p0_b]

        #Initial values and bounds for J0218
        else:
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.5
            p0_sigma_0 = .3
            p0_sigma_1 = .3
            p0_x0_0 = self.peak_center()[0]+1.0
            p0_x0_1 = self.interpulse_center()[0]+1.0

            bounds =([0,      0.9, 0, 0,      1.4, 0, 0],
                     [np.inf, 1.1, 1, np.inf, 1.6, 1, max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, p0_a_1, p0_x0_1, p0_sigma_1, p0_b]


        try:
            popt, pcov = curve_fit(two_gaus, phasebins_fitting, 
                                   counts_fitting, 
                                   p0=p0,
                                   bounds=bounds)
        except:
            print("Fit failed due to invalid bounds")
            for i in range(len(p0)):
                print(bounds[0][i], p0[i], bounds[1][i])
            raise ValueError


        fit_counts = two_gaus(phasebins_fitting, *popt)

        n_phase = 3
        fit_counts_extended = np.array([])
        for i in range(n_phase):
            fit_counts_extended = np.append(fit_counts_extended, fit_counts)
        phasebins_fitting_extended = np.array([round(b,4) - 1.0 
            for b in np.arange(phase_min, phase_max+n_phase-1, self.bs) ])
        ax.set_xlim(0,2)

        plt.setp(ax.get_xticklabels()[0], visible=False)
        plt.setp(ax.get_xticklabels()[-1], visible=False)

        ax.plot(phasebins_fitting, counts_fitting)
        ax.plot(phasebins_fitting_extended,
                fit_counts_extended,
                 color='xkcd:azure', lw=6, zorder=1, alpha=.4)
        ax.plot(self.phasebins_extended, self.counts_extended, 
                marker='.', color='xkcd:violet', zorder=2)
        if annotate:
            ax.text(.95, .95, f"{round(popt[1], 4)}, {round(popt[4], 4)}", 
                    fontsize=20, transform=ax.transAxes, 
                    ha='right', va='top')

        PoptTup = collections.namedtuple('PoptTup',
                ['primary_amplitude', 'primary_position', 
                 'primary_sigma', 'interpulse_amplitude', 
                 'interpulse_position', 'interpulse_sigma',
                 'vertical_shift'])

        popt_tup = PoptTup(*popt)
                   


        if created_fig:
            plt.show()
            return popt_tup
        else:
            return ax, popt_tup

