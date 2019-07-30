#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import astropy.fitting
import collections
from collections.abc import Iterable
import math
from math import log10, floor
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib import rc
import numpy as np
import pickle
import os
from astropy.table import Table
from fuzzywuzzy import process
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy import exp
import niutils

#Dom Rowan and Lauren Lugo 2019

desc="""
Class for pulsar profiles and light curves
"""

rc('text', usetex=True)

class LightCurve:
    def __init__(self, evtfile):
        assert(type(evtfile) == str)
        if not os.path.isfile(evtfile):
            raise FileNotFoundError(
                    "Cmon Lauren use the right evt path \nhttps://tinyurl.com/yylzpd92")

        self.tab = Table.read(evtfile, hdu=1)
        self.pi = self.tab['PI']
        self.ph = self.tab['PULSE_PHASE']
        self.piratio = self.tab['PI_RATIO']
        self.counts = None # Initialize counts to none
        self.name = process.extract(evtfile, 
                                    ['PSR B1821-24', 'PSR B1937+21', 
                                     'PSR J0218+4232'],
                                    limit=1)[0][0]
 
        self.filters = {
                        'energy':[],
                        'phase':[],
                        'trumpet':[]
                       }

    # Apply energy mask
    def mask(self, lower_pi=0, upper_pi=10, lower_ph=0, upper_ph=1):
        en_mask = (self.pi > lower_pi) & (self.pi < upper_pi)
        ph_mask = (self.ph > lower_ph) & (self.ph < upper_ph)
        full_mask = en_mask & ph_mask
        self.pi = self.pi[full_mask]
        self.ph = self.ph[full_mask]

        self.filters['energy'].append([lower_pi, upper_pi])
        self.filters['phase'].append([lower_ph, upper_ph])

    def test_trumpet_cut(self, fconst, fastsig=1200, fastquart=0, 
                         n=1, ax=None, plot=False):
        
        #Define trumpet function
        def trumpet_cut(pi, c, s, q, n=1):
            return c + (s/10)/pi**n + q*pi**3

        #Create mask
        mask = [ (self.piratio[i] < trumpet_cut(self.pi[i], 
                                                fconst, 
                                                fastsig, 
                                                fastquart, 
                                                n=n))
                 for i in range(len(self.piratio)) ]

        mask_flip = [ not l for l in mask ]

        #Store the number of photons cut
        self.n_cut = len(mask) - sum(mask)
        
        if self.n_cut != 0:
            log.info("Appling trumpet cut")

        #Apply the mask
        self.pi = self.pi[mask]
        self.ph = self.ph[mask]
        self.piratio = self.piratio[mask]

        if plot:
            if ax is None:
                fig, ax = plt.subplots(1, 1, figsize=(16, 8))
                created_fig=True
            else:
                created_fig=False
            ax.scatter(self.pi, self.piratio, color='xkcd:blue', marker='.',
                        alpha=.5)

            ax = niutils.plotparams(ax)
            if created_fig:
                plt.show()
            else:
                return ax

        else:
            return ax

    
    # Apply Trumpet Cut using this mask and changing the threshold	
    def TrumpetMask(self,fastconst = 1.1, fileName = 'newFile.fits'):
        oldList = []
        #fastconst = 1.1 is the overall ratio threshold (normal ratio is 1.0 with a tolerance of 0.1=10%)
        # I recommend changing the threshold number blc.ploty small amounts to see the difference in info
        if fastconst == 1.18:
            print('hello1')
            t = np.arange(0,1201,1)
            newLine = []
            for pi in t:
                newLine.append(fastconst +((1200/10)/(pi**1.5))+ (5e-12)*pi**3)
        else:
            t = np.arange(0,1201,1)
            newLine = []
            for pi in t:
                newLine.append(fastconst +((1200/10)/(pi)))
		
	# Creating temporary storage for data that will be lost once the mask is        # applied
        # Creating temporary storage for data that will be 
        # lost once the mask is applied
        oldData =[]
        oldEnergy =[]
        #Creating a temp array to make a mask
        t_mask = []
        erase = []

        # loop goes through the rows one by one
        for py in range(len(self.piratio)):

            #Creating a mask by deciding on boolean
            #newLine holds a value for a given energy
            #The orginal boolean is piratio less than newLine
            t_mask.append(self.piratio[py] < newLine[self.pi[py]]) 
            #Saving all "false" information so it can be plotted before
            #The original boolean is piratio greater than newLine
            if (self.piratio[py] > newLine[self.pi[py]]):
                #print (self.pi[py])
                erase.append(py)
                oldData.append(self.piratio[py])
                oldEnergy.append(self.pi[py])
         
        oldList.append(oldEnergy)
        oldList.append(oldData)

        #Applying the mask
        self.piratio = self.piratio[t_mask]
        self.pi = self.pi[t_mask]
        self.ph = self.ph[t_mask]
        
        #removing the rows that ly outside the trumpet cut        
        self.tab.remove_rows(erase)
        #print(len(self.tab['PI']))
                  # num = num-1

        self.newFile = self.tab.write(fileName, format = 'fits', overwrite = True) 
        if fastconst == 1.18:
            oldList.append(self.pi)
            oldList.append(self.piratio)
            return oldList
        else:
            return oldList
        #return  self.newFile
        #run loop to create differnet cuts then return a list of evt files
        return self.newFile
        # run loop to create differnet cuts then return a list of evt files
        #call them in newCreateSpecs.py
        plt.scatter(oldEnergy,oldData, s=1)
        #colormag = np.vstack([self.pi,self.piratio])
        #z = gaussian_kde(colormag)(colormag)
        plt.ylim(bottom = 0)
        plt.scatter(self.pi,self.piratio, s= 1, label = 'newTrumpetCut')
        plt.show()  
	
    #unfinished
    def multiTrumpetCut(self):
	#The first cuts that were made were 1.075,1.05,1.025,1.015
        fastConst = [1.075, 1.05, 1.025, 1.015, 1.18]
        oldStuff = []
        #This for loop takes the length of the list of fast constants and makes a cut for each of the trumpet cut changes
        #the original trupet cut has a fast constant of 1.1
        for i in range(len(fastConst)):
            #We append a list of lists that contain the data of the trumpet cut once it has been preformed.
            
           oldStuff.append(self.TrumpetMask(fastConst[i],f'newfile_{i}.fits'))
           subprocess.call(['mv', f'newfile_{i}.fits', f'newfilecray1821_{i}.evt'])
           print(f'newfilecray1821_{i}.evt')
           subprocess.call(['mv',f'newfilecray1821_{i}.evt', 'evtFiles'])
        fig, ax= plt.subplots(1,1, figsize = (8,6))
        for j in range(len(oldStuff)):
            print(j)
            if j == 4:
                ax.scatter(oldStuff[j][0],oldStuff[j][1], s =1, label = fastConst[j])
                ax.scatter(oldStuff[j][2],oldStuff[j][3], s= 1, label = 'remainder', color = 'navy')
            else:
                ax.scatter(oldStuff[j][0],oldStuff[j][1],s =1,label = fastConst[j])
        ax.set_ylim(0.5,2.5)
        ax.set_title('1821 Trumpet Cuts')
        ax.set_xlabel('PI')
        ax.set_ylabel('PI_RATIO')
        plt.legend()
        fig.savefig('1821AllTrumpet.pdf')


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
        ax = niutils.plotparams(ax)
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
            if cutofftup.min_phase_p2 is not None:
                for i in range(n_phase):
                    ax.axvspan(cutofftup.min_phase_p2+i, 
                               cutofftup.max_phase_p2+i,
                               **default_span)

            if cutofftup.min_phase_p1 > cutofftup.max_phase_p1:
                for i in range(n_phase):
                    ax.axvspan(i, cutofftup.max_phase_p1+i, **default_span)
                    ax.axvspan(cutofftup.min_phase_p1+i, i+1, **default_span)
            else:
                for i in range(n_phase):
                    ax.axvspan(cutofftup.min_phase_p1+i, 
                               cutofftup.max_phase_p1+i, 
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
    def peak_cutoff(self, l1, l2, nsigma=3, 
                    interpulse_lower=0.3, verbose=False):

        if self.counts is None:
            self.generate()
        
        assert(type(l1) == type(l2))

        #If we have one bkgd range
        if not isinstance(l1, Iterable):
            assert(l1 < l2)
            if 0 <= l1 <= l2 <= 1:
                off_pc = [ self.counts[i] for i in
                           range(len(self.phasebins))
                           if l1 <= self.phasebins[i] <= l2 ]
            else:
                # ex l1=.85 l2 =1.15
                off_pc = [ self.counts[i] for i in 
                           range(len(self.phasebins))
                           if ( ( self.phasebins[i] >= l1) or
                                ( self.phasebins[i] <= (l2-1.0) ) ) ]
        #Using two bkgd ranges
        else: 
            assert(l1[1] > l1[0])
            assert(l2[1] > l2[0])

            if l1[1] <= 1:
                off_pc = [ self.counts[i] for i in 
                           range(len(self.phasebins))
                           if l1[0] <= self.phasebins[i] <= l1[1] ]
            else:
                off_pc = [ self.counts[i] for i in 
                           range(len(self.phasebins))
                           if ( ( self.phasebins[i] >= l1[0] ) or
                                ( self.phasebins[i] <= (l1[1]-1) ) ) ]

            if l2[1] <= 1:
                off_pc.extend([self.counts[i] for i in
                               range(len(self.phasebins))
                               if l2[0] <= self.phasebins[i] <= l2[1]])
            else:
                off_pc.extend([self.counts[i] for i in
                               range(len(self.phasebins))
                               if ( (self.phasebins[i] >= l2[0]) or
                                    (self.phasebins[i] <= (l2[1]-1) ) ) ])


        #collect phases where counts are greater than nsigma away
        on_peak_phases = self.phasebins[np.where(
            self.counts >= (np.median(off_pc) + nsigma*np.std(off_pc)))[0]]

        #Determine if pulse overlaps 0 phase
        if (0.0 in on_peak_phases) and (1.0 - self.bs in on_peak_phases):
            self.wraps_zero = True
        else:
            self.wraps_zero = False

        groups = [ [on_peak_phases[0]] ]
        for i in range(1, len(on_peak_phases)):
            if on_peak_phases[i] <= groups[-1][-1] + 2*self.bs:
                groups[-1].append(on_peak_phases[i])
            else:
                groups.append([on_peak_phases[i]])

        if self.wraps_zero:
            groups[0] = groups[-1]+groups[0]
        while len(groups) > 2:
            lengths = [ len(g) for g in groups ]
            for i in range(len(groups)):
                if len(groups[i]) == min(lengths):
                    del groups[i]
                    break

        phase_max = self.phasebins[ np.where(
            self.counts == max(self.counts))[0][0] ]
        if groups[1][0] <= phase_max <= groups[1][-1]:
            groups = groups[::-1]

        """
        if self.wraps_zero:
            print(groups[0])
            groups[0] = [ groups[0][i] + 1 if groups[0][i] < groups[0][0] 
                          else groups[0][i]
                          for i in range(len(groups[0])) ]
            print(groups[0])
        """

        #Define return tuple
        CutoffTup = collections.namedtuple('CutoffTup', 
                ['min_phase_p1', 'max_phase_p1', 
                 'min_phase_p2', 'max_phase_p2', 
                 'median', 'nsigma', 'n'])
        tup = CutoffTup(groups[0][0], groups[0][-1],
                        groups[1][0], groups[1][-1],
                        np.median(off_pc), 
                        np.median(off_pc)+nsigma*np.std(off_pc),
                        nsigma)
        
        if verbose:
            print(f"Min phase primary: {tup.min_phase_p1}")
            print(f"Max phase primary: {tup.max_phase_p1}")
            print(f"Min phase interpulse: {tup.min_phase_p2}")
            print(f"Max phase interpulse: {tup.max_phase_p2}")

        return tup

    def pulse_centers(self, l1=None, l2=None):
        if l1 or l2 is None:
            if self.name == 'PSR B1821-24':
                l1 = (.85, 1.15)
                l2 = (.4, .6)
            elif self.name == 'PSR B1937+21':
                l1 = (.90, 1.20)
                l2 = (.45, .75)
            elif self.name == 'PSR J0218+4232':
                l1 = .25
                l2 = .35

        cutoff_tup = self.peak_cutoff(l1, l2)
            
        #Pulse 1
        idx_p1 = np.where(np.array(self.counts_extended) == 
                          max(self.counts_extended))[0][0]
        p1_center = self.phasebins_extended[idx_p1]

        #Pulse 2
        pulse2_counts = [ self.counts_extended[i] for i in 
                          range(len(self.counts_extended))
                          if (cutoff_tup.min_phase_p2 
                              <= self.phasebins_extended[i] 
                              <= cutoff_tup.max_phase_p2) ]

        idx_p2 = np.where(np.array(self.counts) == max(pulse2_counts))[0][0]
        p2_center = self.phasebins[idx_p2]

        return p1_center, p2_center

    def fit_two(self, model, ax=None, annotate=True,
                      output_fit=False, label=False):
        if self.counts is None:
            self.generate()
        

        phasebins_fitting = self.phasebins
        counts_fitting = self.counts

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 4))
            created_fig = True

        else:
            created_fig = False

        #ax.set_xlim(left=phase_min-.025, right=phase_max+.025)
        ax = niutils.plotparams(ax)

        if self.name == 'PSR B1821-24':
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.5
            p0_sigma_0 = .01
            p0_sigma_1 = .01
            p0_x0_0 = self.pulse_centers()[0]
            p0_x0_1 = self.pulse_centers()[1]

            bounds = ([0, self.pulse_centers()[0]-.1, 0, 0, 
                       self.pulse_centers()[1]-.1, 0, 0],
                      [np.inf, self.pulse_centers()[0]+.1, 1, np.inf,
                       self.pulse_centers()[1]+.1, 1, max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, 
                 p0_a_1, p0_x0_1, p0_sigma_1, p0_b]

            phase_min = 0
            phase_max = phase_min + 1

        #Initial values and bounds for 1937
        elif self.name == 'PSR B1937+21':
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.1
            p0_sigma_0 = .01
            p0_sigma_1 = .01
            p0_x0_0 = self.pulse_centers()[0]
            p0_x0_1 = self.pulse_centers()[1]

            bounds =([0, self.pulse_centers()[0]-.1, 0, 0,
                      self.pulse_centers()[1]-.1, 0, 0],
                     [np.inf, self.pulse_centers()[0]+.1, 1, np.inf, 
                      self.pulse_centers()[1]+.1, 1, max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, 
                 p0_a_1, p0_x0_1, p0_sigma_1, p0_b]

            phase_min = 0
            phase_max = phase_min + 1

        #Initial values and bounds for J0218
        else:
            p0_b = min(counts_fitting)
            p0_a_0 = max(counts_fitting) - p0_b
            p0_a_1 = p0_a_0 *0.5
            p0_sigma_0 = .3
            p0_sigma_1 = .3
            p0_x0_0 = self.pulse_centers()[0]+1.0
            p0_x0_1 = self.pulse_centers()[1]

            bounds =([0,      p0_x0_0-.1, 0, 0,      p0_x0_1-.1, 0, 0],
                     [np.inf, p0_x0_0+.1, 1, np.inf, p0_x0_1+.1, 1, 
                      max(counts_fitting)])

            p0= [p0_a_0, p0_x0_0, p0_sigma_0, 
                 p0_a_1, p0_x0_1, p0_sigma_1, p0_b]

            phase_min = 0.3 
            phase_max = phase_min + 1.0
            phasebins_fitting = np.array([ p for p in self.phasebins_extended 
                                           if phase_min <= p < phase_max ])
            counts_fitting = np.array([ 
                self.counts_extended[i] for i in range(len(self.counts_extended)) 
                if phase_min <= self.phasebins_extended[i] < phase_max ])


        valid = []
        for i in range(len(p0)):
            valid.append( bounds[0][i] <= p0[i] <= bounds[1][i] )
        if not all(valid):
            print("Fit failed due to invalid bounds")
            for i in range(len(p0)):
                print(bounds[0][i], p0[i], bounds[1][i])
            print(self.filters)
            raise ValueError

        #Perform scipy curve fit 
        model = process.extract(model, ['gaussian', 'lorentzian'],
                                limit=1)[0][0]
                        
        if model == 'gaussian':
            popt, pcov = curve_fit(niutils.two_gaus, phasebins_fitting, 
                                   counts_fitting, 
                                   p0=p0,
                                   bounds=bounds)
            fit_counts = niutils.two_gaus(phasebins_fitting, *popt)

        elif model == 'lorentzian':

            p0[2] = p0[2]*2.355
            p0[5] = p0[5]*2.355
            popt, pcov = curve_fit(niutils.two_lorentzians, 
                                   phasebins_fitting, 
                                   counts_fitting, 
                                   p0=p0, 
                                   bounds=bounds)
            
            fit_counts = niutils.two_lorentzians(phasebins_fitting, *popt)



        fit_counts_extended = np.append(fit_counts, fit_counts)
        fit_counts_extended = np.append(fit_counts_extended, fit_counts)
        phasebins_fitting_extended = np.array(
                [round(b,4) for b in np.arange(0, 3, self.bs) ])
        if phase_min != 0:
            phasebins_fitting_extended = np.array(
                [ pb - (1-phase_min) for pb in phasebins_fitting_extended ])

        ax.plot(phasebins_fitting_extended, fit_counts_extended, 
                color='xkcd:azure', lw=6, zorder=1, alpha=.4)

        ax.plot(self.phasebins_extended, self.counts_extended, 
                marker='.', color='xkcd:violet', zorder=2)

        plt.setp(ax.get_xticklabels()[0], visible=False)
        plt.setp(ax.get_xticklabels()[-1], visible=False)
        ax.set_xlim(0,2)

        """
        n_phase = 3
        fit_counts_extended = np.array([])
        for i in range(n_phase):
            fit_counts_extended = np.append(fit_counts_extended, fit_counts)
        phasebins_fitting_extended = np.array([round(b,4) - 1.0 
            for b in np.arange(phase_min, phase_max+n_phase-1, self.bs) ])
        """

        if annotate:
            ax.text(.95, .95, f"{round(popt[1], 4)}, {round(popt[4], 4)}", 
                    fontsize=20, transform=ax.transAxes, 
                    ha='right', va='top')


        PoptTup = collections.namedtuple('PoptTup',
                ['primary_amplitude', 'primary_position', 
                 'primary_sigma', 'secondary_amplitude', 
                 'secondary_position', 'secondary_sigma',
                 'vertical_shift'])

        popt_tup = PoptTup(*popt)
                   
        #Add component labels
        if label:
            if self.name == 'PSR B1821-24':
                p1_coords = (popt_tup.primary_position-8*
                             popt_tup.primary_sigma,
                             popt_tup.primary_amplitude*.8
                             + popt_tup.vertical_shift)
                p2_coords = (popt_tup.secondary_position-4*
                              popt_tup.secondary_sigma,
                              popt_tup.secondary_amplitude*1.25 
                              + popt_tup.vertical_shift)
            elif self.name == 'PSR B1937+21':
                p1_coords = ((popt_tup.primary_position)+10*
                             popt_tup.primary_sigma,
                             popt_tup.primary_amplitude * 
                             0.75 + popt_tup.vertical_shift)
                p2_coords = ((popt_tup.secondary_position)+8*
                              popt_tup.secondary_sigma,
                              popt_tup.secondary_amplitude*1.35 
                              + popt_tup.vertical_shift)

            else:

                p1_coords = ((popt_tup.primary_position-1)+2*
                             popt_tup.primary_sigma,
                             popt_tup.primary_amplitude * 
                             0.8 + popt_tup.vertical_shift)
                p2_coords = ((popt_tup.secondary_position)+1.5*
                              popt_tup.secondary_sigma,
                              popt_tup.secondary_amplitude*1.1
                              + popt_tup.vertical_shift)

            ax.text(p1_coords[0], p1_coords[1], "P1", 
                    fontsize=23, ha='center', va='center')

            ax.text(p2_coords[0], p2_coords[1], "P2", 
                    fontsize=23, ha='center', va='center')

        if output_fit:
            return phasebins_fitting_extended, fit_counts_extended, popt_tup
        elif created_fig:
            plt.show()
            return popt_tup
        else:
            return ax, popt_tup

    def height_ratio(self):
        popt = self.fit_two('lorentz')
        return popt.primary_amplitude / popt.secondary_amplitude

    def astropy_fit(self):
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
        counts_fitting = np.array([
            self.counts_extended[i] for i in range(len(self.counts_extended))
            if phase_min <= self.phasebins_extended[i] < phase_max ])

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 4))
            created_fig = True

        else:
            created_fig = False

        ax.set_xlim(left=phase_min-.025, right=phase_max+.025)
        ax = niutils.plotparams(ax)
        p0_b = min(counts_fitting)
        p0_a_0 = max(counts_fitting) - p0_b
        p0_a_1 = p0_a_0 *0.5
        p0_sigma_0 = .3
        p0_sigma_1 = .3
        if self.peak_center()[0] > .5:
            p0_x0_0 = self.peak_center()[0]
        else:
            p0_x0_0 = self.peak_center()[0]+1.0
        p0_x0_1 = self.interpulse_center()[0]+1.0

        m_init = niutils.two_gauss()
        fit = astropy.modeling.fitting.LevMarLSQFitter()
        model = fit(m_init, phasebins_fitting, counts_fitting)

        ax.plot(phasebins_fitting, counts_fitting)
        ax.plot(phasebins_fitting, model(phasebins_fitting),
                color='xkcd:azure', lw=6, zorder=1, alpha=.4)
        ax.plot(self.phasebins_extended, self.counts_extended,
                marker='.', color='xkcd:violet', zorder=2)

        plt.show()

