#!/usr/bin/env python
import ast
import collections
import numpy as np
import os
import pandas as pd

import niutils

#Dom Rowan 2019
desc="""
Classes for parsing log and wdata files from Xspec

logfile: Handles collection of parameters, errors, and contours
"""


#Class for parsing the log output of Xspec
class logfile:
    
    #Initialize with path to log file
    def __init__(self, fname):

        if not os.path.isfile(fname):
            raise FileNotFoundError

        self.fname = fname
        self.lines = open(self.fname, 'r').readlines()

    #Returns photon index and column density values
    def get_params(self):

        #Find line where fit ended
        idx_params = [ i for i in range(len(self.lines))
                       if 'par  comp' in self.lines[i] ]
        if len(idx_params) == 0:
            return "No parameter information found in logfile"
        else:
            idx_params = idx_params[-1]

        #Find the photon index 
        phot_index = float([ self.lines[idx_params+1].split()[i+1] 
                             for i in range(len(self.lines[idx_params+1].split())) 
                             if 'PhoIndex' in self.lines[idx_params+1].split()[i] 
                            ][0])

        #Find the Hydrogen column density
        nH = float([ self.lines[idx_params+3].split()[i+1]
                     for i in range(len(self.lines[idx_params+3].split()))
                     if '10^22' in self.lines[idx_params+3].split()[i] ][0])

        return phot_index, nH

    #Returns errors on photon index and column density
    def get_errors(self):

        #Find the confidence interval on the parameters
        idx_conf = [ i for i in range(len(self.lines))
                     if 'error 1' in self.lines[i] ]
        if len(idx_conf) == 0:
            print("No error calculation found in logfile")
        else:
            idx_conf = idx_conf[-1]

        if "Cannot do error calc" in self.lines[idx_conf+1]:
            print("Error in XSPEC error calculation")
            phot_err= float('NaN')
            nH_err = float('NaN')
        else:
            phot_tuple = ast.literal_eval(self.lines[idx_conf+2].split()[-1])
            
            if "Warning" not in self.lines[idx_conf+3]:
                nH_tuple = ast.literal_eval(self.lines[idx_conf+3].split()[-1])
            else:
                nH_tuple = ast.literal_eval(self.lines[idx_conf+4].split()[-1])

            phot_err = np.mean( [ abs(val) for val in phot_tuple ] )
            nH_err = np.mean( [ abs(val) for val in nH_tuple ] )

        return phot_err, nH_err

    def phot_string(self):
        phot_index, _ = self.get_params()
        phot_err, _ = self.get_errors()
    
        val_rounded, err_rounded = niutils.round_sigfigs(phot_index, phot_err)

        return f"{val_rounded}" + r'$\pm$' + f"{err_rounded}"

    def nH_string(self):
        _, nH = self.get_params()
        _, nH_err = self.get_errors()

        val_rounded, err_rounded = niutils.round_sigfigs(nH, nH_err)

        return f"{val_rounded}" + r'$\pm$' + f"{err_rounded}"


    #Returns fit model chi2 and degrees of freedom
    def get_chi2(self):
        idx_chisq = [ i for i in range(len(self.lines)) if 'Reduced' in self.lines[i] ][-1]
        chi2 = float([ self.lines[idx_chisq].split()[i+1] 
                       for i in range(len(self.lines[idx_chisq].split()))
                       if '=' in self.lines[idx_chisq].split()[i] ][0])

        dof = float([ self.lines[idx_chisq].split()[i+1]
                      for i in range(len(self.lines[idx_chisq].split()))
                      if 'for' in self.lines[idx_chisq].split()[i] ][0])

        return chi2, dof

    #Returns contour array and parameter bins
    def get_contour(self):
        #Isolate steppar runs
        idx_steppar = [ i for i in range(len(self.lines)) if 'steppar' in self.lines[i] ][-1]

        #Empty lists for data frame
        chi2 = []
        n_phot = []
        phot = []
        n_nH = []
        nH = []
        #Iterate through steppar output and append each column
        for i in range(idx_steppar+5, len(self.lines)):
            if len(self.lines[i].split()) != 7:
                break
            else:
                chi2.append(float(self.lines[i].split()[1]))
                n_phot.append(int(self.lines[i].split()[3]))
                phot.append(float(self.lines[i].split()[4]))
                n_nH.append(int(self.lines[i].split()[5]))
                nH.append(float(self.lines[i].split()[6]))

        #Range of photon indicies and column densities
        gamma_bins = list(set(phot))
        gamma_bins.sort()
        nH_bins = list(set(nH))
        nH_bins.sort()

        #Fill np 2d array
        contourarray = np.zeros((len(nH_bins), len(gamma_bins)))
        for i in range(len(chi2)):
            contourarray[n_nH[i], n_phot[i]] = chi2[i]

        #Output as named tuple
        contour_log_tuple = collections.namedtuple(
                'contour_log_tuple',
                ['array', 'gamma_bins', 'nH_bins'])

        tup = contour_log_tuple(contourarray, gamma_bins, nH_bins)
        return tup


#This class reads in the txt file output from ipl/wdata
class xspecdata:
    def __init__(self, filename):
        #Read in file and find where table breaks
        assert(os.path.isfile(filename))
        with open(filename) as h:
            lines = h.readlines()
        breakidx = np.where(np.array(lines) == 'NO NO NO NO NO\n')[0][0]
        #First table is the spectra
        df0 = pd.read_csv(filename, skiprows=3, delimiter=" ", 
                          header=None, nrows=breakidx-3)
        #Second table gives delchi
        df1 = pd.read_csv(filename, skiprows=breakidx+1, 
                          delimiter=" ", header=None)
        df0.columns = ['energy', 'energy_err', 
                       'counts', 'counts_err', 'model']
        df1.columns = ['energy', 'energy_err', 
                       'delchi', 'delchi_err', 'model']
        self.data = df0
        self.residuals = df1

        self.width = None
        self.lower = None
        self.upper = None
        self.component = None
        self.mincounts = None
        self.phot_index = None
        self.filename = filename
    def set_component(self, component):
        self.component = process.extract(component, ['primary', 'interpulse'],
                                         limit=1)[0][0]
    def set_phaserange(self, p1, p2):
        self.width=p2-p1
        self.lower = round(p1, 2)
        self.upper = round(p2, 2)

    def set_counts(self, n):
        self.mincounts = n


    def phot_index_from_log(self, fname):
        assert(os.path.isfile(fname))
        log = logfile(fname)
        self.phot_index, self.phot_index_err = niutils.round_sigfigs(
                log.get_params()[0], log.get_errors()[0])


    def get_label(self):
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            #label = f"Phase: {self.lower} -- {self.upper}"
            label = f"Phase: {self.lower}" + r'$-$' + str(self.upper)

        if self.mincounts is not None:
            label = label + f", Mincounts: {self.mincounts}"

        if self.phot_index is not None:
            label = label + r', $\Gamma=$'+'{:.2f}'.format(self.phot_index)
            if self.phot_index_err is not None:
                label += r'$\pm$' + '{:.2f}'.format(self.phot_index_err)

        return label
