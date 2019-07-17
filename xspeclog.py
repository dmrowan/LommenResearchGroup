#!/usr/bin/env python
import ast
import collections
import numpy as np
import os
import pandas as pd
from fuzzywuzzy import process

import niutils

#Dom Rowan 2019
desc="""
Classes for parsing log and wdata files from Xspec

logfile: Handles collection of parameters, errors, and contours
"""

class modelparam:
    def __init__(self, value=None, error=None, num=None, unit=None):
        self.value = value
        self.error = error
        self.num = num
        self.unit= unit
    
    def set_name(self, name):
        self.name = name
        self.label = name

    def set_label(self, label):
        self.label = label


#Class for parsing the log output of Xspec
class logfile:
    
    #Initialize with path to log file
    def __init__(self, fname):

        if not os.path.isfile(fname):
            raise FileNotFoundError

        self.fname = fname
        self.lines = open(self.fname, 'r').readlines()
        self.params = None

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
        phot_index_val = float([ v[0] for v in [ [ self.lines[j].split()[i+1]
                                               for i in range(len(self.lines[j].split()))
                                               if 'PhoIndex' in self.lines[j].split()[i] ]
                                            for j in range(idx_params, idx_params+5) ]
                             if len(v) != 0 ][0])


        phot_index_number = int([ v[0] for v in [ [ self.lines[j].split()[i-3]
                                                    for i in range(len(self.lines[j].split()))
                                                    if 'PhoIndex' in self.lines[j].split()[i] ]
                                                 for j in range(idx_params, idx_params+5) ]
                                  if len(v) != 0 ][0])

        self.phot_index = modelparam(value=phot_index_val, num=phot_index_number)
        self.phot_index.set_name("phot_index")
        self.phot_index.set_label(r'Photon Index $\mathrm{\Gamma}$')
                                                   
                      
        #Find the Hydrogen column density
        nH_val = float([ v[0] for v in [ [self.lines[j].split()[i+1]
                                      for i in range(len(self.lines[j].split()))
                                      if '10^22' in self.lines[j].split()[i] ]
                                    for j in range(idx_params, idx_params+5) ]
                     if len(v) != 0 ][0])

        nH_number = int([ v[0] for v in [ [self.lines[j].split()[i-4]
                                           for i in range(len(self.lines[j].split()))
                                           if '10^22' in self.lines[j].split()[i] ]
                                         for j in range(idx_params, idx_params+5) ]
                          if len(v) != 0 ][0])
        
        self.nH = modelparam(value=nH_val, num=nH_number, unit='10^22')
        self.nH.set_name("nH")
        self.nH.set_label(r'Column Density $N_H$ $(10^{22}$ cm$^{-2}$)')
        
        #Find normalization 
        norm_val = float([ v[0] for v in [ [ self.lines[j].split()[i+1]
                                         for i in range(len(self.lines[j].split()))
                                         if 'norm' in self.lines[j].split()[i] ]
                                      for j in range(idx_params, idx_params+5) ]
                       if len(v) != 0 ][0])

        norm_number = int([ v[0] for v in [ [ self.lines[j].split()[i-3]
                                              for i in range(len(self.lines[j].split()))
                                              if 'norm' in self.lines[j].split()[i] ]
                                           for j in range(idx_params, idx_params+5) ]
                            if len(v) != 0 ][0])
       
        self.norm = modelparam(value=norm_val, num=norm_number)
        self.norm.set_name("norm")
        self.norm.set_label("Normalization")
        
        self.params = {'phot_index':self.phot_index, 
                     'nH':self.nH,
                     'norm':self.norm}

        return [ self.params[k].value for k, p in self.params.items() ]

    #Returns errors on photon index and column density
    def get_errors(self):
        if self.params is None:
            self.get_params()

        #Find the confidence interval on the parameters
        idx_conf = [ i for i in range(len(self.lines))
                     if 'Confidence Range' in self.lines[i] ]
        if len(idx_conf) == 0:
            print("No error calculation found in logfile")
        else:
            idx_conf = idx_conf[-1]

    
        err_list = [ [int(self.lines[i].strip("#").split()[0]),
                      np.mean( [ abs(val) 
                                 for val in ast.literal_eval(
                                 self.lines[i].strip("#").split()[-1]) ])] 
                      for i in range(idx_conf+1, idx_conf+len(self.params)) 
                      if 'Warning' not in self.lines[i] ]
        
        for pair in err_list:
            for key, param in self.params.items():
                if param.num == pair[0]:
                    param.error = pair[1]

        return [ self.params[k].error for k, p in self.params.items() ]
        
    def param_string(self, name):
        self.get_params()
        self.get_errors()

        name = process.extract(name, list(self.params.keys()), 
                               limit=1)[0][0]

        val_rounded, err_rounded = niutils.round_sigfigs(self.params[name].value,
                                                       self.params[name].error)

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
        self.get_params()
        #Isolate steppar runs
        idx_steppar = [ i for i in range(len(self.lines)-1) if 'Delta' in self.lines[i]
                        and 'Chi-Squared' in self.lines[i+1] ][-1]

        params_in_contour = [ int(self.lines[idx_steppar+1].split()[-2]),
                              int(self.lines[idx_steppar+1].split()[-1]) ]

        param_x_label = [ param.label for key, param in self.params.items() 
                          if param.num == params_in_contour[0] ][0]
        param_y_label = [ param.label for key, param in self.params.items()
                          if param.num == params_in_contour[1] ][0]

        #Empty lists for data frame
        chi2 = []
        n_param_x = []
        param_x = []
        n_param_y = []
        param_y = []

        #Iterate through steppar output and append each column
        for i in range(idx_steppar+3, len(self.lines)):
            if len(self.lines[i].strip("#").split()) != 6:
                break
            else:
                chi2.append(float(self.lines[i].split()[0]))
                n_param_x.append(int(self.lines[i].split()[2]))
                param_x.append(float(self.lines[i].split()[3]))
                n_param_y.append(int(self.lines[i].split()[4]))
                param_y.append(float(self.lines[i].split()[5]))

        #Range of parameter values
        param_x_bins = list(set(param_x))
        param_x_bins.sort()
        param_y_bins = list(set(param_y))
        param_y_bins.sort()


        #Fill np 2d array
        contourarray = np.zeros((len(param_y_bins), len(param_x_bins)))
        for i in range(len(chi2)):
            contourarray[n_param_y[i], n_param_x[i]] = chi2[i]

        #Output as named tuple
        contour_log_tuple = collections.namedtuple(
                'contour_log_tuple',
                ['array', 'param_x_bins', 'param_y_bins',
                 'xlabel', 'ylabel'])

        tup = contour_log_tuple(contourarray, param_x_bins, param_y_bins,
                                param_x_label, param_y_label)
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
