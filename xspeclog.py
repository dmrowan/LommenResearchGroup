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

def warning_filter(line):
    if "Parameter pegged" in line:
        return
    else:
        print(line)


class modelcomp:

    def __init__(self, name, param):

        self.name = name
        self.params = { param.name:param }

    def add_param(self, new):

        self.params[new.name] = new

    def set_num(self, component_number):
        self.n = component_number


class modelparam:
    def __init__(self, name=None, value=None, error=None, 
                 param_num=None, comp_num=None, unit=None, frozen=False):
        self.name = name
        self.value = value
        self.error = error
        self.param_num = param_num
        self.comp_num = comp_num
        self.unit= unit
        self.frozen = frozen
    
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
        self.model = {}

    #Returns photon index and column density values
    def get_params(self):

        #Find line where fit ended
        idx_params = [ i for i in range(len(self.lines))
                       if 'par  comp' in self.lines[i] ]
        if len(idx_params) == 0:
            return "No parameter information found in logfile"
        else:
            idx_params = idx_params[-1]


        iiter = idx_params+1
        while "___________" not in self.lines[iiter]:
            current_row = self.lines[iiter].strip('#').split()
            parnumber = int(current_row[0])
            compnumber = int(current_row[1])
            compname = current_row[2]
            paramname = current_row[3]
            try:
                value = float(current_row[4])
                unit = None
            except:
                unit = current_row[4]
                value = float(current_row[5])

            if current_row[-1] == 'frozen':
                frozen=True
            else:
                frozen=False

            param = modelparam(name=paramname, 
                               param_num=parnumber, comp_num=compnumber, 
                               value=value, unit=unit, frozen=frozen)

            if compname in self.model:
                if param not in self.model[compname].params:
                    self.model[compname].add_param(param)
            else:
                self.model[compname] = modelcomp(compname, param)
                self.model[compname].set_num(compnumber)

            iiter += 1

        return self.model

    #Returns errors on photon index and column density
    def get_errors(self, print_warnings=True):
        if self.model is None:
            self.get_params()

        #Find the confidence interval on the parameters
        idx_conf = [ i for i in range(len(self.lines))
                     if 'Confidence Range' in self.lines[i] ]
        if len(idx_conf) == 0:
            print("No error calculation found in logfile")
        else:
            idx_conf = idx_conf[-1]


        iiter = idx_conf+1
        while ( (self.lines[iiter].strip('#') != '\n') and
                ( not self.lines[iiter].strip('#').startswith('XSPEC'))):
            if not self.lines[iiter].strip('#').split()[0].isdigit():
                if print_warnings:
                    warning_filter(self.lines[iiter].strip('#').strip('\n'))
            else:
                param_number = int(self.lines[iiter].strip('#').split()[0])   
                error = np.mean( [ abs(val) for val in ast.literal_eval(
                                   self.lines[iiter].strip("#").split()[-1]) ] )

                for comp in self.model.values():
                    for param in comp.params.values():
                        if param.param_num == param_number:
                            param.error = error

            iiter+=1

        
    def param_string(self, comp, param):
        self.get_params()
        self.get_errors()

        comp = process.extract(comp, list(self.model.keys()),
                              limit=1)[0][0]
        param = process.extract(param, list(self.model[comp].params.keys()),
                                limit=1)[0][0]

        value = self.model[comp].params[param].value
        error = self.model[comp].params[param].error
        if error is not None:
            val_rounded, err_rounded = niutils.round_sigfigs(value, error)
            return f"{val_rounded}" + r'$\pm$' + f"{err_rounded}"
        else:
            return f"{value}"

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


    def get_commands(self):

        commands = [ self.lines[i] for i in range(len(self.lines))
                     if self.lines[i][0] == '!' ]

        commands = [ c.replace('!XSPEC12>', '').lstrip() 
                     for c in commands ]

        return commands

#This class reads in the txt file output from ipl/wdata
class xspecdata:
    def __init__(self, filename):
        #Read in file and find where table breaks
        assert(os.path.isfile(filename))
        with open(filename) as h:
            lines = h.readlines()

        try:
            breakidx = np.where(np.array(lines) == 'NO NO NO NO NO\n')[0][0]
            irregular=False
        except:
            irregular=True
            breakidx = np.flatnonzero(
                    np.core.defchararray.find(np.array(lines), 'NO NO NO NO NO\n')!=-1)[0]
        #First table is the spectra
        df0 = pd.read_csv(filename, skiprows=3, delimiter=" ", 
                          header=None, nrows=breakidx-3)
        #Second table gives delchi
        df1 = pd.read_csv(filename, skiprows=breakidx+1, 
                          delimiter=" ", header=None)
        if not irregular:
            df0.columns = ['energy', 'energy_err', 
                           'counts', 'counts_err', 'model']
            df1.columns = ['energy', 'energy_err', 
                           'delchi', 'delchi_err', 'model']
        else:
            df0.columns = ['energy', 'energy_err', 
                           'counts', 'counts_err', 'model', '5', '6']
            df1.columns = ['energy', 'energy_err', 
                           'delchi', 'delchi_err', 'model', '5', '6']
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
        self.phot_string = log.param_string("powerlaw", 'phot')

    def get_label(self):
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            #label = f"Phase: {self.lower} -- {self.upper}"
            label = f"Phase: {self.lower}" + r'$-$' + str(self.upper)

        if self.mincounts is not None:
            label = label + f", Mincounts: {self.mincounts}"

        if self.phot_string is not None:
            label += r', $\Gamma=$ ' + self.phot_string

        return label
