#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Plots phase separation of both mainPulse and interPulse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize':11})
import matplotlib.pyplot as plt
from pylab import show,xlabel,ylabel
import matplotlib.ticker as mticker
import sys
from scipy.optimize import curve_fit
from scipy import stats
from scipy.interpolate import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import argparse
import math
from matplotlib import rc
rc('text', usetex=True)
#Dom-Paper's plot formatting standards
def plotparams(ax):
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=15)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
    return ax
# My slope formatting standards lol
def scientificFormat(unit1, unit2):
    f= mticker.ScalarFormatter(useOffset=False,useMathText=True)
    g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
    fmt = mticker.FuncFormatter(g)
    err1 = str(format(fmt(unit1)))
    err2 = str(format(fmt(unit2)))
    a,b = err1.split('{\\times}')
    l,a= a.split("$")
    a = '{:.2f}'.format(float(a))
    err1= '$'+str(a)+'{\\times}'+b
    c,d=err2.split('{\\times}')
    f,c= c.split("$")
    c= '{:.2f}'.format(float(c))
    err2= '$'+str(c)+'{\\times}'+d
    return err1,err2
def slope_format(slope,error,sf):
    #slope_fields = list(str(slope))
    slope_fields=list('{:.11f}'.format(round_sf(slope,sf)))
    error_fields=list('{:.11f}'.format(round_sf(error,2)))
    f=len(slope_fields)-1
    while ( f>=0 and slope_fields[f]== '0' ):
        slope_fields.pop(f)
        f-=1
    if '-' in slope_fields :
        slope_fields.remove('-')
    if '-' in error_fields:
        error_fields.remove('-')
    lindex=0
    findex=0
    i=0
    while i< len(slope_fields):
        if (slope_fields[i] !='0'):
            lindex=i
        i+=1
    if (len(error_fields) >= (lindex+2)):
        if (int(error_fields[lindex+1]) >= 5):
            error_fields[lindex] = str (int(error_fields[lindex+1])+1)
    if (len(error_fields) <= lindex):
        times = lindex+1-len(error_fields)
        k=0
        while(k < times):
            error_fields.append('0')
            k+=1
    j=0
    while j< len(error_fields):
        if ((error_fields[j] != '0') and (error_fields[j] != '.') ):
            findex=j
            break
        j+=1
    if (lindex < findex):
        onces = len(error_fields)-lindex-1
        m=0
        while( m<onces):
            slope_fields.append('0')
            m+=1
        lindex= lindex+onces

    digits = error_fields[findex:lindex+1]
    answer = ''.join(slope_fields)+'('+''.join(digits)+')'
    if (len(digits) >= 3 ):
        return slope_format(slope,error,1)
    return answer
def round_sf(num,sig_fig):
    if num!=0:
        return round(num, -int(math.floor(math.log10(abs(num)))-(sig_fig-1)))
    else:
        return 0
#----MAIN----
#MP phase
energy_file = str(sys.argv[1])
energy_data= pd.read_csv(energy_file, header=None,names=['energy_bound1','energy_bound2','MP_phase'],delimiter=' ',engine='python')
energy_mid=(np.array(energy_data.energy_bound1,dtype='float')+np.array(energy_data.energy_bound2,dtype='float'))/2
#MP error
main_errorFile = str(sys.argv[2])
main_error = pd.read_csv(main_errorFile,header=None,names=['errors'])
#IP phase
interPhase_file = str(sys.argv[3])
inter_phase = pd.read_csv(interPhase_file,header=None  ,names=['obs1','obs2','interPhase'], delimiter=' ', engine='python')
#IP error
intererrorFile = str(sys.argv[4])
intererror = pd.read_csv(intererrorFile,header=None,names=['errors'])
#ANALYSIS
pulse_separation = np.array(energy_data.MP_phase,dtype='float')-np.array(inter_phase.interPhase,dtype='float')
MP_err = np.array(main_error.errors,dtype='float')
IP_err = np.array(intererror.errors,dtype='float')
propag_err = np.sqrt(MP_err**2+IP_err**2)
#PLOTTING
fig,ax1 = plt.subplots(1,1, figsize = (8,5))
plt.errorbar(energy_mid,pulse_separation,fmt='o',elinewidth=2,capsize=10,capthick=2,yerr=propag_err,color='xkcd:violet')
ax1.text(0.040,0.92,'PSR B1821'+r'$-$'+'24',size='x-large',transform=ax1.transAxes)
plt.xticks(energy_mid,('$1.0 -2.0$', '$2.0 - 3.0$', '$3.0 - 4.0$', '$4.0 - 5.0$', '$5.0 - 6.0$'),rotation=0)
ax1.set_xlabel('Energy (KeV)',fontsize = 15)
ax1.set_ylabel('Pulse Separation',fontsize = 15)
ax1 = plotparams(ax1)
#LINEAR REGRESSION
def line(x,a,b):
    return a*x+b
#perform weighted linear regression
opt,cov = curve_fit(line,energy_mid,pulse_separation,sigma=propag_err)
line_label='Slope = '
slope = opt[0]
error = np.sqrt(cov[0][0])
if (slope < 0 ):
    line_label+=r'$-$'
unit1,unit2=scientificFormat(slope,error)
line_label+= slope_format(slope,error,2)
print(line_label)
plt.plot(energy_mid,line(energy_mid,*opt),color='xkcd:magenta',linestyle='--',label=line_label)
custom_lines= [Line2D([0],[0],marker ='o',color='xkcd:violet', lw=2, linestyle='none'),Line2D([0],[0],color='xkcd:magenta', lw=2, linestyle='--')]
plt.legend(custom_lines,['Pulse separation',line_label],edgecolor='black',loc=1)

plt.show()
