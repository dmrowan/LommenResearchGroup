#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Plots phase of both mainPulse and interPulse

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
#significant-figures rounding function
def round_sf(num,sig_fig):
    if num!=0:
        return round(num, -int(math.floor(math.log10(abs(num)))-(sig_fig-1)))
    else:
        return 0
#Dom-Paper's plot formatting standards
def plotparams(ax):
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=15)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
    return ax
#command-Line Flags
parser = argparse.ArgumentParser(description ='Specify plotting parameters by c.l flags.')
parser.add_argument('mainFWHM.gauss')
parser.add_argument('mainFWHM_errors.gauss')
parser.add_argument('interFWHM.gauss')
parser.add_argument('interFWHM_errors.gauss')

parser.add_argument('-PSR','--pulsar', help=' name of pulsar', required=True)
parser.add_argument('-EBIN','--energy_binning', help=' energy_binning used', required=True)
args= vars(parser.parse_args())

#---------------Main_Pulse--------------
#input energy and main_phase data from <<mainPhase.gauss>>
energy_file = str(sys.argv[1])
energy_data= pd.read_csv(energy_file, header=None,names=['energy_bound1','energy_bound2','MP_FWHM'],delimiter=' ',engine='python')
energy_mid=(np.array(energy_data.energy_bound1)+np.array(energy_data.energy_bound2))/2
#-------PLOTTING-----
fig,ax1 = plt.subplots(1,1, figsize = (8,5))

#--chose pulsar--
if args['pulsar'] == 'B1821' :
    #plt.title('PSR B1821-24', size=25)
    ax1.text(0.067,0.94,'PSR B1821-24',size='x-large',transform=ax1.transAxes)
if args['pulsar'] == 'B1937' :
    #plt.title('PSR B1937+21', size=25)

    #ax1.set_yticks([x for x in np.arange(0.0270, 0.0300, 0.0005)])
    ax1.text(0.040,0.92,'PSR B1937+21',size='x-large',transform=ax1.transAxes)
if args['pulsar'] == 'J0218' :
    #plt.title('PSR J0218+4232', size=25)
    ax1.text(0.040,0.92,'PSR J0218+4232',size='x-large',transform=ax1.transAxes)


ax1.set_xlabel('Energy (KeV)',fontsize = 15)
ax1.set_ylabel('FWHM (phase units)',fontsize = 15)
ax1.scatter(energy_mid, energy_data.MP_FWHM, marker='x',color='purple', linewidth = 2, label= 'FWHM of P1', s = 60)
#add error-bars
main_errorFile = str(sys.argv[2])
main_error = pd.read_csv(main_errorFile,header=None,names=['errors'])
ax1.errorbar(energy_mid, energy_data.MP_FWHM,fmt='none', ecolor='purple',capzise=11, capthick = 2, yerr = main_error.errors)
#plt.show()

#apply LS-linear regression ==> plot line of best-fit
fit_slopesM,covMatrixM= np.polyfit(energy_mid,energy_data.MP_FWHM,1,w=(1/main_error.errors),cov=True)
def linearRegressionM(x):
    return x*fit_slopesM[0]+fit_slopesM[1]
ax1.plot(energy_mid,list(map(linearRegressionM,energy_mid)),'--',color='magenta',linewidth = 2)
#--chose energy range ticks--
if args['energy_binning'] == '1' :
    plt.xticks(energy_mid,('$1.0 -2.0$', '$2.0 - 3.0$', '$3.0 - 4.0$', '$4.0 - 5.0$', '$5.0 - 6.0$'),rotation=0)
if args['energy_binning'] == '0' :
    plt.xticks(energy_mid,('$0.5 -1.5$', '$1.5 - 2.5$', '$2.5 - 3.5$', '$3.5 - 4.5$', '$4.5 - 5.5$'))
if args['energy_binning'] == '8' :
    plt.xticks(energy_mid,('$0.2 -1.2$', '$1.2 - 2.2$', '$2.2 - 3.2$', '$3.2 - 4.2$', '$4.2 - 5.2$','$5.2 - 6.2$'))
#---------------Inter_Pulse---------
#input inter_phase data from <<interPhase.gauss>>
interPhase_file = str(sys.argv[3])
inter_phase = pd.read_csv(interPhase_file,header=None  ,names=['obs1','obs2','interFWHM'], delimiter=' ', engine='python')
ax1.scatter(energy_mid, inter_phase.interFWHM,marker='o',color='blue',linewidth= 2, label= 'FWHM of P2', s = 60)

#add error-bars
intererrorFile = str(sys.argv[4])
intererror = pd.read_csv(intererrorFile,header=None,names=['errors'])
print(intererror.errors)
ax1.errorbar(energy_mid,inter_phase.interFWHM,fmt='none',ecolor='blue',capsize=11 ,capthick = 2, yerr = intererror.errors)
#apply LS-linear regression ==> plot line of best-fit
fit_slopesI,covMatrixI= np.polyfit(energy_mid,inter_phase.interFWHM,1, w=(1/intererror.errors),cov=True)
def linearRegressionI(x):
    return x*fit_slopesI[0]+fit_slopesI[1]
ax1.plot(energy_mid,list(map(linearRegressionI,energy_mid)),'--',color='cyan', linewidth = 2)
#legends
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
main_slope, main_interc= scientificFormat(np.sqrt(covMatrixM[0][0]),np.sqrt(covMatrixM[1][1]))
inter_slope,inter_interc=scientificFormat(np.sqrt(covMatrixI[0][0]),np.sqrt(covMatrixI[1][1]))

custom_lines= [Line2D([0],[0],marker ='x',color='purple', lw=2, linestyle='none'),Line2D([0],[0],marker ='o', color='blue', lw=2, linestyle ='none'),Line2D([0],[0],linestyle='-',color='magenta',lw=2),Line2D([0],[0],linestyle='--',color='cyan',lw=2)]
#---------Helper_function to format legend:
def slope_format(slope,error, sf):
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
            error_fields[lindex] = str (int(error_fields[lindex])+1)
            problemIdx = lindex
            prev = problemIdx-1
            while(int(error_fields[problemIdx]) > 9 and problemIdx>=0):
                #skip '.'
                if (error_fields[prev] == '.'):
                    prev = problemIdx-2
                error_fields[problemIdx] = '0'
                error_fields[prev] = str(int(error_fields[prev])+1)
                problemIdx=prev
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

#label-equations
mainFit_label = 'Line of best fit for width of main-pulse\ny = ({:.5f} $\pm$ '.format(fit_slopesM[0])+  str(main_slope)+') x + {:.5f} $\pm$ '.format(fit_slopesM[1])+main_interc

interFit_label='Line of best fit for width of interpulse\ny = ({:.5f} $\pm$ '.format(fit_slopesI[0])+  str(inter_slope)+') x + {:.5f} $\pm$ '.format(fit_slopesI[1])+inter_interc
main_label='Slope P1 = '
inter_label='Slope P2   = '
if (fit_slopesM[0] < 0 ):
    main_label+='-'
if (fit_slopesI[0] <0):
    inter_label+='-'
main_label+=slope_format(round_sf(fit_slopesM[0],3),round_sf(np.sqrt(covMatrixM[0][0]),3),2)
inter_label+=slope_format(round_sf(fit_slopesI[0],3),round_sf(np.sqrt(covMatrixI[0][0]),3),2)
print(mainFit_label,interFit_label)

plt.legend(custom_lines,['P1 FWHM','P2 FWHM', main_label,inter_label], handletextpad=2, prop={'size':10},loc=1)
fig.tight_layout()
ax1= plotparams(ax1)
print(slope_format(0.00409,0.00109,2))
print(fit_slopesM[0],round_sf(fit_slopesM[0],2),np.sqrt(covMatrixM[0][0]),round_sf(np.sqrt(covMatrixM[0][0]),2))


plt.show()
