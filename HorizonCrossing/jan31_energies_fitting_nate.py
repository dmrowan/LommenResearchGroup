#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.optimize import curve_fit
import itertools
import math
#########################################################################################

#reading the data table
tab = Table.read('cleanfilt.evt',hdu=1)
tab2 = Table.read('ni2200300101.mkf',hdu=1)
elevArray = np.array(tab2['ELV'])
timeArray = np.array(tab['TIME'])
enArray = np.array(tab['PI'])

#########################################################################################

#function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return timeArray[index[0]]

#########################################################################################

#function that deduces the number of counts per bin size 
def countRate(Time,binSize):
  binCounts=[]
  for i in np.arange(min(Time),max(Time)+binSize,binSize):
    desind=np.where((Time >= i) & (Time < i + binSize))
    binCounts.append(np.size(desind[0]))
  return np.array(binCounts)

##########################################################################################

#bin size and energy cutoffs
binSize_all = 1
lowEn = [35,200]
highEn = [200,800]


#times at which each energy occurs
lowTime = enSplit(lowEn)
highTime = enSplit(highEn)


#number of counts in each bin of time (counts per second)
lowBinRate = countRate(lowTime,binSize_all)
highBinRate = countRate(highTime,binSize_all)

#x-axis arrays
lowaxis=np.arange(0,len(lowBinRate),binSize_all)
highaxis=np.arange(0,len(highBinRate),binSize_all)

#########################################################################################

#plot the data
plt.plot(lowaxis,lowBinRate,'r|',label='{} keV-{} keV'.format(lowEn[0]/100,lowEn[1]/100))
plt.plot(highaxis,highBinRate,'b|',label='{} keV-{} keV'.format(highEn[0]/100,highEn[1]/100))

#plt.xlim(0,400) #first spike
plt.xlim(5650,5800) #second spike


plt.title('Counts Per Second vs Time (V464 Sagittarius, Jan 31)')
plt.xlabel('Time (s)')
plt.ylabel('Counts Per Second')
plt.legend()
plt.show()

print(elevArray)

#########################################################################################
'''
##fit the data##

#subsets of each time array for the 1st spike
startTime = 5650
stopTime = 5800

def timeFit(Time):
    timeFit_spike = []
    for i in range(len(Time)):
        if ((Time[i]>=startTime)&(Time[i]<stopTime)):
            timeFit_spike.append(Time[i])
    return np.array(timeFit_spike)
#(another way to index a list is 'for indx,val in enumerate(binTime):')

#find time in fit window
lowTime_spike = timeFit(lowTime)
highTime_spike = timeFit(highTime)

#find count rate in fit window

#def removeZero(Rate_spike):
#	c=Rate_spike.tolist()
#	while 0 in c: c.remove(0)
#	return np.array(c)

lowRate_spike = countRate(lowTime_spike,binSize_all)
highRate_spike = countRate(highTime_spike,binSize_all)

#x-axis for fits
lowaxis_fit=np.linspace(0,400,len(lowRate_spike))
highaxis_fit=np.linspace(0,400,len(highRate_spike))


#defines function for nth-Order Polynomial
def SecOr(x,a,b,c):
        return(a*x**2+b*x+c)

def ThirdOr(x,a,b,c,d):
	return(a*x**3+b*x**2+c*x+d)

def FourthOr(x,a,b,c,d,e):
	return(a*x**4+b*x**3+c*x**2+d*x+e)

def FifthOr(x,a,b,c,d,e,f):
	return(a*x**5+b*x**4+c*x**3+d*x**2+e*x+f)

def SeventhOr(x,a,b,c,d,e,f,g,h):
	return(a*x**7+b*x**6+c*x**5+d*x**4+e*x**3+f*x**2+g*x+h)


fig,ax = plt.subplots()
yerrLow=np.sqrt(lowRate_spike)
yerrHigh=np.sqrt(highRate_spike)

#plot the trendlines and analyze error

#new method to calculate uncertainty
def paramUnc(Popt,Pcov,Xint):
    Popt.tolist()
    fVal = SeventhOr(startTime+Xint,*Popt)
    frac_unc_params = []
    added_frac_unc = 0

    for paramIn in range(len(Popt)):
      Popt[paramIn]=Popt[paramIn]+math.sqrt(abs(pcov[paramIn][paramIn]))
      fNew = SeventhOr(startTime+Xint,*Popt)
      frac_unc = abs(fNew-fVal)/fVal
      frac_unc_params.append((frac_unc)**2)

    for i in range(len(frac_unc_params)):
      added_frac_unc += frac_unc_params[i]

    return math.sqrt(added_frac_unc) 


#Low energy
popt, pcov = curve_fit(SeventhOr,lowaxis_fit,lowRate_spike)
plt.plot(lowaxis_fit,SeventhOr(lowaxis_fit,*popt),'r')
chisqLow = sum(((lowRate_spike-SeventhOr(lowaxis_fit,*popt))/yerrLow)**2)
print(f'The low energy fit parameters: {popt}(highest power first)')
print(f'Low Energy chi-squared: {chisqLow}')
lowDiagonal = [math.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaLow = math.sqrt(sum(lowDiagonal))
print(f'The overall fractional sigma is {sigmaLow}')
low_frac_unc = paramUnc(popt,pcov,9.57002)
print(f'The overall fractional uncertainty is {low_frac_unc}')
print('....')

#High energy
popt, pcov = curve_fit(SeventhOr,highaxis_fit,highRate_spike)
plt.plot(highaxis_fit,SeventhOr(highaxis_fit,*popt),'b')
chisqHigh= sum(((highRate_spike-SeventhOr(highaxis_fit,*popt))/yerrHigh)**2)
print(f'The high enery fit parameters: {popt}')
print(f'High Energy chi-squared: {chisqHigh}')
highDiagonal = [math.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaHigh = math.sqrt(sum(highDiagonal))
print(f'The overall fractional sigma is {sigmaHigh}')
high_frac_unc = paramUnc(popt,pcov,2.13372)
print(f'The overall fractional uncertainty is {high_frac_unc}')
print('....')

#plot the data with errorbars
ax.errorbar(lowaxis_fit,lowRate_spike,fmt='r.',yerr=yerrLow,label='{} keV-{} keV'.format(lowEn[0]/100,lowEn[1]/100))
ax.errorbar(highaxis_fit,highRate_spike,fmt='b.',yerr=yerrHigh,label='{} keV-{} keV'.format(highEn[0]/100,highEn[1]/100))

#display
plt.title('Counts Per Second vs Time (V464 Sagittarius, Feb 3)')
plt.xlabel('Time (s)')
plt.ylabel('Counts Per Second')
plt.legend()
plt.show()
'''