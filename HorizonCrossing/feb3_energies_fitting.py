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
  mintime=min(Time)
  maxtime=max(Time)
  binCounts=[]
  for i in np.arange(mintime,maxtime+binSize,binSize):
    desind=np.where((Time >= i) & (Time < i + binSize))
    binCounts.append(np.size(desind[0]))
  return np.array(binCounts)

##########################################################################################

#bin size and energy cutoffs
binSize_all = 1
lowEn = [30,100]
midEn = [100,200]
highEn = [200,1000]



#times at which each energy occurs
lowTime = enSplit(lowEn)
midTime = enSplit(midEn)
highTime = enSplit(highEn)


#number of counts in each bin of time (counts per second)
lowBinRate = countRate(lowTime,binSize_all)
midBinRate = countRate(midTime,binSize_all)
highBinRate = countRate(highTime,binSize_all)

#x-axis arrays
lowaxis=np.arange(0,len(lowBinRate),binSize_all)
midaxis=np.arange(0,len(midBinRate),binSize_all)
highaxis=np.arange(0,len(highBinRate),binSize_all)

#########################################################################################

#plot the data
#plt.plot(lowaxis,lowBinRate,'r|',label='{} keV-{} keV'.format(lowEn[0]/100,lowEn[1]/100))
#plt.plot(midaxis,midBinRate,'g|',label='{} keV-{} keV'.format(midEn[0]/100,midEn[1]/100))
#plt.plot(highaxis,highBinRate,'b|',label='{} keV-{} keV'.format(highEn[0]/100,highEn[1]/100))

#plt.xlim(450,570) #first spike
#plt.xlim(6000,6120) #second spike
#plt.xlim(11540,11660) #third spike, HORIZON CROSSING!!!
#plt.xlim(17150,17270) #fourth spike


#plt.title('Counts Per Second vs Time (V464 Sagittarius, Feb 3)')
#plt.xlabel('Time (s)')
#plt.ylabel('Counts Per Second')
#plt.legend()
#plt.show()

#########################################################################################

##fit the data##

#subsets of each time array for the 3rd spike
startTime = 390+1.92224*10**8
stopTime = 500+1.92224*10**8

def timeFit(Time):
    timeFit_spike = []
    for i in range(len(Time)):
        if ((Time[i]>=startTime)&(Time[i]<stopTime)):
            timeFit_spike.append(Time[i])
    return np.array(timeFit_spike)
#(another way to index a list is 'for indx,val in enumerate(binTime):')

#find time in fit window
lowTime_spike = timeFit(lowTime)
midTime_spike = timeFit(midTime)
highTime_spike =timeFit(highTime)

#find count rate in fit window

def removeZero(Rate_spike):
	c=Rate_spike.tolist()
	while 0 in c: c.remove(0)
	return np.array(c)

lowRate_spike = removeZero(countRate(lowTime_spike,binSize_all))
midRate_spike = removeZero(countRate(midTime_spike,binSize_all))
highRate_spike= removeZero(countRate(highTime_spike,binSize_all))

#x-axis for fits
lowaxis_fit=np.linspace(0,110,len(lowRate_spike))
midaxis_fit=np.linspace(0,110,len(midRate_spike))
highaxis_fit=np.linspace(0,110,len(highRate_spike))


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
yerrMid=np.sqrt(midRate_spike)
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

#Mid energy
popt, pcov = curve_fit(SeventhOr,midaxis_fit,midRate_spike)
plt.plot(midaxis_fit,SeventhOr(midaxis_fit,*popt),'g')
chisqMid= sum(((midRate_spike-SeventhOr(midaxis_fit,*popt))/yerrMid)**2)
print(f'The mid energy fit parameters: {popt}')
print(f'Mid Energy chi-squared: {chisqMid}')
midDiagonal = [math.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaMid = math.sqrt(sum(midDiagonal))
print(f'The overal fractional sigma is {sigmaMid}')
mid_frac_unc = paramUnc(popt,pcov,4.54309)
print(f'The overall fractional uncertainty is {mid_frac_unc}')
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
ax.errorbar(midaxis_fit,midRate_spike,fmt='g.',yerr=yerrMid,label='{} keV-{} keV'.format(midEn[0]/100,midEn[1]/100))
ax.errorbar(highaxis_fit,highRate_spike,fmt='b.',yerr=yerrHigh,label='{} keV-{} keV'.format(highEn[0]/100,highEn[1]/100))

#display
plt.title('Counts Per Second vs Time (V464 Sagittarius, Feb 3)')
plt.xlabel('Time (s)')
plt.ylabel('Counts Per Second')
plt.legend()
plt.show()


'''
def poptArray(i,sigma):
	return [popt[i]+sigma[i],popt[i]-sigma[i]]

t0=poptArray(0,sigmaLow)
t1=poptArray(1,sigmaLow)
t2=poptArray(2,sigmaLow)
t3=poptArray(3,sigmaLow)
t4=poptArray(4,sigmaLow)
t5=poptArray(5,sigmaLow)
t6=poptArray(6,sigmaLow)
t7=poptArray(7,sigmaLow)

combos=list(itertools.product(t0,t1,t2,t3,t4,t5,t6,t7))

values=[]
for i in range(len(combos)):
	a,b,c,d,e,f,g,h=combos[i]
	values.append(SeventhOr(lowaxis,a,b,c,d,e,f,g,h))

values=np.array(values)
fitError=np.std(values,axis=0)
'''
