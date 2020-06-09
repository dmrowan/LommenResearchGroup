#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.optimize import curve_fit
import itertools
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

fig,ax= plt.subplots()
yerrLow=np.sqrt(lowRate_spike)
yerrMid=np.sqrt(midRate_spike)
yerrHigh=np.sqrt(highRate_spike)

def paramUnc(paramIn,Popt,Pcov):
	Popt.tolist()
	fVal= SeventhOr(startTime,*Popt)
	Popt[paramIn]=Popt[paramIn]+np.sqrt(abs(pcov[paramIn][paramIn]))
	fNew=SeventhOr(startTime,*Popt)
	frac_unc=abs(fNew-fVal)/fVal
	return (frac_unc)**2


#defines the curve fit and plots the model over the data

popt, pcov = curve_fit(SeventhOr,lowaxis_fit,lowRate_spike)
print("The low energy popt array is",popt)
plt.plot(lowaxis_fit,SeventhOr(lowaxis_fit,*popt),'r-')
chisqLow= sum(((lowRate_spike-SeventhOr(lowaxis_fit,*popt))/yerrLow)**2)
print("The low energy chi-squared is",chisqLow)
lowDiagonal= [pcov[i][i] for i in range(len(np.diagonal(pcov)))]
sigmaLow=np.sqrt(sum(lowDiagonal))
print("The overall fractional sigma is",sigmaLow)


low_unc_a=paramUnc(0,popt,pcov)
low_unc_b=paramUnc(1,popt,pcov)
low_unc_c=paramUnc(2,popt,pcov)
low_unc_d=paramUnc(3,popt,pcov)
low_unc_e=paramUnc(4,popt,pcov)
low_unc_f=paramUnc(5,popt,pcov)
low_unc_g=paramUnc(6,popt,pcov)
low_unc_h=paramUnc(7,popt,pcov)

low_frac_unc=np.sqrt(low_unc_a+low_unc_b+low_unc_c+low_unc_d+low_unc_e+low_unc_f+low_unc_g+low_unc_h)
print("The overal fractional uncertainity is", low_frac_unc)

print("----------------------")


popt, pcov= curve_fit(SeventhOr,midaxis_fit,midRate_spike)
print("The mid energy popt array is",popt)
plt.plot(midaxis_fit,SeventhOr(midaxis_fit,*popt),'g-')
chisqMid= sum(((midRate_spike-SeventhOr(midaxis_fit,*popt))/yerrMid)**2)
print("The mid energy chi-squared is",chisqMid)
midDiagonal= [pcov[i][i] for i in range(len(np.diagonal(pcov)))]
sigmaMid=np.sqrt(sum(midDiagonal))
print("The overall fractional sigma is",sigmaMid)

med_unc_a=paramUnc(0,popt,pcov)
med_unc_b=paramUnc(1,popt,pcov)
med_unc_c=paramUnc(2,popt,pcov)
med_unc_d=paramUnc(3,popt,pcov)
med_unc_e=paramUnc(4,popt,pcov)
med_unc_f=paramUnc(5,popt,pcov)
med_unc_g=paramUnc(6,popt,pcov)
med_unc_h=paramUnc(7,popt,pcov)

med_frac_unc=np.sqrt(med_unc_a+med_unc_b+med_unc_c+med_unc_d+med_unc_e+med_unc_f+med_unc_g+med_unc_h)
print("The overal fractional uncertainity is", med_frac_unc)

print("----------------------")



popt, pcov= curve_fit(SeventhOr,highaxis_fit,highRate_spike)
print("The high energy popt array is",popt)
plt.plot(highaxis_fit,SeventhOr(highaxis_fit,*popt),'b-')
chisqHigh= sum(((highRate_spike-SeventhOr(highaxis_fit,*popt))/yerrHigh)**2)
print("The high energy chi-squared is",chisqHigh)
highDiagonal= [np.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaHigh=np.sqrt(sum(highDiagonal))
print("The overall fractional sigma is",sigmaHigh)

high_unc_a=paramUnc(0,popt,pcov)
high_unc_b=paramUnc(1,popt,pcov)
high_unc_c=paramUnc(2,popt,pcov)
high_unc_d=paramUnc(3,popt,pcov)
high_unc_e=paramUnc(4,popt,pcov)
high_unc_f=paramUnc(5,popt,pcov)
high_unc_g=paramUnc(6,popt,pcov)
high_unc_h=paramUnc(7,popt,pcov)

high_frac_unc=np.sqrt(high_unc_a+high_unc_b+high_unc_c+high_unc_d+high_unc_e+high_unc_f+high_unc_g+high_unc_h)
print("The overal fractional uncertainity is", high_frac_unc)


print("----------------------")


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
