#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import scipy.signal
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import signal

#Time range around the horizon crossing
startTime = 390+1.92224*10**8
stopTime = 500+1.92224*10**8

startTimeIndex = 311883
stopTimeIndex = 352360

#read in the data files

tab_ni = Table.read('ni2200300102.mkf',hdu=1)
timeArray = np.array(tab_ni['TIME'])
elevArray = np.array(tab_ni['ELV'])
enArray_low = np.array(tab_ni['FPM_XRAY_PI_0035_0200'])
enArray_mid = np.array(tab_ni['FPM_XRAY_PI_0800_1200'])

tab_evt = Table.read('cleanfilt.evt',hdu=1)
eventTime = np.array(tab_evt['TIME'][startTimeIndex:stopTimeIndex])
enArray = np.array(tab_evt['PI'][startTimeIndex:stopTimeIndex])

#bin size and energy band cutoffs
binSize_all = 1
lowEn = [30,100]
midEn = [100,200]
highEn = [200,1000]

#interpolate the times.evt to go over the range of elevations.mkf
f = interpolate.interp1d(timeArray,elevArray,kind='linear')
elev_evt = f(eventTime)

#calculate altitude based on elevation angle
R = 6378
H = 410
theta = np.arcsin(R/(R+H))
altArray = []
for val in elev_evt:
  altArray.append((R+H)*np.sin(theta+val*(np.pi/180))-R)
altArray=np.array(altArray)

#function that splits the altitudes based on energy bands
def altSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return altArray[index[0]]

#function that deduces the number of counts per bin size 
def countRate(Time,alt_array,binSize):
  binCounts=[]
  altitude = []
  for i in np.arange(min(Time),max(Time)+binSize,binSize):
    desind=np.where((Time >= i) & (Time < i + binSize))
    if len(alt_array[desind[0]]) !=0:
      binCounts.append(np.size(desind[0]))
      altitude.append(np.mean(alt_array[desind[0]]))
  return np.array(binCounts),np.array(altitude)

#function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return eventTime[index[0]]

def percTrans(Alt,Rate):
  plateau = np.where(((Alt>200)&(Alt<250)))
  avg = np.mean(Rate[plateau[0]])
  return (Rate/avg)*100

def SeventhOr(x,a,b,c,d,e,f,g,h):
	return(a*x**7+b*x**6+c*x**5+d*x**4+e*x**3+f*x**2+g*x+h)

def curveFit(newAlt,Rate):
  popt,pcov = curve_fit(SeventhOr,newAlt,Rate)
  return popt,pcov

class EnergyBands:

  def __init__(self,energy_band,bin_size):
    self.energy_band = energy_band
    self.bin_size = bin_size
    self.time = enSplit(energy_band)
    self.alt = altSplit(energy_band)
    self.rate,self.new_alt = countRate(self.time,self.alt,bin_size)
    self.perc_trans = percTrans(self.new_alt,self.rate)
    self.popt,self.pcov = curveFit(self.new_alt,self.rate)
    self.popt_perc,self.pcov_perc = curveFit(self.new_alt,self.perc_trans)

low_en = EnergyBands(lowEn,binSize_all)
mid_en = EnergyBands(midEn,binSize_all)
high_en = EnergyBands(highEn,binSize_all)

################################################################################
#plot the data

plt.figure(1)
plt.plot(low_en.new_alt,low_en.rate,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
#plt.plot(low_en.new_alt,SeventhOr(low_en.new_alt,*low_en.popt),'r-')

plt.plot(mid_en.new_alt,mid_en.rate,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
#plt.plot(mid_en.new_alt,SeventhOr(mid_en.new_alt,*mid_en.popt),'g-')

plt.plot(high_en.new_alt,high_en.rate,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
#plt.plot(high_en.new_alt,SeventhOr(high_en.new_alt,*high_en.popt),'b-')

plt.title("Counts/Sec vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("X-Ray Photon Counts/Sec")
plt.legend()
plt.show()

plt.figure(2)

plt.plot(low_en.new_alt,low_en.perc_trans,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
#plt.plot(low_en.new_alt,SeventhOr(low_en.new_alt,*low_en.popt_perc),'r-')
plt.plot(low_en.new_alt,scipy.signal.savgol_filter(low_en.perc_trans,51,3),'r--')

plt.plot(mid_en.new_alt,mid_en.perc_trans,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
#plt.plot(mid_en.new_alt,SeventhOr(mid_en.new_alt,*mid_en.popt_perc),'g-')
plt.plot(mid_en.new_alt,scipy.signal.savgol_filter(mid_en.perc_trans,51,3),'g--')

plt.plot(high_en.new_alt,high_en.perc_trans,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
#plt.plot(high_en.new_alt,SeventhOr(high_en.new_alt,*high_en.popt_perc),'b-')
plt.plot(high_en.new_alt,scipy.signal.savgol_filter(high_en.perc_trans,51,3),'b--')


plt.title("Percent Transmittance vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("Percent Transmittance of X-Rays")
plt.legend()
plt.show()

#################################################################
