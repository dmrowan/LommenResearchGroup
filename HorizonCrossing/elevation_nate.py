#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy import interpolate

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
    binCounts.append(np.size(desind[0]))
    altitude.append(np.mean(alt_array[desind[0]]))
  return np.array(binCounts),np.array(altitude)

#function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return eventTime[index[0]]


class EnergyBands:

  def __init__(self,energy_band,bin_size):
    self.energy_band = energy_band
    self.bin_size = bin_size
    self.time = enSplit(energy_band)
    self.alt = altSplit(energy_band)
    self.rate,self.new_alt = countRate(self.time,self.alt,bin_size)


low_en = EnergyBands(lowEn,binSize_all)
mid_en = EnergyBands(midEn,binSize_all)
high_en = EnergyBands(highEn,binSize_all)

################################################################################
#plot the data
plt.figure(1)
plt.plot(low_en.new_alt,low_en.rate,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
plt.plot(mid_en.new_alt,mid_en.rate,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
plt.plot(high_en.new_alt,high_en.rate,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
plt.title("Counts/Sec vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("X-Ray Photon Counts/Sec")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(low_en.new_alt,(low_en.rate/max(low_en.rate))*100,'r--',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
plt.plot(mid_en.new_alt,(mid_en.rate/max(mid_en.rate))*100,'g--',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
plt.plot(high_en.new_alt,(high_en.rate/max(high_en.rate))*100,'b--',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
plt.title("Percent Transmittance vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("Percent Transmittance of X-Rays")
plt.legend()
plt.show()

