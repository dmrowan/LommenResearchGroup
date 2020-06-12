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
def countRate(Time,binSize):
  binCounts=[]
  for i in np.arange(min(Time),max(Time)+binSize,binSize):
    desind=np.where((Time >= i) & (Time < i + binSize))
    binCounts.append(np.size(desind[0]))
  return np.array(binCounts)

#function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return eventTime[index[0]]

#time axis from count rate vs time plot
def createAxis(alt_split,bin_size):
  return np.arange(0,len(alt_split),bin_size)

#interpolate the count rate to match the range of elevations
def newAlt(Axis,oldAlt,CountRate):
  g = interpolate.interp1d(Axis,oldAlt,kind='linear')
  return g(np.arange(0,len(CountRate),binSize_all))

class EnergyBands:

  def __init__(self,energy_band,bin_size):
    self.energy_band = energy_band
    self.bin_size = bin_size
    self.time = enSplit(energy_band)
    self.rate = countRate(self.time,bin_size)
    self.alt = altSplit(energy_band)
    self.axis = createAxis(self.alt,bin_size)
    self.new_alt = newAlt(self.axis,self.alt,self.rate)


low_en = EnergyBands(lowEn,binSize_all)
mid_en = EnergyBands(midEn,binSize_all)
high_en = EnergyBands(highEn,binSize_all)

################################################################################
#plot the data

plt.plot(low_en.new_alt,low_en.rate,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
plt.plot(mid_en.new_alt,mid_en.rate,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
plt.plot(high_en.new_alt,high_en.rate,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
plt.title("Counts/Sec vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("X-Ray Photon Counts/Sec")
plt.legend()
plt.show()