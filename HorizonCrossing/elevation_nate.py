#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy import interpolate
#s = np.where(elevArray==min(elevArray))
#e = np.where(elevArray==max(elevArray))

########################################################
startTime = 390+1.92224*10**8
stopTime = 500+1.92224*10**8

startTimeIndex = 311883
stopTimeIndex = 352360

########################################################

tab_ni = Table.read('ni2200300102.mkf',hdu=1)
timeArray = np.array(tab_ni['TIME'])
elevArray = np.array(tab_ni['ELV'])
enArray_low = np.array(tab_ni['FPM_XRAY_PI_0035_0200'])
enArray_mid = np.array(tab_ni['FPM_XRAY_PI_0800_1200'])

tab_evt = Table.read('cleanfilt.evt',hdu=1)
eventTime = np.array(tab_evt['TIME'])
enArray = np.array(tab_evt['PI'][startTimeIndex:stopTimeIndex])


########################################################
#interpolation for the elevation values from mkf, applied to evt times

f = interpolate.interp1d(timeArray,elevArray,kind='linear')

elev_evt = f(eventTime[startTimeIndex:stopTimeIndex])



##########################################################
R = 6378
H = 410
theta = np.arcsin(R/(R+H))
altArray = []
for val in elev_evt:
  altArray.append((R+H)*np.sin(theta+val*(np.pi/180))-R)

def elevSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return elev_evt[index[0]]
##########################################################
#function that deduces the number of counts per bin size 
def countRate(Time,binSize):
  mintime=min(Time)
  maxtime=max(Time)
  binCounts=[]
  for i in np.arange(mintime,maxtime+binSize,binSize):
    desind=np.where((Time >= i) & (Time < i + binSize))
    binCounts.append(np.size(desind[0]))
  return np.array(binCounts)

#function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
  index=np.where((enArray>=energy_level[0])&(enArray<energy_level[1]))
  return eventTime[index[0]]

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

#splits the elevation values into three separaate arrays for each energy band
lowEnergyElev = elevSplit(lowEn)
midEnergyElev = elevSplit(midEn)
highEnergyElev = elevSplit(highEn)

print(len(lowBinRate))
print(len(lowEnergyElev))
################################################################################
'''
#indices corresponding to start and stop of crossing -- I used np.where(timeArray==starttime+Xint), np.where(timeArray==starttime+Ratemax)
print(elevArray[2265])
print(elevArray[2365])
print(altArray[2265])
print(altArray[2365])

plt.plot(altArray[2265:2365],enArray_low[2265:2365],'b.')
plt.plot(altArray[2265:2365],enArray_mid[2265:2365],'r.')
plt.title('NICER Elevation During Feb 3 Crossing')
plt.xlabel('Altitude of X-Ray Transparency (km)')
plt.ylabel('Count Rate')
plt.show()
##########################################################
'''
