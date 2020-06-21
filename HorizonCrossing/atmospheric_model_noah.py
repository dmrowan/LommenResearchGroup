#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.integrate import quad

#Time range around the horizon crossing
startTime = 390+1.92224*10**8
stopTime = 500+1.92224*10**8

startTimeIndex = 311883
stopTimeIndex = 352360

#read in the data files


tab_ni = Table.read('ni2200300102.mkf',hdu=1)
timeArray = np.array(tab_ni['TIME'])
elevArray = np.array(tab_ni['ELV'])

tab_evt = Table.read('cleanfilt.evt',hdu=1)
eventTime = np.array(tab_evt['TIME'][startTimeIndex:stopTimeIndex])

#interpolate the times.evt to go over the range of elevations.mkf
f = interpolate.interp1d(timeArray,elevArray,kind='linear')
elev_evt = f(eventTime)

#calculate altitude based on elevation angle, h(elevation)
R = 6378 #km
H = 410 #km
theta = np.arcsin(R/(R+H))
altArray = []
for val in elev_evt:
  altArray.append((R+H)*np.sin(theta+val*(np.pi/180))-R)
altArray=np.array(altArray) #h

'''
Variable Dictionary

R = Earth's radius
H = approximate altitude of ISS
h = tangential altitude of line of sight
z = radial altitude from Earth to line of sight
dtot = length of line of sight through atmosphere
'''
#Definition of Parameters
sigmaN = 3.38*.1**(-3) #cm^2/g
z0 = 80 #km
l = 8.5 #km, scalar height of atmosphere
m = 1 #nondimensional factor
p0 = .0012 * np.exp(z0/l) #g/cm^3
dtot = 2 * np.sqrt((R+H)**2-(R+altArray)**2)

for hi in range(len(altArray)):
  X = np.linspace(0,.5*dtot[hi]*10**5,100)


#Functions that build toward creating an array for percent transmittance

def rho(x,elevation):
  def z(x,elevation):
    def h(elevation):
      h_array = []
      for e in elevation:
        h_array.append((R+H)*np.sin(theta+e*(np.pi/180))-R)
      return np.array(h_array)
    z_array = []
    for h in h(elevation):
      z_array.append(np.sqrt((R+h)**2+x**2)-R)
    return np.array(z_array)
  rho_array = []
  for z in z(x,elevation):
    rho_array.append(p0*np.exp((-m*(z-z0))/l))
  return np.array(rho_array)

RHO = np.array([i[0] for i in rho(X,elev_evt)])


trans_array = []
for hind in range(len(altArray)):
  trans_array.append(np.exp(sigmaN*-2*RHO[hind]*(len(elev_evt)/(2*len(dtot))))*100)
trans_array = np.array(trans_array)

plt.plot(altArray,trans_array)
plt.show()

