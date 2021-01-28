#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 2021

@author: Nathaniel Ruhl and Noah Schwab

This script contains a piecewise model for transmittance vs. time for a circular 
sattelite orbit. The model uses the Beer-Lambert law and MSIS data. This script uses
"msis_419.txt", "08-10keV.dat", "cleanfilt.evt", and "ni2200300102.mkf" all from feb3.
"""
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from scipy import interpolate

#%% Define constants

R = 6371
H = 420
T = 5574 ## 5569.5 is the exact value from the convolution on feb3
omega = 2 * np.pi / T  
theta = np.arcsin(R/(R+H))

#%% Create arrays of orbital phase angles, at one second per index

eclipse_array1 = np.arange(0, theta, omega)
hc_array = np.arange(theta, np.pi/2, omega)
unattenuated_array = np.arange(np.pi/2, 3*np.pi/2, omega)
occultation_array = np.arange(3*np.pi/2, 2*np.pi-theta, omega)
eclipse_array2 = np.arange(2*np.pi-theta, 2*np.pi, omega)

full_angle_array = list(eclipse_array1) + list(hc_array) + list(unattenuated_array) + list(occultation_array) + list(eclipse_array2)

full_angle_array = np.array(full_angle_array)


#%% Confirmation of full orbit

full_angle_array = list(eclipse_array1) + list(hc_array) + list(unattenuated_array) + list(occultation_array) + list(eclipse_array2)

full_angle_array = np.array(full_angle_array)

#%% Create tangent altitude arrays for the crossings

tangentAlt_hc = (R+H)*np.sin(hc_array)-R
tangentAlt_occultation = -((R+H)*np.sin(occultation_array)+R)

#%% read in MSIS model data

data = ascii.read("msis_419.txt")

height = np.array(data['h'])
density = np.array(data['dens'])
temp = np.array(data['T'])

def msisSync_hc(Y_msis, tangentAlt):
    height[0] = tangentAlt[0]
    height[len(height)-1] = tangentAlt[len(tangentAlt)-1]
    func = interpolate.interp1d(height, Y_msis)
    return np.array(func(tangentAlt))

def msisSync_occultation(Y_msis, tangentAlt):
    height[len(height)-1] = tangentAlt[0]
    height[0] = tangentAlt[len(tangentAlt)-1]
    func = interpolate.interp1d(height, Y_msis)
    return np.array(func(tangentAlt))

msisRho_hc = msisSync_hc(density, tangentAlt_hc)
msisT_hc = msisSync_hc(temp, tangentAlt_hc)

msisRho_occultation = msisSync_occultation(density, tangentAlt_occultation)
msisT_occultation = msisSync_occultation(temp, tangentAlt_occultation)

#%% Percent Transmit function

def TransmitModel(tanAlt, density):
    tau = []
    halfDist = np.sqrt((R+H)**2-(R+tanAlt)**2)
    for hi in range(len(tanAlt)):
        g = 0
        intStepSize = 1
        upperBound1 = halfDist[hi]
        # upperBound2 = 1000  # km
        X1 = np.arange(0, upperBound1, intStepSize)
        neg3=float(-3)
        sigma = 400*10**neg3
        
        for n in X1:
            g += density[hi]*intStepSize*10**5

        tau.append(2*sigma*g)
    tau = np.array(tau)
    trans = np.exp(-tau)
    return trans

#%% Create transmittance arrays

percent_trans_eclipse1 = [0 for i in eclipse_array1]
percent_trans_eclipse2 = [0 for i in eclipse_array2]
percent_trans_unattenuated = [1 for i in unattenuated_array]
transmit_model_hc = TransmitModel(tangentAlt_hc, msisRho_hc)
transmit_model_occultation = TransmitModel(tangentAlt_occultation, msisRho_occultation)

full_transmit_model = list(percent_trans_eclipse1) + list(transmit_model_hc) + list(percent_trans_unattenuated) + list(transmit_model_occultation) + list(percent_trans_eclipse2)

full_transmit_model = np.array(full_transmit_model)

#%% Plot the the Transmit vs Time model

plt.figure(1)
plt.title('Percent Transmittance vs. Time')
plt.plot(eclipse_array1/omega, percent_trans_eclipse1, 'blue')
plt.plot(hc_array/omega, transmit_model_hc,'orange')
plt.plot(unattenuated_array/omega, percent_trans_unattenuated, 'blue')
plt.plot(occultation_array/omega, transmit_model_occultation,'orange')
plt.plot(eclipse_array2/omega, percent_trans_eclipse2, 'blue')
plt.xlabel('Time (s)')
plt.ylabel('Transmittance')

#%% Read in data from the february 3rd '08-10 keV' and compare to model

'''We are able to line up the transmission curve to the data via the elevation
angle in the mkf file. The data file "08-10keV.dat" starts at an elevation
and tangent altitude of zero and spans to 'theta' and 420 km respectively'''

data = ascii.read('08-10keV.dat')

crossing_elev = np.deg2rad(np.array(data['elev (mkf)']))
crossing_phaseAngle = crossing_elev+theta

feb3_trans_hc = np.array(data['%Transmission (evt)']/100)
feb3_trans_hc = list(feb3_trans_hc) + [1 for i in range(41)]
feb3_trans_occultation = feb3_trans_hc[::-1]

plt.figure(2)
plt.title('Percent Transmittance vs. Time')
plt.plot(eclipse_array1/omega, percent_trans_eclipse1, 'blue')
plt.plot(hc_array/omega, transmit_model_hc,'orange',label='model')
plt.plot(hc_array/omega, feb3_trans_hc,'green',label='data')
#plt.plot(unattenuated_array/omega, percent_trans_unattenuated, 'blue')
#plt.plot(occultation_array/omega, feb3_trans_occultation,'green')
#plt.plot(occultation_array/omega, transmit_model_occultation,'orange')
#plt.plot(eclipse_array2/omega, percent_trans_eclipse2, 'blue')
plt.xlabel('Time (s)')
plt.ylabel('Transmittance')
plt.legend()

#%% Read in the EVT file and calculate count rate

tabEVT = Table.read('cleanfilt.evt', hdu=1)
eventTime = np.array(tabEVT['TIME'])
enArray = np.array(tabEVT['PI'])
binSize = 1
fullTransmit = 600

def calcCountRate(timeArray):
    binCounts = []
    binTime = []
    for time_bin in np.arange(min(timeArray), max(timeArray)+binSize, binSize):
        desind = np.where((timeArray >= time_bin) & (timeArray < time_bin + binSize))
        binCounts.append(np.size(desind[0])/fullTransmit)
        binTime.append(time_bin)
    return np.array(binCounts), np.array(binTime)

binnedRate, binnedTime = calcCountRate(eventTime)

#%% Manually find an index before the third peak
#plt.plot(binnedRate,'.')
time0_indx = 11400

#%% Start EVT data at 0 sec

time0_offset = binnedTime[time0_indx]

binnedTime = binnedTime - binnedTime[time0_indx]

#%% Find the index after a full period from start point
index_fullPeriod = np.where(binnedTime > T)[0][0]

#%% Plot data against model - time on x-axis

plt.figure(3)
plt.title('Data vs Model')
plt.plot(binnedTime[time0_indx:(time0_indx+index_fullPeriod)], binnedRate[time0_indx:(time0_indx+index_fullPeriod)],'.')                                            # plot data
for n in range(1):                                                              # plot model
    plt.plot((eclipse_array1+n*2*np.pi)/omega, percent_trans_eclipse1, 'blue')
    plt.plot((hc_array+n*2*np.pi)/omega, transmit_model_hc,'orange')
    plt.plot((unattenuated_array+n*2*np.pi)/omega, percent_trans_unattenuated, 'blue')
    plt.plot((occultation_array+n*2*np.pi)/omega, transmit_model_occultation,'orange')
    plt.plot((eclipse_array2+n*2*np.pi)/omega, percent_trans_eclipse2, 'blue')
plt.xlabel('Time (s)')
plt.ylabel('Transmittance')

#%% Plot data against model - orbital phase on x-axis

plt.figure(4)
plt.title('Data vs Model')
plt.plot(binnedTime[0:index_fullPeriod]*omega, binnedRate[0:index_fullPeriod],'.')
for n in range(1):
    plt.plot((eclipse_array1+n*2*np.pi), percent_trans_eclipse1, 'orange')
    plt.plot((hc_array+n*2*np.pi), transmit_model_hc,'orange')
    plt.plot((unattenuated_array+n*2*np.pi), percent_trans_unattenuated, 'orange')
    plt.plot((occultation_array+n*2*np.pi), transmit_model_occultation,'orange')
    plt.plot((eclipse_array2+n*2*np.pi), percent_trans_eclipse2, 'orange')
plt.xlabel('Phase from Arbitrary Line (rad)')
plt.ylabel('Transmittance')

#%% Calculate time difference between horizon crossing curves @ 50% point

plt.figure(5)
plt.title('Data vs Model')
plt.plot(full_angle_array,full_transmit_model)
plt.plot(binnedTime[time0_indx:(time0_indx+index_fullPeriod)]*omega, binnedRate[time0_indx:(time0_indx+index_fullPeriod)], '.')

deltaT_50 = np.where(full_transmit_model>0.5)[0][0]-np.where(binnedRate[time0_indx:(time0_indx+index_fullPeriod)]>0.5)[0][0]

# If our bin sizes were not one second, we would need to say delta_50 = full_angle_array[deltaT_50_indx]/omega
# This is not an exact way to calculate the time in between - the first data point above 05 in the data
# is 0.6 - a better way may be to use the same 50% method but with an interpolated/exponenential trendline binnedRate, or a convolution.

t50_indx = np.where(binnedRate[time0_indx:(time0_indx+index_fullPeriod)]> 0.5)[0][0]

#%% Read in the mkf file

tab_ni = Table.read('ni2200300102.mkf', hdu=1)
timeMKF = np.array(tab_ni['TIME'])
elevMKF = np.array(tab_ni['ELV'])
ramMKF = np.array(tab_ni['RAM_ANGLE'])
position = tab_ni['POSITION']
position_mag = np.array(np.sqrt((position[:,0])**2+(position[:,1])**2+(position[:,2])**2))
positionX = np.array(position[:,0])/position_mag
positionY = np.array(position[:,1])/position_mag
positionZ = np.array(position[:,2])/position_mag
velocity = tab_ni['VELOCITY']
velocity_mag = np.array(np.sqrt((velocity[:,0])**2+(velocity[:,1])**2+(velocity[:,2])**2))

#%% Verify RA of ISS @ time0_offset - delta_T

startMKF_indx = np.where(timeMKF == round(time0_offset))[0][0]                  # index of t0

# This step assumes source is in plane of orbit
#opposite_point_indx_mkf = np.where(timeMKF == round(time0_offset-deltaT_50))[0][0]
t50_indx_mkf = np.where(timeMKF == round(time0_offset+t50_indx))[0][0]

RA_ISS = np.arctan2(positionY,positionX)
DEC_ISS = np.arcsin(positionZ)

RA_source = np.deg2rad(275)
DEC_source = np.deg2rad(-25)

print('RA of 50% point: ' + str(np.rad2deg(RA_ISS[t50_indx_mkf])))

#print('The RA of the ISS at the point opposite the source: ' + str(np.rad2deg(RA_ISS[opposite_point_indx_mkf])))  
print('The RA of the source: ' + str(np.rad2deg(RA_source)))
print('180 degree difference is expected')
print('-----')
# You will notice that the line below causes an error because the desired time is not in the mkf
#print('The DEC of the ISS at the point opposite the source: ' + str(np.rad2deg(RA_ISS[opposite_point_indx_mkf])))
print('The DEC of the source: ' + str(np.rad2deg(DEC_source)))

#%% Find RA/DEC of 50% point

# Rough comparison of the 50% point due to precession. Note that this is not accurate because
# the 50% point is not well defined for the not-horizon crossings. In order for this code to work,
# you need to read in the EVT file again without doing the time offset step.

#plt.plot(binnedRate,'.')
index_mkf_504 = np.where(timeMKF == round(binnedTime[17217]))
index_mkf_503 = np.where(timeMKF == round(binnedTime[11613]))
index_mkf_502 = np.where(timeMKF == round(binnedTime[6130]))
index_mkf_501 = np.where(timeMKF == round(binnedTime[484]))

RA_ISS = np.arctan2(positionY,positionX)

print(np.rad2deg(RA_ISS[index_mkf_501]))
print(np.rad2deg(RA_ISS[index_mkf_502]))
print(np.rad2deg(RA_ISS[index_mkf_503]))
print(np.rad2deg(RA_ISS[index_mkf_504]))


