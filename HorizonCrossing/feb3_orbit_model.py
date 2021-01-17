#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 2021
Last updated: Jan 16 2021

@author: nathanielruhl

This script contains a piecewise model for transmittance vs. time for a circular 
sattelite orbit. It uses the Beer-Lambert law and MSIS data.
"""
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from scipy import interpolate

#%% Define constants

R = 6371;
H = 420;
T = 5574;
omega = 2 * np.pi / T;  
theta = np.arcsin(R/(R+H))

#%% Create arrays of orbital phase angles, at one second per index

eclipse_array1 = np.arange(0, theta+omega, omega)
hc_array = np.arange(theta, np.pi/2+omega, omega)
unattenuated_array = np.arange(np.pi/2, 3*np.pi/2+omega, omega)
occultation_array = np.arange(3*np.pi/2, 2*np.pi-theta+omega, omega)
eclipse_array2 = np.arange(2*np.pi-theta, 2*np.pi+omega, omega)


#%% Confirmation of full orbit

full_angle_array = list(eclipse_array1) + list(hc_array) + list(unattenuated_array) + list(occultation_array) + list(eclipse_array2)
'''
plt.plot(-np.cos(full_angle_array), np.sin(full_angle_array), ',')
plt.plot(-np.cos(hc_array), np.sin(hc_array), 'o')
plt.plot(-np.cos(occultation_array), np.sin(occultation_array), 'o')
'''

#%% Create tangent altitude arrays

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

#%% Read in data from the february 3rd and compare to model

data = ascii.read('08-10keV.dat')

crossing_elev = np.deg2rad(np.array(data['elev (mkf)']))
crossing_phaseAngle = crossing_elev+theta

feb3_trans_hc = np.array(data['%Transmission (evt)']/100)
feb3_trans_hc = list(feb3_trans_hc) + [1 for i in range(42)]
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

#%% Plot data against model
binnedTime = binnedTime - binnedTime[0]

plt.figure(3)
plt.title('Data vs Model')
plt.plot(binnedTime, binnedRate,'.')
plt.plot(eclipse_array1/omega, percent_trans_eclipse1, 'blue')
plt.plot(hc_array/omega, transmit_model_hc,'orange')
plt.plot(unattenuated_array/omega, percent_trans_unattenuated, 'blue')
plt.plot(occultation_array/omega, transmit_model_occultation,'orange')
plt.plot(eclipse_array2/omega, percent_trans_eclipse2, 'blue')
plt.xlabel('Time (s)')
plt.ylabel('Transmittance')



        



