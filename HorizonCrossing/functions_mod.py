#!/usr/bin/env python
from astropy.table import Table
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
from astropy.io import ascii
# from scipy.integrate import quad

# Time range around the horizon crossing
startTime = 390+1.92224*10**8
stopTime = 500+1.92224*10**8

# indices for the event file
startTimeIndex = 311883
stopTimeIndex = 352360

# read in the data files
tab_ni = Table.read('ni2200300102.mkf', hdu=1)
timeArray = np.array(tab_ni['TIME'])
elevArray = np.array(tab_ni['ELV'])
azArray = np.array(tab_ni['RAM_ANGLE'])
enArray_low = np.array(tab_ni['FPM_XRAY_PI_0035_0200'])
enArray_mid = np.array(tab_ni['FPM_XRAY_PI_0800_1200'])

tab_evt = Table.read('cleanfilt.evt', hdu=1)
eventTime = np.array(tab_evt['TIME'][startTimeIndex:stopTimeIndex])
enArray = np.array(tab_evt['PI'][startTimeIndex:stopTimeIndex])

# read in MSIS model data
data = ascii.read("msis_model.txt")

height = np.array(data['km'])
density = np.array(data['g/cm^3'])
temp = np.array(data['K'])


# bin size and energy band cutoffs
binSize_all = 1
lowEn = [30, 70]
lowMidEn = [70, 100]
midEn = [100, 200]
highEn = [200, 600]
allEn = [0, 600]


class EnergyBands:

    def __init__(self, energy_band, bin_size):
        self.energy_band = energy_band
        self.bin_size = bin_size
        self.time, self.energies = enSplit(energy_band)
        self.alt = altSplit(energy_band)
        self.rate, self.new_alt = countRate(self.time, self.alt, bin_size)
        self.time_axis = Axis(self.rate, bin_size)
        self.perc_trans = percTrans(self.new_alt, self.rate)
        self.sigmaN = Sigma(self.energies)
        self.trans_model = Transmit(self.alt, self.sigmaN, rho0, L)


# interpolate the times.evt to go over the range of elevations.mkf
f = interpolate.interp1d(timeArray, elevArray, kind='linear')
elev_evt = f(eventTime)

g = interpolate.interp1d(timeArray, azArray, kind='linear')
az_evt = g(eventTime)

# calculate altitude based on elevation angle
R = 6378
H = 410
theta = np.arcsin(R/(R+H))
altArray = []
for indx, val in enumerate(elev_evt):
    h = ((R+H)*np.sin(theta+val*(np.pi/180)))-R
    altArray.append(h)
altArray = np.array(altArray)


# Time axis
def Axis(Rate, binSize):
    return np.arange(0, len(Rate), binSize)


# function that splits the altitudes based on energy bands
def altSplit(energy_level):
    index = np.where((enArray >= energy_level[0]) & (enArray < energy_level[1]))
    return altArray[index[0]]


# function that deduces the number of counts per bin size
def countRate(Time, alt_array, binSize):
    binCounts = []
    altitude = []
    for i in np.arange(min(Time), max(Time)+binSize, binSize):
        desind = np.where((Time >= i) & (Time < i + binSize))
        if len(alt_array[desind[0]]) != 0:
            binCounts.append(np.size(desind[0]))
            altitude.append(np.mean(alt_array[desind[0]]))
    return np.array(binCounts), np.array(altitude)


# function that makes a list of times corresponding to each energy range
def enSplit(energy_level):
    index = np.where((enArray >= energy_level[0]) & (enArray < energy_level[1]))
    return eventTime[index[0]], enArray[index[0]]


def percTrans(Alt, Rate):
    plateau = np.where(((Alt > 200) & (Alt < 250)))
    avg = np.mean(Rate[plateau[0]])
    return (Rate/avg)*100


# functions to make the atmospheric model
# altArray=h in mathematica
k = 1.38064852*10**-23
T = 500
mu = 28
mp = 1.6726219*10**-27
g = 9.8
L = (k*T)/(1000*mu*mp*g)
z0 = 135
rho0 = 0.0012*np.exp(-z0/L)


def Sigma(energy):
    c = np.float(-3)
    return (3.31*10**3)*(np.mean(energy)/100)**c


# i is the index in altArray
def Z(x, i, Alt):
    return np.sqrt(x**2+(R+Alt[i])**2)-R


def Rho(x, i, Alt, p0, l):
    return p0*np.exp(-(Z(x, i, Alt)-z0)/l)


# numerical integration
def Transmit(Alt, sigma, p0, l):
    elem = 100
    tau = []
    dist = 2*np.sqrt((R+H)**2-(R+Alt)**2)
    for hi in range(len(Alt)):
        g = 0
        x2 = (dist[hi]*10**5)/2
        X = np.linspace(0, x2, elem)
        for n in X:
            dx = x2/elem
            g += Rho(n, hi, Alt, p0, l)*dx
        tau.append(-2*sigma*g)
    tau = np.array(tau)
    trans = 100*np.exp(tau)
    return trans


# calculating fit uncertrainty based on parameter uncertainties at the point with x=X
def paramUnc(Popt, Pcov, X):
    Popt.tolist()
    fVal = SeventhOr(startTime + X, *Popt)
    frac_unc_params = []
    added_frac_unc = 0

    for paramIn in range(len(Popt)):
        Popt[paramIn] = Popt[paramIn] + np.sqrt(abs(Pcov[paramIn][paramIn]))
        fNew = SeventhOr(startTime + X, *Popt)
        frac_unc = abs(fNew-fVal)/fVal
        frac_unc_params.append((frac_unc)**2)

    for i in range(len(frac_unc_params)):
        added_frac_unc += frac_unc_params[i]

    return np.sqrt(added_frac_unc)


# best fit polynomials and scipy.curve_fit
def SeventhOr(x, a, b, c, d, e, f, g, h):
    return(a*x**7+b*x**6+c*x**5+d*x**4+e*x**3+f*x**2+g*x+h)


def curveFit(X, Y):
    popt, pcov = curve_fit(SeventhOr, X, Y)
    return popt, pcov


# Error analysis
def chiSq(Y, X, Popt, Yerr):
    return sum(((Y - SeventhOr(X, *Popt)) / Yerr)**2)
