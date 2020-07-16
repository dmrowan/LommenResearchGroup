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
location = np.array(tab_ni['POSITION'])
locationX = location[:, 0]
locationY = location[:, 1]
locationZ = location[:, 2]

tab_evt = Table.read('cleanfilt.evt', hdu=1)
eventTime = np.array(tab_evt['TIME'][startTimeIndex:stopTimeIndex])
enArray = np.array(tab_evt['PI'][startTimeIndex:stopTimeIndex])


def Interpolator(col_mkf):
    f = interpolate.interp1d(timeArray, col_mkf, kind='linear')
    return f(eventTime)


elev_evt = Interpolator(elevArray)
az_evt = Interpolator(azArray)
Xcoord = Interpolator(locationX)
Ycoord = Interpolator(locationY)
Zcoord = Interpolator(locationZ)

# calculate altitude based on elevation angle
R = 6378
H = 410
theta = np.arcsin(R/(R+H))
altArray = []
for indx, val in enumerate(elev_evt):
    h = ((R+H)*np.sin(theta+val*(np.pi/180)))-R
    altArray.append(np.float(h))
altArray = np.array(altArray)


# read in MSIS model data
data = ascii.read("msis_model.txt")

height = np.array(data['km'])
density = np.array(data['g/cm^3'])
temp = np.array(data['K'])


def msisSync(Y_msis):
    height[0] = altArray[0]
    height[len(height)-1] = altArray[len(altArray)-1]
    func = interpolate.interp1d(height, Y_msis)
    return np.array(func(altArray))


msisRho = msisSync(density)
msisT = msisSync(temp)


# constants
binSize_all = 1
k = 1.38064852e-23
mu = 28
mp = 1.6726219e-27
g = 9.8
# L = (k*T)/(1000*mu*mp*g)
z0 = 135
# p0 = 0.0012*np.exp(-z0/L)

# Energy Band cutoffs
# bin size and energy band cutoffs
binSize_all = 1
lowEn = [30, 70]
lowMidEn = [70, 100]
lowMidEn1 = [70, 80]
lowMidEn2 = [80, 90]
lowMidEn3 = [90, 100]
midEn = [100, 200]
highEn = [200, 600]
allEn = [0, 600]
newEnRange = [80, 100]


class EnergyBands:

    def __init__(self, energy_band, bin_size):
        self.energy_band = energy_band
        self.bin_size = bin_size
        self.time, self.energies = EnergyBands.enSplit(self)
        self.alt = EnergyBands.altSplit(self)
        self.rate, self.new_alt, self.binTime, self.binnedEnergy, self.binnedX, self.binnedY, self.binnedZ = EnergyBands.countRate(self)
        self.T_pre = EnergyBands.msisSplit(self, msisT)
        self.rho_pre = EnergyBands.msisSplit(self, msisRho)
        self.rho_msis, self.T_msis = EnergyBands.countRateSync(self)
        self.perc_trans = EnergyBands.percTrans(self)
        self.sigmafit_popt, self.sigmafit_pcov = EnergyBands.modelFit_sigma(self)
        self.trans_model = EnergyBands.Transmit2(self)
        self.modelDensity = EnergyBands.hydroStaticModel(self)

    # function that splits the altitudes based on energy bands
    def altSplit(self):
        index = np.where((enArray >= self.energy_band[0]) & (
            enArray < self.energy_band[1]))
        return altArray[index[0]]

    # function that deduces the number of counts per bin size
    def countRate(self):
        binCounts = []
        binTime = []
        altitude = []
        binnedEnergy = []
        binnedX = []
        binnedY = []
        binnedZ = []
        for i in np.arange(min(self.time), max(self.time)+self.bin_size, self.bin_size):
            desind = np.where((self.time >= i) & (self.time < i + self.bin_size))
            if len(self.alt[desind[0]]) != 0:
                binCounts.append(np.size(desind[0]))
                altitude.append(np.mean(self.alt[desind[0]]))
                binTime.append(np.mean(self.time[desind[0]]))
                binnedEnergy.append(np.mean(self.energies[desind[0]]))
                binnedX.append(np.mean(Xcoord[desind[0]]))
                binnedY.append(np.mean(Ycoord[desind[0]]))
                binnedZ.append(np.mean(Zcoord[desind[0]]))
        return np.array(binCounts), np.array(altitude), np.array(binTime), np.array(binnedEnergy), np.array(binnedX), np.array(binnedY), np.array(binnedZ)

    def msisSplit(self, msis_col):
        index = np.where((enArray >= self.energy_band[0]) & (enArray < self.energy_band[1]))
        return msis_col[index[0]]

    def countRateSync(self):
        rho = []
        temp = []
        for i in np.arange(min(self.time), max(self.time)+self.bin_size, self.bin_size):
            desind = np.where((self.time >= i) & (self.time < i + self.bin_size))
            if len(self.alt[desind[0]]) != 0.:
                rho.append(np.mean(self.rho_pre[desind[0]]))
                temp.append(np.mean(self.T_pre[desind[0]]))
        return np.array(rho), np.array(temp)

    # function that makes a list of times corresponding to each energy range
    def enSplit(self):
        index = np.where((enArray >= self.energy_band[0]) & (
            enArray < self.energy_band[1]))
        return eventTime[index[0]], enArray[index[0]]/100

    def percTrans(self):
        plateau = np.where(((self.new_alt > 200) & (self.new_alt < 250)))
        avg = np.mean(self.rate[plateau[0]])
        return (self.rate/avg)*100

    # functions to make the atmospheric model
    # altArray=h
    @property
    def atmHeight(self):
        return np.array((k*self.T_msis)/(1000*mu*mp*g))

    # Transmit for fitting sigma
    def Transmit(self, sigma):
        tau = []
        halfDist = np.sqrt((R+H)**2-(R+self.new_alt)**2)
        for hi in range(len(self.new_alt)):
            g = 0
            intStepSize = 0.05e5  # 0.05km
            upperBound1 = halfDist[hi]*10**5
            upperBound2 = 1000e5  # 1000km
            X1 = np.arange(0, upperBound1, intStepSize)
            X2 = np.arange(0, upperBound2, intStepSize)
            for n in X1:
                g += self.rho_msis[hi]*intStepSize

            for n in X2:
                g += self.rho_msis[hi]*intStepSize

            tau.append(sigma*g)
        tau = np.array(tau)
        trans = 100*np.exp(-tau)

        return np.array(trans)

    # using the equation found from fitting sigma originally (msis_edge.ipynb)
    @property
    def sigmaTrend(self):
        c = np.float(-3)
        return (135.5)*(self.binnedEnergy)**c + 87.6   # 3.31e3E + const

    def Transmit2(self):
        tau = []
        halfDist = np.sqrt((R+H)**2-(R+self.new_alt)**2)
        for hi in range(len(self.new_alt)):
            g = 0
            intStepSize = 0.05e5  # 0.05 km
            upperBound1 = halfDist[hi]*10**5
            upperBound2 = 1000e5  # 1000km
            X1 = np.arange(0, upperBound1, intStepSize)
            X2 = np.arange(0, upperBound2, intStepSize)
            for n in X1:
                g += self.rho_msis[hi]*intStepSize

            for n in X2:
                g += self.rho_msis[hi]*intStepSize

            tau.append(self.sigmaTrend[hi]*g)
        tau = np.array(tau)
        trans = 100*np.exp(-tau)

        return np.array(trans)

    def modelFit_sigma(self):
        popt, pcov = curve_fit(EnergyBands.Transmit, self, self.perc_trans)
        return popt, pcov

    def hydroStaticModel(self):
        numInt = []
        A = []
        T = []
        Rho_z = []
        f = 0
        intStepSize = 0.05  # km
        for alt in np.arange(0, len(self.altExtend), 1):
            i = alt - 1
            A.append((self.altExtend[i]**2)*(k/(g*mu*mp)))
            T.append(self.tempExtend[i])

        upperIndex = len(self.altExtend) - 1
        index_zdz = len(self.altExtend) - 1
        index_z = len(self.altExtend) - 2

        for indx in np.arange(len(self.altExtend), 0, -1):

            if (indx == len(self.altExtend)):
                f = 0
                rho0 = self.rho_msis[len(self.rho_msis)-1]
                loops = upperIndex - index_zdz
                for i in range(loops):
                    f += intStepSize * rho0 * self.altExtend[index_zdz+i]**2
                numInt.append(f)
                Rho_z.append((((numInt[0]/A[0])+T[0]*rho0))/T[0])
                index_zdz -= 1
                index_z -= 1

            elif (index_z > 0):
                f = 0
                loops = upperIndex - index_zdz
                for i in range(loops):
                    rho0 = Rho_z[i]
                    f += intStepSize * rho0 * self.altExtend[index_zdz+i]**2
                numInt.append(f)
                Rho_z.append((((numInt[len(numInt)-1]/A[index_z])+T[index_zdz]*Rho_z[len(Rho_z)-1]))/T[index_z])
                index_zdz -= 1
                index_z -= 1

            else:
                f = 0
                loops = upperIndex - index_zdz
                for i in range(loops):
                    rho0 = Rho_z[i]
                    f += intStepSize * rho0 * self.altExtend[index_zdz+i]**2
                numInt.append(f)
                Rho_z.append((((numInt[len(numInt)-1]/A[len(A)-1])+T[len(T)-1]*Rho_z[len(Rho_z)-1]))/T[len(T)-2])
                index_zdz -= 1
                index_z -= 1

        numInt.reverse()
        Rho_z.reverse()
        return np.array(Rho_z)

    @property
    def altExtend(self):
        intStepSize = 0.05
        X = np.arange(min(self.new_alt), max(self.new_alt), intStepSize)
        return X

    @property
    def tempExtend(self):
        intStepSize = 0.05
        X = np.arange(min(self.altExtend), max(self.altExtend), intStepSize)
        function = interpolate.interp1d(self.new_alt, self.T_msis)
        return function(X)

    # # calculating fit uncertrainty based on parameter uncertainties at the point with x=X -- this will need to change...
    # def paramUnc(self, Popt, Pcov, X):
    #     Popt.tolist()
    #     fVal = SeventhOr(startTime + X, *Popt)
    #     frac_unc_params = []
    #     added_frac_unc = 0
    #
    #     for paramIn in range(len(Popt)):
    #         Popt[paramIn] = Popt[paramIn] + np.sqrt(abs(Pcov[paramIn][paramIn]))
    #         fNew = SeventhOr(startTime + X, *Popt)
    #         frac_unc = abs(fNew-fVal)/fVal
    #         frac_unc_params.append((frac_unc)**2)
    #
    #     for i in range(len(frac_unc_params)):
    #         added_frac_unc += frac_unc_params[i]
    #
    #     return np.sqrt(added_frac_unc)
