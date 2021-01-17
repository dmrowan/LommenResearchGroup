#!/usr/bin/env python
from scipy.optimize import curve_fit
import functions_mod as f
import matplotlib.pyplot as plt
import numpy as np


low_en = f.EnergyBands(f.lowEn, f.binSize_all)
lowMid_en = f.EnergyBands(f.lowMidEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)

plt.figure(1)
# plt.plot(low_en.alt, low_en.trans_model, 'k-', markersize='5', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (expected model)')
plt.plot(low_en.new_alt, low_en.perc_trans, 'r.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (data)')

popt, pcov = curve_fit(lambda Alt, p0, l: f.Transmit(Alt, low_en.sigmaN, p0, l), low_en.new_alt, low_en.perc_trans)

plt.plot(low_en.alt, f.Transmit(low_en.alt, low_en.sigmaN, *popt), 'k--', label='data fit to model')

plt.title('Percent Transmission vs Altitude')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("low energy parameters")

print(f"Mean energy: {np.mean(low_en.energies/100)}")

print("Expected sigmaN: ", low_en.sigmaN)
# print("Measured sigmaN: ", popt[0])

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[0])

print("Expected L: ", f.L)
print("Measured L: ", popt[1])

print('----------------------')

plt.figure(2)
# plt.plot(lowMid_en.alt, lowMid_en.trans_model, 'k-', markersize='5', label=f'{f.lowMidEn[0]/100}keV-{f.lowMidEn[1]/100}keV (expected model)')
plt.plot(lowMid_en.new_alt, lowMid_en.perc_trans, 'r.', label=f'{f.lowMidEn[0]/100}keV-{f.lowMidEn[1]/100}keV (data)')

popt, pcov = curve_fit(lambda Alt, p0, l: f.Transmit(Alt, lowMid_en.sigmaN, p0, l), lowMid_en.new_alt, lowMid_en.perc_trans)

plt.plot(lowMid_en.alt, f.Transmit(lowMid_en.alt, lowMid_en.sigmaN, *popt), 'k--', label='data fit to model')

plt.title('Percent Transmission vs Altitude')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("lowMid energy parameters")

print(f"Mean energy: {np.mean(lowMid_en.energies/100)}")

print("Expected sigmaN: ", lowMid_en.sigmaN)
# print("Measured sigmaN: ", popt[0])

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[0])

print("Expected L: ", f.L)
print("Measured L: ", popt[1])

print('----------------------')

plt.figure(3)
# plt.plot(mid_en.alt, mid_en.trans_model, 'k-', markersize='5', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV (expected model)')
plt.plot(mid_en.new_alt, mid_en.perc_trans, 'g.', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV (data)')

popt, pcov = curve_fit(lambda Alt, p0, l: f.Transmit(Alt, mid_en.sigmaN, p0, l), mid_en.new_alt, mid_en.perc_trans)

plt.plot(mid_en.alt, f.Transmit(mid_en.alt, mid_en.sigmaN, *popt), 'k--', label='data fit to model')

plt.title('Percent Transmission vs Altitude')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("mid energy parameters")

print(f"Mean energy: {np.mean(mid_en.energies/100)}")

print("Expected sigmaN: ", mid_en.sigmaN)
# print("Measured sigmaN: ", popt[0])

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[0])

print("Expected L: ", f.L)
print("Measured L: ", popt[1])

print('----------------------')

plt.figure(4)
# plt.plot(high_en.alt, high_en.trans_model, 'k-', markersize='5', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV (expected model)')
plt.plot(high_en.new_alt, high_en.perc_trans, 'b.', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV (data)')

popt, pcov = curve_fit(lambda Alt, p0, l: f.Transmit(Alt, high_en.sigmaN, p0, l), high_en.new_alt, high_en.perc_trans)

plt.plot(high_en.alt, f.Transmit(high_en.alt, high_en.sigmaN, *popt), 'k--', label='data fit to model')

plt.title('Percent Transmission vs Altitude')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("high energy parameters")

print(f"Mean energy: {np.mean(high_en.energies/100)}")

print("Expected sigmaN: ", high_en.sigmaN)
# print("Measured sigmaN: ", popt[0])

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[0])

print("Expected L: ", f.L)
print("Measured L: ", popt[1])

plt.show()
