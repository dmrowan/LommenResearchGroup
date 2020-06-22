#!/usr/bin/env python
from scipy.optimize import curve_fit
import functions_mod as f
import matplotlib.pyplot as plt


low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)

plt.figure(1)
plt.plot(low_en.alt, low_en.trans_model, 'k-', markersize='5', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (expected)')
plt.plot(low_en.new_alt, low_en.perc_trans, 'r.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (measured)')

plt.title(f'Percent Transmission vs Altitude (Scale Height ={f.L}km)')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(mid_en.alt, mid_en.trans_model, 'k-', markersize='5', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV (expected)')
plt.plot(mid_en.new_alt, mid_en.perc_trans, 'g.', markersize='5', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV (measured)')

plt.title(f'Percent Transmission vs Altitude (Scale Height ={f.L}km)')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

plt.figure(3)
plt.plot(high_en.alt, high_en.trans_model, 'k-', markersize='5', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV (expected)')
plt.plot(high_en.new_alt, high_en.perc_trans, 'b.', markersize='5', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV (measured)')

plt.title(f'Percent Transmission vs Altitude (Scale Height ={f.L}km)')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()



#curvefit for percent transmittance
popt, pcov = curve_fit(f.Transmit, low_en.new_alt, low_en.perc_trans)

plt.figure(4)
plt.plot(low_en.new_alt, low_en.perc_trans, 'r.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (measured)')
plt.plot(low_en.new_alt,f.Transmit(low_en.new_alt,*popt),'k--', label='Atmospheric Model fit to Data')
plt.title("Percent Transmission Fit")
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("low energy parameters")

print("Expected sigmaN: ",low_en.sigmaN)
print("Measured sigmaN: ", popt[0])
print("Difference =", abs(low_en.sigmaN - popt[0]))

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[1])
print("Difference =", abs(f.rho0 - popt[1]))

print("Expected L: ", f.L)
print("Measured L: ", popt[2])
print("Difference =", abs(f.L - popt[2]))

print('----------------------')

popt, pcov = curve_fit(f.Transmit, mid_en.new_alt, mid_en.perc_trans)

plt.figure(5)
plt.plot(mid_en.new_alt, mid_en.perc_trans, 'g.', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV (measured)')
plt.plot(mid_en.new_alt,f.Transmit(mid_en.new_alt,*popt),'k--', label='Atmospheric Model fit to Data')
plt.title("Percent Transmission Fit")
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("mid energy parameters")

print("Expected sigmaN: ",mid_en.sigmaN)
print("Measured sigmaN: ", popt[0])
print("Difference =", abs(mid_en.sigmaN - popt[0]))

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[1])
print("Difference =", abs(f.rho0 - popt[1]))

print("Expected L: ", f.L)
print("Measured L: ", popt[2])
print("Difference =", abs(f.L - popt[2]))

print('----------------------')

popt, pcov = curve_fit(f.Transmit, high_en.new_alt, high_en.perc_trans)

plt.figure(6)
plt.plot(high_en.new_alt, high_en.perc_trans, 'b.', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV (measured)')
plt.plot(high_en.new_alt,f.Transmit(high_en.new_alt,*popt),'k--', label='Atmospheric Model fit to Data')
plt.title("Percent Transmission Fit")
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude (km)')
plt.grid()
plt.legend()

print("high energy parameters")

print("Expected sigmaN: ",high_en.sigmaN)
print("Measured sigmaN: ", popt[0])
print("Difference =", abs(high_en.sigmaN - popt[0]))

print("Expected rho0: ", f.rho0)
print("Measured rho0:", popt[1])
print("Difference =", abs(f.rho0 - popt[1]))

print("Expected L: ", f.L)
print("Measured L: ", popt[2])
print("Difference =", abs(f.L - popt[2]))

plt.show()
