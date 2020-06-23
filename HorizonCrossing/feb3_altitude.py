#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import functions_mod as f
from scipy import signal

# call the function module 'f'
low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)


plt.figure(1)

plt.plot(low_en.new_alt, low_en.rate, 'r.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV')
# plt.plot(low_en.new_alt, f.SeventhOr(low_en.new_alt, f.*low_en.popt_rateAlt), 'r-')

plt.plot(mid_en.new_alt, mid_en.rate, 'g.', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
# plt.plot(mid_en.new_alt, f.SeventhOr(mid_en.new_alt, f.*mid_en.popt_rateAlt), 'g-')

plt.plot(high_en.new_alt, high_en.rate, 'b.', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')
# plt.plot(high_en.new_alt, f.SeventhOr(high_en.new_alt, f.*high_en.popt_rateAlt), 'b-')

plt.title("Counts/Sec vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("X-Ray Photon Counts/Sec")
plt.legend()



plt.figure(2)

plt.plot(low_en.new_alt, low_en.perc_trans, 'r--', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV')
# plt.plot(low_en.new_alt, f.SeventhOr(low_en.new_alt, f.*low_en.popt_perc), 'r-')
#plt.plot(low_en.new_alt, signal.savgol_filter(low_en.perc_trans, 51, 3), 'r--')

plt.plot(mid_en.new_alt, mid_en.perc_trans, 'g--', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
# plt.plot(mid_en.new_alt, f.SeventhOr(mid_en.new_alt, f.*mid_en.popt_perc), 'g-')
#plt.plot(mid_en.new_alt, signal.savgol_filter(mid_en.perc_trans, 51, 3), 'g--')

plt.plot(high_en.new_alt, high_en.perc_trans, 'b--', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')
# plt.plot(high_en.new_alt, f.SeventhOr(high_en.new_alt, f.*high_en.popt_perc), 'b-')
#plt.plot(high_en.new_alt, signal.savgol_filter(high_en.perc_trans, 51, 3), 'b--')


plt.title("Percent Transmittance vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("Percent Transmittance of X-Rays")
plt.xlim(50,250)
plt.legend()



plt.figure(3)




plt.plot(low_en.alt, low_en.trans_model, 'r--', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV')
plt.plot(mid_en.alt, mid_en.trans_model, 'g--', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
plt.plot(high_en.alt, high_en.trans_model, 'b--', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')
plt.title("Simple Atmospheric Model: Percent Transmittance vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("Percent Transmittance of X-Rays")
plt.xlim(50,250)
plt.legend()


plt.show()
