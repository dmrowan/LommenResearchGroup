#!/usr/bin/env python
import matplotlib.pyplot as plt
import functions_mod as f


low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)


plt.plot(low_en.alt, low_en.trans_model, 'r.', markersize='5', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (expected)')
plt.plot(low_en.new_alt, low_en.perc_trans, 'b.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (measured)')
# plt.plot(mid_en.alt, mid_en.trans_model, 'g.', markersize='5', label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
# plt.plot(high_en.alt, high_en.trans_model, 'b.', markersize='5', label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
plt.title(f'Percent Transmission vs Altitude (Scale Height ={f.L}km)')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude')
plt.grid()
plt.legend()
plt.show()

##############################################################################
# plot the data
'''
plt.figure(1)
plt.plot(low_en.new_alt,low_en.rate,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
#plt.plot(low_en.new_alt,SeventhOr(low_en.new_alt,*low_en.popt),'r-')

plt.plot(mid_en.new_alt,mid_en.rate,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
#plt.plot(mid_en.new_alt,SeventhOr(mid_en.new_alt,*mid_en.popt),'g-')

plt.plot(high_en.new_alt,high_en.rate,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
#plt.plot(high_en.new_alt,SeventhOr(high_en.new_alt,*high_en.popt),'b-')

plt.title("Counts/Sec vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("X-Ray Photon Counts/Sec")
plt.legend()
plt.show()

plt.figure(2)

plt.plot(low_en.new_alt,low_en.perc_trans,'r.',label=f'{lowEn[0]/100}keV-{lowEn[1]/100}keV')
#plt.plot(low_en.new_alt,SeventhOr(low_en.new_alt,*low_en.popt_perc),'r-')
plt.plot(low_en.new_alt,scipy.signal.savgol_filter(low_en.perc_trans,51,3),'r--')

plt.plot(mid_en.new_alt,mid_en.perc_trans,'g.',label=f'{midEn[0]/100}keV-{midEn[1]/100}keV')
#plt.plot(mid_en.new_alt,SeventhOr(mid_en.new_alt,*mid_en.popt_perc),'g-')
plt.plot(mid_en.new_alt,scipy.signal.savgol_filter(mid_en.perc_trans,51,3),'g--')

plt.plot(high_en.new_alt,high_en.perc_trans,'b.',label=f'{highEn[0]/100}keV-{highEn[1]/100}keV')
#plt.plot(high_en.new_alt,SeventhOr(high_en.new_alt,*high_en.popt_perc),'b-')
plt.plot(high_en.new_alt,scipy.signal.savgol_filter(high_en.perc_trans,51,3),'b--')


plt.title("Percent Transmittance vs. Altitude for 3 Energy Bands (FEB 3)")
plt.xlabel("Altitude (km)")
plt.ylabel("Percent Transmittance of X-Rays")
plt.legend()
plt.show()
'''
#################################################################
