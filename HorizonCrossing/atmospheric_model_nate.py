#!/usr/bin/env python
import functions_mod as f
import matplotlib.pyplot as plt


low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)


plt.plot(low_en.alt, low_en.trans_model, 'r.', markersize='5', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (expected)')
plt.plot(low_en.new_alt, low_en.perc_trans, 'b.', label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV (measured)')
# plt.plot(mid_en.alt, mid_en.trans_model, 'g.', markersize='5', label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
# plt.plot(f.high_en.alt, high_en.trans_model, 'b.', markersize='5', label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')
plt.title(f'Percent Transmission vs Altitude (Scale Height ={f.L}km)')
plt.ylabel('Percent Transmission (%)')
plt.xlabel('Tangential Altitude')
plt.grid()
plt.legend()
plt.show()
