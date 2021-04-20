#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

# s = np.where(elevArray==min(elevArray))
# e = np.where(elevArray==max(elevArray))

tab_ni = Table.read('ni2200300102.mkf', hdu=1)
timeArray = np.array(tab_ni['TIME'])
elevArray = np.array(tab_ni['ELV'])
enArray_low = np.array(tab_ni['FPM_XRAY_PI_0035_0200'])
enArray_mid = np.array(tab_ni['FPM_XRAY_PI_0800_1200'])
########################################################
R = 6378
H = 410
theta = np.arcsin(R / (R + H))
altArray = []
for val in elevArray:
    altArray.append((R + H) * np.sin(theta + val * (np.pi / 180)) - R)

# indices corresponding to start and stop of crossing -- I used np.where(timeArray==starttime+Xint), np.where(timeArray==starttime+Ratemax)
print(elevArray[2265])
print(elevArray[2365])
print(altArray[2265])
print(altArray[2365])

plt.plot(altArray[2265:2365], enArray_low[2265:2365],
         'b.', label='0.35 keV-2.00 keV')
# plt.plot(altArray[2265:2365],enArray_mid[2265:2365],'r.')
plt.title('Atmospheric Transparency of X-Rays During Feb 3 Crossing')
plt.xlabel('Altitude of X-Ray Transparency (km)')
plt.ylabel('Count Rate')
plt.legend()
plt.show()
