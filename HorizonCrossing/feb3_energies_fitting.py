#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import functions_mod as f

# call the function module 'f'
low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)

# Poisson error bins
yerrLow = np.sqrt(low_en.rate)
yerrMid = np.sqrt(mid_en.rate)
yerrHigh = np.sqrt(high_en.rate)

# plot the data with errorbars
fig, ax = plt.subplots()
ax.errorbar(low_en.time_axis, low_en.rate, fmt='r.', yerr=yerrLow, label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV')
ax.errorbar(mid_en.time_axis, mid_en.rate, fmt='g.', yerr=yerrMid, label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
ax.errorbar(high_en.time_axis, high_en.rate, fmt='b.', yerr=yerrHigh, label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')


# Low energy - trendline
plt.plot(low_en.time_axis, f.SeventhOr(low_en.time_axis, *low_en.popt_rateTime), 'r')
chisqLow = f.chiSq(low_en.rate, low_en.time_axis, low_en.popt_rateTime, yerrLow)
print(f'Low Energy chi-squared: {chisqLow}')
low_frac_unc = f.paramUnc(low_en.popt_rateTime, low_en.pcov_rateTime, 9.57002)
print(f'The overall fractional uncertainty is {low_frac_unc}')
print('....')

# Mid energy - trendline
plt.plot(mid_en.time_axis, f.SeventhOr(mid_en.time_axis, *mid_en.popt_rateTime), 'g')
chisqMid = f.chiSq(mid_en.rate, mid_en.time_axis, mid_en.popt_rateTime, yerrMid)
print(f'Low Energy chi-squared: {chisqMid}')
mid_frac_unc = f.paramUnc(mid_en.popt_rateTime, mid_en.pcov_rateTime, 4.54309)
print(f'The overall fractional uncertainty is {mid_frac_unc}')
print('....')

# High energy - trendline
plt.plot(high_en.time_axis, f.SeventhOr(high_en.time_axis, *high_en.popt_rateTime), 'r')
chisqHigh = f.chiSq(high_en.rate, high_en.time_axis, high_en.popt_rateTime, yerrHigh)
print(f'Low Energy chi-squared: {chisqHigh}')
high_frac_unc = f.paramUnc(high_en.popt_rateTime, high_en.pcov_rateTime, 2.13372)
print(f'The overall fractional uncertainty is {high_frac_unc}')
print('....')

plt.title('Counts Per Second vs Time (V464 Sagittarius, Feb 3)')
plt.xlabel('Time (s)')
plt.ylabel('Counts Per Second')
plt.legend()
plt.show()


'''
def poptArray(i,sigma):
    return [popt[i]+sigma[i],popt[i]-sigma[i]]


t0=poptArray(0,sigmaLow)
t1=poptArray(1,sigmaLow)
t2=poptArray(2,sigmaLow)
t3=poptArray(3,sigmaLow)
t4=poptArray(4,sigmaLow)
t5=poptArray(5,sigmaLow)
t6=poptArray(6,sigmaLow)
t7=poptArray(7,sigmaLow)

combos=list(itertools.product(t0,t1,t2,t3,t4,t5,t6,t7))

values=[]
for i in range(len(combos)):
    a, b, c, d, e, f, g, h = combos[i]
    values.append(f.SeventhOr(low_en.time_axis, a, b, c, d, e, f, g, h))

values = np.array(values)
fitError = np.std(values, axis = 0)
'''
