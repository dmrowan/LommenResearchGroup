#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import functions_mod as f

# call the function module 'f'
low_en = f.EnergyBands(f.lowEn, f.binSize_all)
mid_en = f.EnergyBands(f.midEn, f.binSize_all)
high_en = f.EnergyBands(f.highEn, f.binSize_all)
all_en = f.EnergyBands(f.allEn, f.binSize_all)


# plot the trendlines and analyze error

fig, ax = plt.subplots()
yerrLow = np.sqrt(low_en.rate)
yerrMid = np.sqrt(mid_en.rate)
yerrHigh = np.sqrt(high_en.rate)


# Low energy
popt, pcov = f.curveFit(low_en.time_axis, low_en.rate)
plt.plot(low_en.time_axis, f.SeventhOr(low_en.time_axis, *popt), 'r')
chisqLow = sum(((low_en.rate - f.SeventhOr(low_en.time_axis, *popt)) / yerrLow)**2)
print(f'The low energy fit parameters: {popt}(highest power first)')
print(f'Low Energy chi-squared: {chisqLow}')
lowDiagonal = [np.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaLow = np.sqrt(sum(lowDiagonal))
print(f'The overall fractional sigma is {sigmaLow}')
low_frac_unc = f.paramUnc(popt, pcov, 9.57002)
print(f'The overall fractional uncertainty is {low_frac_unc}')
print('....')

# Mid energy
popt, pcov = f.curveFit(mid_en.time_axis, mid_en.rate)
plt.plot(mid_en.time_axis, f.SeventhOr(mid_en.time_axis, *popt), 'g')
chisqMid = sum(((mid_en.rate-f.SeventhOr(mid_en.time_axis, *popt))/yerrMid)**2)
print(f'The mid energy fit parameters: {popt}')
print(f'Mid Energy chi-squared: {chisqMid}')
midDiagonal = [np.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaMid = np.sqrt(sum(midDiagonal))
print(f'The overal fractional sigma is {sigmaMid}')
mid_frac_unc = f.paramUnc(popt, pcov, 4.54309)
print(f'The overall fractional uncertainty is {mid_frac_unc}')
print('....')

# High energy
popt, pcov = f.curveFit(high_en.time_axis, high_en.rate)
plt.plot(high_en.time_axis, f.SeventhOr(high_en.time_axis, *popt), 'b')
chisqHigh = sum(((high_en.rate - f.SeventhOr(high_en.time_axis, *popt)) / yerrHigh)**2)
print(f'The high enery fit parameters: {popt}')
print(f'High Energy chi-squared: {chisqHigh}')
highDiagonal = [np.sqrt(abs(pcov[i][i]))/abs(popt[i]) for i in range(len(np.diagonal(pcov)))]
sigmaHigh = np.sqrt(sum(highDiagonal))
print(f'The overall fractional sigma is {sigmaHigh}')
high_frac_unc = f.paramUnc(popt, pcov, 2.13372)
print(f'The overall fractional uncertainty is {high_frac_unc}')
print('....')

plt.title('Counts Per Second vs Time (V464 Sagittarius, Feb 3)')
plt.xlabel('Time (s)')
plt.ylabel('Counts Per Second')
plt.legend()
plt.show()


# plot the data with errorbars
# ax.errorbar(low_en.time_axis, low_en.rate, fmt='r.', yerr=yerrLow, label=f'{f.lowEn[0]/100}keV-{f.lowEn[1]/100}keV')
# ax.errorbar(mid_en.time_axis, mid_en.rate, fmt='g.', yerr=yerrMid, label=f'{f.midEn[0]/100}keV-{f.midEn[1]/100}keV')
# ax.errorbar(high_en.time_axis, high_en.rate, fmt='b.', yerr=yerrHigh, label=f'{f.highEn[0]/100}keV-{f.highEn[1]/100}keV')

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
