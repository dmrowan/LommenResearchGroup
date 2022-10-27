#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
Plots histogram of the integrated intensities of pulse profiles for N pulses per profile
"""

def int_hist(N):

    intint = pd.read_csv('intdata/crabintdata_%s.txt' %N, header = None)
    intint = list(intint[0])
    intint = [x for x in intint if x > 0.001]
    intint = [x for x in intint if x < 3]
    intint = np.array(intint)

    print('The total number of profiles in the %s pulse histogram is '%N, len(intint))

    binwidths = list(np.linspace(0, 0.12, 100))
    plt.hist(intint, bins=binwidths)

    plt.xlabel('Integrated intensity (counts/pulse)')
    plt.ylabel('# of Profiles')
    plt.title('Integrated intensity distribution for %s pulses/profile'%N)
    plt.show()
    #plt.savefig('crab_%s.png' % timewidth)
    #plt.clf()
    return()

N_pulses = [15]
for N in N_pulses:
    int_hist(N)


