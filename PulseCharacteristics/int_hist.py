#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
Plots histogram of the integrated intensities of pulse profiles for N pulses per profile
"""

def int_hist(N): #input is the number of pulses per profile

    #Read in intensities from file 
    intint = pd.read_csv('intdata/crabintdata_%s.txt' %N, header = None) 
    intint = list(intint[0])

    #Remove outliers that are not in the main distribution and wont affect modeling (currently set up for 15 pulses per profile)
    intint = [x for x in intint if x > 0.001] #remove values that are too small
    intint = [x for x in intint if x < 3] #remove any values that are too big
    intint = np.array(intint)

    print('The total number of profiles in the %s pulse histogram is '%N, len(intint))

    #Plot the integrated intensity distribution as a histogram
    binwidths = list(np.linspace(0, 0.12, 100)) #bins currently set up for 15 pulses per profile
    plt.hist(intint, bins=binwidths)

    plt.xlabel('Integrated intensity (counts/pulse)')
    plt.ylabel('# of Profiles')
    plt.title('Integrated intensity distribution for %s pulses/profile'%N)
    #plt.savefig('intplots/crab_%s.png'%N) #can save the plot
    plt.show()
    #plt.clf()
    return()

#Plot the intensity distribution for N pulses; must have a corresponding text file containing intensity values
N_pulses = [15]
for N in N_pulses:
    int_hist(N)


