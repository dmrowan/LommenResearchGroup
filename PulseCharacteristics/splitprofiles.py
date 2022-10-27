#!/usr/bin/env python
import numpy as np

def splitprofile(n_rotations, timetab, phase, times):

    #split data into pulses
    pulses = [] #list of lists of phases in each pulse for each interval
    counter = 0 #counter to keep track of index of each phase/time
   
    print(len(phase))
    for i in range(len(timetab)): #does each GTI interval separately
        pulsestemp = [] #temporary list of pulses
        templist = [] #temporary list of phases in each pulse
        while timetab['START'][i] <= times[counter] <= timetab['STOP'][i] and counter+1 < len(phase)-1: #while the time is in the GTI range and the counter is still in range
            templist = [] #resets temporary list after each pulse is finished
            while phase[counter] < phase[counter+1] and counter+1 < len(phase)-1: #while same pulse (phase between 0 and 1) and counter is still in range
                templist.append(phase[counter]) #append phase into temporary list of phases for that pulse
                counter += 1 #next index
            counter += 1 #next index
            pulsestemp.append(templist) #append list of phases for the pulse into the list of pulses
        pulsestemp = pulsestemp[1:-1] #remove first and last pulse in GTI interval to make sure there is an integer number of pulses
        pulses.append(pulsestemp) #append list of pulses for that GTI interval into pulses

    pulses2 = [] #list of lists of phases in each pulse
    for n in pulses: 
        if len(n) != 0: #if not an empty list
            for i in n:
                pulses2.append(i) #append pulses into pulses2, no longer split by interval

    #split into profiles
    n_profiles = int(len(pulses2)/n_rotations) #number of full profiles that can be made
    print(n_profiles, len(pulses2))
    profiles = [] 
    while len(profiles) < n_profiles: #create profiles until all full profiles that can be made are filled
        profiles.append(pulses2[:n_rotations]) #append list of first n_rotations pulses into list profiles
        pulses2 = pulses2[n_rotations:] #remove first n_rotations pulses from list pulses2

    #flatten each profile so that phases is a list of lists of phases in each profile
    phases = []
    for i in profiles:
        phases.append(list(np.concatenate(i).flat)) #flatten profile 

    return(phases)
