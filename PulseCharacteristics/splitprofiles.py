#!/usr/bin/env python
import numpy as np

"""
Function used by ampread.py; splits a single ObsID into pulse profiles containing N pulses each.
I think it's a little too easy to get confused with the naming of all the different lists here, especially as many are nested, so here are all the lists that will be used, with what each looks like.

Note a pulse=rotation of the pulsar; (pulse)phase=the PULSE_PHASE value from the event file; (pulse)profile=list of N pulses that will be plotted

pulses: list of lists of GTI intervals, where each GTI intervals list is a list of pulses, where each list of pulses contains the PULSE_PHASE values
     Ex:  pulses = [[GTI interval 1], [GTI interval 2]...]
          GTI interval 1 = [[pulse 1], [pulse 2]...]
          pulse 1 = [PULSE_PHASE value 1, PULSE_PHASE value 2...]

pulses2: list of pulses, no longer split by GTI interval 
     Ex:  pulses2 = [[pulse 1], [pulse 2]...]
          pulse 1 = [PULSE_PHASE value 1, PULSE_PHASE value 2...]

profiles: list of list of profiles, where each profile contains N pulses
     Ex:  profiles = [[pulse profile 1], [pulse profile 2]...]
          pulse profile 1 = [[pulse 1], [pulse 2] .... [pulse N]]
          pulse 1 = [PULSE_PHASE value 1, PULSE_PHASE value 2...]

phases: list of profiles, where each profile contains all the phase values for all pulses 
     Ex:  phases = [[pulse profile 1], [pulse profile 2]...]
          pulse profile 1 = [PULSE_PHASE value 1, PULSE_PHASE value 2...]

"""

def splitprofile(n_rotations, timetab, phase, times): #need N rotations per profile, table of GTI values, list of PULSE_PHASE, and list of TIME

    # Splits a full ObsID into pulses (lists of phases from one pulse), organized by GTI interval
    pulses = [] #phase values for each pulse for each GTI interval (triple nested list)
    counter = 0 #counter to keep track of index of each phase and time
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
        pulses.append(pulsestemp) #append list of pulses for that GTI interval into the list named "pulses"

    # Currently, the list "pulses" is split up by GTI interval; need to just have a list of pulses, not split up by interval
    pulses2 = [] #list of lists of phases in each pulse
    for n in pulses: 
        if len(n) != 0: #if not an empty list
            for i in n:
                pulses2.append(i) #append pulses into pulses2, no longer split by interval

    # Split into profiles that contain n_rotations pulses each
    n_profiles = int(len(pulses2)/n_rotations) #the total number of full profiles that can be made from this ObsID
    profiles = [] 
    while len(profiles) < n_profiles: #create profiles until all full profiles that can be made are filled
        profiles.append(pulses2[:n_rotations]) #append list of first n_rotations pulses into list profiles
        pulses2 = pulses2[n_rotations:] #remove first n_rotations pulses from list pulses2 so next time we can again take the first n_rotations pulses

    #flatten each profile so that phases is a list of lists of phases in each profile
    phases = []
    for i in profiles:
        phases.append(list(np.concatenate(i).flat)) #flatten profile 

    return(phases)
