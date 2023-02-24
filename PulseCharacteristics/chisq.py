#!/usr/bin/env python
import numpy as np
import pandas as pd
from models import *

"""
Calculates the chi squared for two models given input parameters and compares them using the p-value.

Note: uncertain if works correctly
"""

def chisquare(n_pulses, model1, model2, p1, p2, s): #use F-test p-values to compare two models given parameters for each
    #n_pulses is N pulses per profile, model1 and model 2 are the two model names, p1 and p2 are the list of parameters for model, s is the number of sections the intensity data is split into
    
    #read in intensity data, remove outliers, split into s sections 
    intint = pd.read_csv('intdata/crabintdata_%s.txt' %n_pulses, header = None)
    intint = list(intint[0])
    intint = [x for x in intint if x > 0]
    intint = [x for x in intint if x < 0.2]
    intint = np.array(intint)
    ints = np.array_split(intint, s)

    #for each section, use function model to fit the model to the data (intensity histogram) 
    #everything with '1' is for model 1, everything with '2' is model 2
    pvals = []
    for section in ints: 
        #get x and y values of data and fitted model for model 1 and 2; no plotting
        x1, y1, y_r1, y_err1, y_f = model(n_pulses, model1, p1, plot=False, modelplot=False, split=[section,s])
        x2, y2, y_r2, y_err2, y_f = model(n_pulses, model2, p2, plot=False, modelplot=False, split=[section,s])

        #use function chisq to calculate chi squared for each model
        c1 = chisq(y1, y_r1)
        c2 = chisq(y2, y_r2)
        print(model1, c1, model2, c2)

        #calculate degrees of freedom on each model 
        n1 = len(y1)-n1
        n2 = len(y2)-n2
        
        #using the model with the lower chi squared as the numerator, calculate F
        if c1 <= c2:
            F = (c1/n1)/(c2/n2)
            dfn = n1
            dfd = n2
        else:
            F = (c2/n2)/(c1/n1)
            dfn = n2
            dfd = n1

        #from F and degrees of freedom, ge the F-test p-value
        p_val = 1 - stats.f.cdf(F, dfn, dfd)
        pvals.append(p_val)

    return(pvals)

#Compare all three models
m1m2 = chisquare(15, 'log normal', 'gaussian', [1,-1.3,0.8], [1, 0.1, 0.01], 1)
m3m2 = chisquare(15, 'power law', 'gaussian', [1, 5.25], [1,0.1,0.01], 1)
m3m1 = chisquare(15, 'log normal', 'power law', [1, -1.3, 0.8], [1, 5.25], 1)

print(m1m2, m3m2, m3m1)

