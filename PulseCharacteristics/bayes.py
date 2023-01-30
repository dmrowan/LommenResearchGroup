#!/usr/bin/env python
import numpy as np
import pandas as pd
from models import *

"""
Uses Bayes’ theorem on the chi squared values for the models or the log likelihood for the models to compare the models. Highest output value indicates the best model

Note: uncertain if working correctly, if Bayes’ can be done using log likelihoods, or if it can (or needs to) be done using chi squared
"""

def chisquare(n_pulses, model1, p1, n1, s): #same as in chisq.py, but set up to only do one model
    intint = pd.read_csv('intdata/crabintdata_%s.txt' %n_pulses, header = None)
    intint = list(intint[0])
    intint = [x for x in intint if x > 0]
    intint = [x for x in intint if x < 0.2]
    intint = np.array(intint)
    ints = np.array_split(intint, s)

    pvals = []
    for section in ints:
        x1, y1, y_r1, y_err1, y_f = model(n_pulses, model1, p1, plot=False, modelplot=False, split=[section, s])

        c1 = chisq(y1, y_r1)
        c2 = 1
        n2 = 1

        if c1 <= c2:
            F = (c1/n1)/(c2/n2)
        else:
            F = (c2/n2)/(c1/n1)

        dfn = n1
        dfd = n2
        p_val = 1 - stats.f.cdf(F, dfn, dfd)
        pvals.append(p_val)

    return(pvals)

def likelihood(n_pulses, model1, p1, s): #reads in intensity data, uses the function model to fit model to data, calculates log likelihood
    
    #read in intensity data, remove outliers, split data into s sections
    intint = pd.read_csv('intdata/crabintdata_%s.txt' %n_pulses, header = None)
    intint = list(intint[0])
    intint = [x for x in intint if x > 0]
    intint = [x for x in intint if x < 0.2]
    intint = np.array(intint)
    ints = np.array_split(intint, s)

    #use function model to fit the model to the data (intensity histogram), then calculate log likelihood
    lhs = []
    for section in ints:
        x1, y1, y_r1, y_err1, y_f = model(n_pulses, model1, p1, plot=False, modelplot=False, split=[section, s])

        #normalize
        total =  np.trapz(y1)
        model_p = y_r1/total

        #calculate log likelihood
        lh = 0
        for i in range(len(y_f)):
            lh += y_f[i]*np.log10(model_p[i])
        
        lhs.append(lh)

    return(lhs)

#uses Bayes' to compare all three models using either (log) likelihood or chi squared
#we don't think we can do it with LOG likelihood (need just likelihood but we can't convert it), or chi squared
def bayes(N, method): #N is number of trials (or "sections" in functions chisquare and likelihood), method is either chisquare or likelihood
    print("Bayes' theorem for %s trials"%N)

    #calculate probability of each model
    if method == chisquare:
        pm1 = chisquare(15, 'log normal', [1,-1.3,0.8], 2, N)
        pm2 = chisquare(15, 'gaussian', [1,0.1,0.01], 1, N)
        pm3 = chisquare(15, 'power law', [1, 5.25], 1, N)

    if method == likelihood:
        pm1 = likelihood(15, 'log normal', [1,-1.3, 0.8], N)
        pm2 = likelihood(15, 'gaussian', [1, 0.1, 0.001], N)
        pm3 = likelihood(15, 'power law', [1, 5.25], N)

    #set priors as having equal probability
    priorm1 = 0.5
    priorm2 = 0.5
    
    #use Bayes' to compare model 1 and model 2
    for i in range(len(pm1)):
        postm1 = pm1[i]*priorm1/(pm1[i]*priorm1 + pm2[i]*priorm2)
        postm2 = pm2[i]*priorm2/(pm1[i]*priorm1 + pm2[i]*priorm2)

        if (postm1 + postm2 > 1.00001) or (postm1 + postm2 < 0.99999):
            print(postm1+postm2)
            raise ValueError

        priorm1 = postm1
        priorm2 = postm2

    print("M1 vs M2", postm1, postm2)

    priorm2 = 0.5
    priorm3 = 0.5

    #use Bayes' to compare model 2 and model 3
    for i in range(len(pm2)):
        postm2 = pm2[i]*priorm2/(pm2[i]*priorm2 + pm3[i]*priorm3)
        postm3 = pm3[i]*priorm3/(pm2[i]*priorm2 + pm3[i]*priorm3)

        if (postm2 + postm3 > 1.00001) or (postm2 + postm3 < 0.99999):
            print(postm2+postm3)
            raise ValueError

        priorm2 = postm2
        priorm3 = postm3

    print("M2 vs M3", postm2, postm3)

N = [1,2,4,8,16, 32]
for n in N:
    bayes(n, likelihood)
