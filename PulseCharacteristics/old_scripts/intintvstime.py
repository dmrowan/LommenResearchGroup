#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math

def intints(timewidth):

    intint = pd.read_csv('crabintdataf5_%s.txt' %timewidth, header = None)
    intint = list(intint[0])
 
    print(len(intint)) 
    plt.plot(intint, ".")
    plt.show()

intints(10)
