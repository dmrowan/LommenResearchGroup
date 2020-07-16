#!/usr/bin/env python

from astropy.table import Table
from astropy import log
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math

fnames = pd.read_csv('crabfilenames.txt', header = None)
fnames = list(fnames[0])

print(fnames)
