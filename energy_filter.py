#!/usr/bin/env python
from __future__ import print_function, division
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
from astropy.table import Table
import sys
from nicer.values import *
import argparse

desc = "This will filter an already existing event file by the energy range specified"

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("--infile",  help = "The evt file you'd like to filter", default = none)
parser.add_argument("--min",  help = "Minimum energy range to be considered", type = float, default = .25)
parser.add_argument("--max",  help = "Maximum energy range to be considered", type = float, default = 12.)
parser.add_argument("--outfile",  help = "The name of the output file", type = int, default =200)

args = parser.parse_args()


evtdata = Table.read(args.infile,hdu=1)
print (evtdata.info)
pi = evtdata['PI']
print (np.min(pi), np.max(pi))
# So for example if you want to keep all the TOAs between 5 and 10kEv
# that's between 500 and 1000 PI
lowend= args.min
hiend= args.max
mask = (evtdata['PI'] > int(lowend)) & (evtdata['PI'] < int(hiend))
filtered_table = evtdata[mask]
filtered_table.write(args.outfile, format='fits')
