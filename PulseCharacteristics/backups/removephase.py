#!/usr/bin/env python
import sys

import astropy.io.fits as pyfits
import astropy.units as u
import numpy as np
from astropy import log
import fitsio
import pint.polycos as polycos
import os.path
import os

obsid = 3020580108
filename = '%s_pipe/cleanfilt2.evt'%obsid

h = fitsio.read_header(filename, ext=1)

columns = []
for i in range(1, h['TFIELDS']+1):
    if (h["TTYPE%s" %i] != 'DEADTIME') & (h["TTYPE%s" %i] != 'MPU_A_TEMP') & (h["TTYPE%s" %i] != 'PULSE_PHASE'):
        columns.append(h["TTYPE%s" %i])

print(columns)

tab = fitsio.read(filename, columns=columns , ext=1)
tab2 = fitsio.read(filename, ext=2)

print(tab[0])

fitsio.write('%s_pipe/cleanfilt3.evt'%obsid, tab, header=h)
fitsio.write('%s_pipe/cleanfilt3.evt'%obsid, tab2)
