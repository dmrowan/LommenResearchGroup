#!/usr/bin/env python
from __future__ import print_function, division
from astropy.table import Table
from astropy import log
import datetime
import os
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import shutil
import subprocess
import time
from bs4 import BeautifulSoup
import getpass
import pint
import requests
import niutils

obsID = '2013010104'
clobber = False
user = 'nicer_team'
passwd = 'sextant'
silent_curl = False
outdir = "./"
sourcename = 'PSR_B0531+21'

if outdir == "./":
if not os.getcwd().endswith(sourcename):
    if not outdircheck(sourcename):
        return 0

cmd = ['custom_data_download.py', sourcename, user, pwd, '--outdir', outdir, '--decryptkey', decryptkey, '--unzip', '--obsIDs', obsID]

subprocess.call(cmd)
