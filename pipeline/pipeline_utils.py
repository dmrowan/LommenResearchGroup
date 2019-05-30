#!/usr/bin/env python
import os
import argparse
from astropy.io import fits
import subprocess
import pint
from fuzzywuzzy import process
import pandas as pd
import dateutil

#Dom Rowan 2019

desc="""
Handful of pipeline utilities for debugging default_pipe2.7
"""

def clean():
	for d in os.listdir("./"):
		if os.path.isdir(d):
			if d.endswith("_pipe"):
				if not os.path.isfile(d+"/cleanfilt_cut.evt"):
					print(d, d[:-5])
					subprocess.call(['rm', '-rf', d])
					subprocess.call(['rm', '-rf', d[:-5]])

def version_check():
	print "Heasoft version: "
	print subprocess.check_output('echo $HEADAS', shell=True)
	print("Last nicersoft git pull: ")
	print subprocess.check_output('stat -c %y /homes/pipeline/nicersoft/.git/FETCH_HEAD', shell=True)
	print "pint version: " 
	print pint.__version__

def get_dates():
	obsIDs = []
	dates = []
	for d in os.listdir("./"):
		if os.path.isdir(d):
			try:
				int_d = int(d)
				obsID_fits = fits.open(d+"/auxil/ni"+d+".att")
				obs_date = obsID_fits[1].header['DATE-OBS']
				obsIDs.append(d)
				dates.append(dateutil.parser.parse(obs_date))
			except:
				continue

	df = pd.DataFrame({"obsID":obsIDs, "date":dates})
	df = df.sort_values(by=['date'])
	df = df.reset_index(drop=True)
	print(df)

					


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("function", help="Options: clean, version_check, get_dates", type=str)
	args=parser.parse_args()

	function = process.extract(args.function, ['clean', 'version_check', 'get_dates'],
							   limit=1)[0][0]

	if function == 'clean':
		clean()
	elif function == 'version_check':
		version_check()
	elif function == 'get_dates':
		get_dates()
	else:
		print("Invalid function input")



