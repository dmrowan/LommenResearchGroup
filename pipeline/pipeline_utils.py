#!/usr/bin/env python
import os
import argparse
from astropy.io import fits
from astropy.table import Table
import subprocess
import sys
import pint
from fuzzywuzzy import process
import pandas as pd
import dateutil
from tqdm import tqdm

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
	print "Python version: "
	print sys.version

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

def get_exposure():
	obsIDs = []
	clean_exp = []
	raw_exp = []
	cut_exp = []
	for d in tqdm(os.listdir("./")):
		if os.path.isdir(d):
			try:
				int_d = int(d)
				obsIDs.append(d)
				good_path=True
			except:
				continue
			
			if good_path:
				if os.path.isfile(d+"_pipe/cleanfilt.evt"):
					clean_time = Table.read(d+"_pipe/cleanfilt.evt", hdu=2).meta['EXPOSURE']
				else:
					clean_time = 0

				raw_time = Table.read(d+"/xti/event_cl/ni"+d+"_0mpu7_ufa.evt", hdu=1).meta["EXPOSURE"]

				if os.path.isfile(d+"_pipe/cleanfilt_cut.evt"):
					cut_time = Table.read(d+"_pipe/cleanfilt_cut.evt", hdu=1).meta['EXPOSURE']
				else:
					cut_time = 0

				clean_exp.append(clean_time)
				raw_exp.append(raw_time)
				cut_exp.append(cut_time)

	df = pd.DataFrame({'obsID':obsIDs, 'raw':raw_exp, 'clean':clean_exp, 'cut':cut_exp})
	df.to_csv("exposures.csv", index=False)

def calc_exposure(csv, source, to_latex=None):
    if type(evt) != str:
        raise ValueError("filename must be string")
    if not os.path.isfile(evt):
        raise FileNotFoundError

    df = pd.read_csv(calc_exposure)

    name = process.extract(source, ['PSR B1821-24', 'PSR B1937+21'], 
                           limit=1)[0][0]

    if name == 'PSR_B1821-24':
        command_head = '1821'
    else:
        command_head = '1937'

    total_raw = np.sum(df['raw'])
    total_clean = np.sum(df['clean'])
    total_cut = np.sum(df['cut'])
    
    if to_latex is not None:
        with open(to_latex, 'w') as f:
            f.write("\\newcommand{\\"+command_head+"_raw_exp}{"+total_raw+"}")
            f.write("\\newcommand{\\"+command_head+"_clean_exp}{"+total_clean+"}")
            f.write("\\newcommand{\\"+command_head+"_cut_exp}{"+total_cut+"}")

            f.write("\\newcommand{\\"+command_head+"_raw_exp}{"+total_raw+"}")
            f.write("\\newcommand{\\"+command_head+"_clean_exp}{"+total_clean+"}")
            f.write("\\newcommand{\\"+command_head+"_cut_exp}{"+total_cut+"}")
    

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("function", help="Options: clean, version_check, get_dates", type=str)
	args=parser.parse_args()

	function = process.extract(args.function, 
							   ['clean', 'version_check', 'get_dates', 'get_exposure'],
							   limit=1)[0][0]

	if function == 'clean':
		clean()
	elif function == 'version_check':
		version_check()
	elif function == 'get_dates':
		get_dates()
	elif function == 'get_exposure':
		get_exposure()
	else:
		print("Invalid function input")



