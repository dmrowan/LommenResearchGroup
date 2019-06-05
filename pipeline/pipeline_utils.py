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
import numpy as np

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

def get_exposure_table():
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

def get_exposure(tex_out):
	df = print_nicer_segment(username='nicer_team', password='sextant')
	source = df[df['Target Name']==sourcename]

	n_obdIDs = len(source)
	good_exp = np.sum(source['Good Expo[s]'])

	heads = ['1821', '1937', '0218']
	exps = []
	for sourcename in ['PSR_B1821-24', 'PSR_B1937+21', 'PSR_J0218+4232']:
		source = df[df['Target Name']==sourcename]
		n_obdIDs = len(source)
		good_exp = np.sum(source['Good Expo[s]'])
		exps.append(good_exp)

	with open(tex_out, 'w') as f:
		f.write("\\newcommand{\\"+command_head+"_raw_exp}{"+total_raw+"}\n")
		f.write("\\newcommand{\\"+command_head+"_clean_exp}{"+total_clean+"}\n")
		f.write("\\newcommand{\\"+command_head+"_cut_exp}{"+total_cut+"}\n")

		f.write("\\newcommand{\\"+command_head+"_raw_exp}{"+total_raw+"}\n")
		f.write("\\newcommand{\\"+command_head+"_clean_exp}{"+total_clean+"}\n")
		f.write("\\newcommand{\\"+command_head+"_cut_exp}{"+total_cut+"}")


def print_nicer_segment(url = 'https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html',
                        username = None, password=None):
    """
    This prints out the segment detail table in text format

    usage: % print_nicer_segment(username = "nicer_user_name" password = "nicer password")

    outputs: prints the nicer segment table to the terminal

    :param url: location of the segment detail page
    :param username: nicer team username
    :param password: nicer team password
    :return:
    """
    from bs4 import BeautifulSoup
    import requests
    if (not username) or (not password):
        raise ValueError("must supply username and password to access the NICER obs page")
    req = requests.get(url, auth=(username, password))
    if req.status_code != 200:
        raise ValueError('Problem accessing {0} with ({1}, {2}) \nReturn code: {3}'.format(
            url, username, password, req.status_code))

    soup = BeautifulSoup(req.text, 'lxml')
    tabs = soup.find_all('table')[1]
    df = pd.read_html(str(tabs))
    return df[0]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("function", help="Options: clean, version_check, get_dates", type=str)
	parser.add_argument("--source", help="Sourcename for get exposure", type=str, default=None)
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
		get_exposure(args.source)
	else:
		print("Invalid function input")



