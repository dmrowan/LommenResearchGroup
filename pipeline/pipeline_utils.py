#!/usr/bin/env python
import os
import argparse
from astropy.io import fits
from astropy import log
from astropy.table import Table
import matplotlib.pyplot as plt
import subprocess
import sys
import pint
from fuzzywuzzy import process
import pandas as pd
import dateutil
from tqdm import tqdm
import numpy as np
import pickle

#Dom Rowan 2019

desc="""
Handful of pipeline utilities for debugging default_pipe2.7
"""


def plotparams(ax):
	'''
	Basic plot params 

	:param ax: axes to modify

	:type ax: matplotlib axes object

	:returns: modified matplotlib axes object
	'''
	ax.minorticks_on()
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	ax.tick_params(direction='in', which='both', labelsize=15)
	ax.tick_params('both', length=8, width=1.8, which='major')
	ax.tick_params('both', length=4, width=1, which='minor')
	for axis in ['top', 'bottom', 'left', 'right']:
		ax.spines[axis].set_linewidth(1.5)
	return ax
		

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
	parser.add_argument("--tex", help="Output tex name", type=str, default=None)
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
		get_exposure(args.tex)
	else:
		print("Invalid function input")

def gti_length_hist(directory, ax=None):

	table_lengths = []
	dirnames = os.listdir(directory)
	for d in dirnames:
		if d.endswith("_pipe"):
			t = Table.read(d+"/tot.gti", hdu=1)
			print(d, len(t))
			table_lengths.append(len(t))

	
	idx = np.where(np.array(table_lengths) == max(table_lengths))[0][0]
	print(dirnames[idx])
	if ax is None:
		fig, ax = plt.subplots(1, 1, figsize=(12, 6))
		ax = plotparams(ax)
		ax.hist(table_lengths, color='xkcd:darkblue', edgecolor='black', 
				bins=30)

		plt.show()

	else:
		ax.hist(table_lengths, color='xkcd:darkblue', edgecolor='black', 
				bins=30)
		return ax

def pickle_from_table(table, outname, hdu=1):

	assert(os.path.isfile(table))

	tab = Table.read(table, hdu=hdu)
	with open(outname, 'wb') as f:
		pickle.dump(tab, f)

def mkf_hist(table, key, hdu=1, bins=50, pickle_file=None, save=None):

	if pickle is None:
		tab = Table.read(table, hdu=hdu)
	else:
		tab = pickle.load(open(pickle_file, 'rb'))
	fig, ax = plt.subplots(1, 1, figsize=(12, 6))
	ax = plotparams(ax)

	assert(key in tab.colnames)

	ax.hist(tab[key], bins=bins, edgecolor='black', color='xkcd:darkblue')
	ax.set_ylabel("N Rows", fontsize=20)
	ax.set_xlabel(key, fontsize=20)

	if save is not None:
		fig.savefig(save)
	else:
		plt.show()

def make_br_gtis(merged_mkf):
	assert(os.path.isfile(merged_mkf))

	br_ranges = ['--40', '40--60', '60--80', '80--180', '180--']
	formatted_br_ranges = []
	for br in br_ranges:
		if br.split('--')[0] == '':
			formatted_br_ranges.append('(BR_EARTH < {})'.format(br.split('--')[1]))
		elif br.split('--')[1] == '':
			formatted_br_ranges.append('(BR_EARTH > {})'.format(br.split('--')[0]))

		else:
			formatted_br_ranges.append('(BR_EARTH > {0}) && (BR_EARTH < {1})'.format(*br.split('--')))

	gti_names = [ 'gti_br_earth_{}.gti'.format(br) for br in br_ranges ]


	for i in range(len(br_ranges)):
		log.info("Running nimaketime for {0}".format(formatted_br_ranges[i]))
		overonly_string = ['FPM_OVERONLY_COUNT<1',
						   'FPM_OVERONLY_COUNT<(1.52*COR_SAX**(-0.633))']
		cor_sax_string = ['(COR_SAX.gt.(1.914*KP**0.684+0.25))']
		kp_string = ['KP.lt.5']
		sunshine_string = ['(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)'.format(60)]
		br_string = [formatted_br_ranges[i]]
		extra_expr = overonly_string+cor_sax_string+kp_string+sunshine_string+br_string

		extra_expr = "("+" && ".join("%s" %expr for expr in extra_expr) + ")"

		cmd = ['nimaketime', 'infile={}'.format(merged_mkf),
			   'outfile={}'.format(gti_names[i]), 'nicersaafilt=YES', 
			   'saafilt=NO', 'trackfilt=YES', 'ang_dist=0.015', 'elv=20', 
			   'br_earth=0', 'min_fpm=38', 'underonly_range=0-200', 
			   'expr={0}'.format(extra_expr), 
			   'outexprfile=br_expr_{}.txt'.format(br_ranges[i]), 
			   'clobber=YES']


		subprocess.call(cmd)

def get_range(fname):
	r = [''.join(filter(str.isdigit, fname.split('--')[0])),
		 ''.join(filter(str.isdigit, fname.split('--')[1]))]
	return r

#Make a function for using GTIs to get 5 EVTs and MKFs
def extract_br_events(gti_list, expr_list):

	#Doing this in chunks to catch and file errors
	for fname in gti_list:
		assert('gti_br' in fname)
		assert(fname.endswith('.gti'))
		assert(os.path.isfile(fname))

	gti_ranges = []
	for fname in gti_list:
		gti_ranges.append(get_range(fname))

	expr_ranges = []
	expr_text_list = []
	for fname in expr_list:
		expr_ranges.append(get_range(fname))
		with open(fname, 'r') as f:
			expr_text = f.read()

		expr_text = '(' + expr_text.replace('\n',' ') + ')'
		expr_text_list.append(expr_text)
				  
	assert(gti_ranges == expr_ranges)


	#Need to run niextract events and ftcopy 
	for i in range(len(gti_list)):
		evt_name = gti_list[i].replace('gti', '')[1:]+'evt'
		log.info("Creating EVT {0} from GTI {1}".format(evt_name, gti_list[i]))
		cmd = ['niextract-events', 'filename=bkgd_combined.evt',
			   'eventsout={0}'.format(evt_name), 
			   'timefile={0}'.format(gti_list[i]),
			   'gti=GTI', 'clobber=yes']

		subprocess.call(cmd)
		mkf_name = evt_name.replace('evt', 'mkf')

		log.info("Copying EXPR for BR_EARTH {0} to MKF {1}".format(expr_ranges[i], mkf_name))
		cmd = ['ftcopy', '{0}[{1}]'.format('bkgd_merged.mkf', expr_text_list[i]), mkf_name,
			   'clobber=yes', 'history=yes']

		subprocess.call(cmd)


#Need to merge all ufa files for nibackgen3c50
#call from /students/pipeline/heasoft6.26
#use the evt list to get all the same files
def merge_ufa(evt_list, verbose=False):
	with open(evt_list, 'r') as f:
		paths = f.readlines()

	paths = [ p.strip('\n') for p in paths ]

	ufa_paths = []
	for p in paths:
		for i in range(7):
			ufa_paths.append('/'.join([p.split('/')[0], p.split('/')[1].strip('_pipe'), 
						  'xti', 'event_cl', 
						  'ni'+p.split('/')[1].strip('_pipe')+'_0mpu{}_ufa.evt'.format(i)]))
				 
	with open('bkgd_all_ufa_list', 'w') as f:
		for p in ufa_paths:
			f.write(p+'\n')


	output='bkgd_combined_ufa.evt'
	cmd = ['nimpumerge', 'infiles=@bkgd_all_ufa_list', 'outfile={}'.format(output)]
	"""
	cmd = ['niextract-events', 'filename=@bkgd_all_ufa_list',
		   'eventsout={}'.format(output), 'clobber=yes', 'gti=GTI']
	if verbose:
		cmd.append('chatter=5')
	"""

#subprocess.call(cmd)
	



