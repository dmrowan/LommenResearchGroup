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

def pickle_from_table(table, outname, hdu=1):

	assert(os.path.isfile(table))

	tab = Table.read(table, hdu=hdu)
	with open(outname, 'wb') as f:
		pickle.dump(tab, f)

def get_range(fname):
	r = [''.join(filter(str.isdigit, fname.split('--')[0])),
		 ''.join(filter(str.isdigit, fname.split('--')[1]))]
	return r

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

