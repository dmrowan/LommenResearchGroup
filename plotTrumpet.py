#!/usr/bin/env python
import genspectra
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from spectraplots import plotparams
from tqdm import tqdm 
from astropy import log
import xspec

def plotTrumpetSpecs(evtFiles):
	mainVals1937 = [0.01,0.07]
	mainVals1821 = [0.98,1.05]
	interVals1821 = [0.51,0.59]
	interVals1937 = [0.54, 0.59]
	#the lower and upper are the two different phase ranges for the main and interpulse
	#"Back" refers to the background region/ off pulse region used as background
	#We use the same background region for each pulsar
	lowBack = [.2,.4]
	upBack = [.7,.9]
	   
    
   	#for num in range(evtplotTrupetFiles): 
	#using old code in genspectra to produce a pha file for a given evt file
	genspectra.gen_spectra('evtFiles/newfile1821_4.evt', mainVals1821, interVals1821, 0, 1200, 1000, save_pha = "ourOnPeakboth.pha",run_xspec = False)
	genspectra.gen_spectra('evtFiles/newfile1821_4.evt', lowBack, upBack, 0, 1200, 1000, save_pha = "ourOffPeak.pha",run_xspec = False)

	s1 = xspec.Spectrum("ourOnPeakboth.pha")
	s1.response = "/packages/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf"
	s1.response.arf = "/packages/caldb/data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf"
	s1.background = "ourOffPeak.pha"
	s1.ignore("**-0.8,10.0-**")
		#this is the command "abund wilm"
	xspec.Xset.abund = "wilm"
		#m1 holds the model that was implemented on the data
	m1 = xspec.Model("tbabs*pow")
		#fitting and plotting the data/creating values for the data of the plot
	xspec.Fit.perform()
	print("error stuff")
	xspec.Fit.error("1-3")
	print("hello")
	xspec.Plot.device ="/xs"
	xspec.Plot("ufspec delchi")

if __name__ == '__main__':
	plotTrumpetSpecs('../PSR_B1927+21_combined.evt')


