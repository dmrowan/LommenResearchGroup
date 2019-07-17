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
	xValues = []
	yValues = []
	xErrors = []
	yErrors = []
	folded = []    
    
    #for num in range(evtplotTrupetFiles): 
		#using old code in genspectra to produce a pha file for a given evt file
	genspectra.gen_spectra(evtFiles, 0.0, 0.06, 0, 1200, 1000, save_pha = "ourOnPeak.pha",run_xspec = False)
	genspectra.gen_spectra(evtFiles, 0.15, 0.45, 0, 1200, 1000, save_iha = "ourOffPeak.pha",run_xspec = False)

				    ## Going through standard xspec proceedure
	s1 = xspec.Spectrum("ourOnPeak.pha")
	s1.background = "ourOffPeak.pha"
	s1.response = "/packages/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf"
	s1.response.arf = "/packages/caldb/data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf"
	s1.ignore("**-0.8,10.0-**")
		#this is the command "abund wilm"
	xspec.Xset.abund = "wilm"
		#m1 holds the model that was implemented on the data
	m1 = xspec.Model("tbabs*pow")
		#fitting and plotting the data/creating values for the data of the plot
	xspec.Fit.perform()
	xspec.Plot.device ="/xs"
	xspec.Plot("uf data")
	xValues.append(xspec.Plot.x())
	yValues.append(xspec.Plot.y())
	folded.append(xspec.Plot.model())
	xErrors.append(xspec.Plot.xErr())
	yErrors.append(xspec.Plot.yErr())

   	##Making a figure 
	fig = plt.figure(figsize = (10,11))
	plt.subplots_adjust(top = .98, right = .98, hspace = .15, left =.15)
	outer = gridspec.GridSpec(2, 1, height_ratios = [1,1])

	inner_p = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec = outer[0], hspace = 0, height_ratios = [3,1])
	inner_i = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],hspace=0, height_ratios=[3,1])

	axp1 = plt.Subplot(fig, inner_p[1])
	axp0 = plt.Subplot(fig, inner_p[0], sharex=axp1)
	axi1 = plt.Subplot(fig, inner_i[1])
	axi0 = plt.Subplot(fig, inner_i[0], sharex=axi1)
	 #A sourcename stuff goes here, may need to put it back later


	 #Labels and colors
	errorbarparams = dict(ls=' ', color='#28145b')
	labels=["Primary Pulse", "Interpulse"]
	 #creating errorbars
	for i in range(len(yValues)):
		axp0.errorbar(xValues[i], yValues[i], yerr= yErrors[i], xerr= xErrors[i], **errorbarparams, marker = 'o', label = 'Data', zorder=2)
		 #not quite sure how to add the plot params
		axp0 = plotparams(ax)
		 # text and scaling added to the plot for 1937
		axp0.text(.98, .35, 'Primary Pulse', transform= axp0.transAxes, fontsize=20, ha='right', va='top') 
		axp0.set_yscale('log')
		axp0.set_xscale('log')
		axp0.set_xlim(right=10)
		axp0.plot(xValues[i], yValues[i],ls='-', lw=3, color='#0da0ff', zorder=1, label="Powerlaw*tbabs")
		#Adds size then adds the figlcure to the plot
		axp0.legend(loc=(.68, .05), fontsize=15, edgecolor='black') #for 1937
	fig.add_subplot(axp0)

	#Plotting residuals/ErrorBars
	residuals = [np.subtract(yValues[i], folded[i])
		          for i in range(len(xVals)) ]
	for j in range(len(xValues)):
		axp1.errorbar(xValues[i], residuals[i],yerr=yErrs, ls=' ', marker='.',color='#0da0ff')
		axp1.axhline(0, ls=':', color='0.8', lw=4)
		axp1 = plotparams(axp1)
		axp1.set_xscale('log')
		axp1.set_xlim(right=10)
		
	axp1.set_ylim(bottom=-.00003, top=.00003)
	fig.add_subplot(axp1)
	#Here we are setting the tick lables 
	plt.setp(axi0.get_xticklabels(), visible=False)
	plt.setp(axp0.get_xticklabels(), visible=False)
	#We are adding text to the fig
	fig.text(.03, .55, "Normalized Flux", ha='center', va='center', rotation='vertical', fontsize=30)
	axi1.set_xlabel("Energy (keV)", fontsize=30)
	#fig.savefig("plotunfolded.jpeg", dpi=300)

	#plt.plot(xVals, yVals, 'ro', xVals, folded)
	#plt.xlabel('channels')
	#plt.ylabel('counts')
	plt.show()
	#plt.savefig('myplot')
