#!/usr/bin/env python
"""
Provide a method for fitting an analytic template to photon phase
using coded pulsar-specific data for MSPs:
B1821-24 -- B1937+21 --J0218+4232
and user input for other pulsars.

Authors:
         Matthew Kerr <matthew.kerr@nrl.navy.mil>
         Paul S. Ray <paul.ray@nrl.navy.mil>
Modifications done by :
         Zaynab Ghazi <zghazi@brynmawr.edu>
"""
from __future__ import division
from __future__ import print_function
from builtins import range
from builtins import object
import numpy as np
import pylab as pl
import os
import astropy.io.fits as pyfits
from pint.templates.lcprimitives import LCGaussian,LCKernelDensity,LCEmpiricalFourier
from pint.templates.lcfitters import LCFitter
from pint.templates.lctemplate import LCTemplate
from optparse import OptionParser
import pickle

def light_curve(phases,weights=None,nbins=25,ec='blue',ls='solid',label=None,axes=None,fignum=1,nmc=100,template=None):
    if axes is None:
        pl.figure(fignum); axes = pl.gca()
    bins = np.linspace(0,1,nbins+1)
    bcs = (bins[1:]+bins[:-1])/2.0
    nph = len(phases)
    cod = axes.hist(phases,bins=bins,weights=weights,density=True,histtype='step',ec=ec)[0]
    if weights is None:
        err = (cod*float(nbins)/len(phases))**0.5
    else:
        err = np.empty([nbins,nmc])
        for i in range(nmc):
            rweights = weights[np.argsort(np.random.rand(nph))]
            err[:,i] = np.histogram(phases,bins=bins,weights=rweights,density=True)[0]
        err = np.std(err,axis=1)
    axes.errorbar(bcs,cod,yerr=err,color=ec,capsize=0,ls=' ',marker=' ')
    if (weights is not None):
        bg_level = 1-(weights**2).sum()/weights.sum()
        axes.axhline(bg_level,color=ec)
    else: bg_level = 0
    if template is not None:
        dom = np.linspace(0,1,101)
        axes.plot(dom,template(dom)*(1-bg_level)+bg_level,color='red')
    axes.axis([0,1,axes.axis()[2],axes.axis()[3]])
    axes.set_ylabel('Profile Amplitude')
    axes.set_xlabel('Pulse Phase')
    if template is not None:
        axes_resid = axes.twinx()
        model = (template(bcs)*(1-bg_level)+bg_level)
        resids = cod/model-1
        axes_resid.errorbar(bcs,resids,yerr=err/model,
                color='green',ls=' ',marker='o',alpha=0.5)
        axes_resid.axhline(0)
        axes_resid.set_ylabel('Fractional Residuals')
        axes_resid.axis([0,1,axes_resid.axis()[2],axes_resid.axis()[3]])
def get_phases(ft1file,get_weights=False,weightcol='WEIGHT',tmax=999999999):
    f = pyfits.open(ft1file)
    phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
    mask = f['events'].data.field("TIME") < tmax
    phases = phases[mask]
    if get_weights:
        weights = np.asarray(f['EVENTS'].data.field(weightcol),dtype=float)
        weights = weights[mask]
    else: weights = None
    f.close()
    return phases,weights

class BiasedFitter(object):
    def init(self):
        self.nbins = 50
        self.weights = None
        self.fignum = 1
        self.errors = False
        self.unbinned = False

    def welcome(self):
        print('\nWelcome ( ͡° ͜ʖ ͡°) . This is an updated non-interactive version of nitemplate.py the NICERsoft likelihood estimation script.')

    def __init__(self,phases,**kwargs):
        self.init()
        self.__dict__.update(**kwargs)
        self.phases = phases
        self.primitives = []
        self.norms = []
        self.dom = np.linspace(0,1,100)
        self.welcome()
        pl.close(self.fignum)
        self.fig = pl.figure(self.fignum)
        self.ax  = pl.gca()
        self.calibrate()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax)
        pl.show(block=False)
        pl.pause(0.00000000000000000000000000000000000001)
        pl.close()

    def do_fit(self,doshow=False):
        ubstr = 'unbinned' if self.unbinned else 'binned'
        template = LCTemplate(self.primitives,norms=self.norms)
        fitter   = LCFitter(template,self.phases,weights=self.weights)
        fitter.fit(estimate_errors=self.errors,use_gradient=True,unbinned=self.unbinned)
        self.fig = pl.figure(self.fignum)
        self.ax = pl.gca()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,template=template)
        if doshow:
            pl.show()
        self.fitter = fitter

    def calibrate(self):
        i=0
        identified = False
        while i<2:
            i+=1
            #Pulsar-specific execution
            if (options.pulsar_name == 'B1821') :
                identified = True
                #main_pulse estimates
                if (i == 1):
                    fwhm = 0.003
                    phase = 0
                    peak = 1.6
                #interpulse estimates:
                else:
                    fwhm = 0.004
                    phase = 0.55
                    peak = 1.4
            if (options.pulsar_name == 'B1937'):
                identified = True
                #main_puse estimates
                if (i==1):
                    fwhm = 0.025
                    phase =0.028
                    peak = 2.3
                #interpulse estimates
                else:
                    fwhm = 0.035
                    phase = 0.56
                    peak = 1.2
            if (options.pulsar_name == 'J0218'):
                identified = True
                #main_puse estimates
                if (i==1):
                    fwhm = 0.13
                    phase =0.001
                    peak = 2.3
                #interpulse estimates
                else:
                    fwhm = 0.12
                    phase = 0.49
                    peak = 1.2
            if (identified == False) :
                fwhm = float(input('Enter FWHM of pulse :'))
                phase = float(input('Enter phase of pulse:'))
                peak = float(input('Enter amplitude of pulse:'))
            #compute standard deviation
            sigma=fwhm/(8*np.log(2))**0.5
            #correct template
            if len(self.primitives)>0:
                template = LCTemplate(self.primitives,norms=self.norms)
                if peak>template(phase):
                    peak-= template(phase)
            ampl = peak*sigma*(2*np.pi)**0.5
            self.primitives.append(LCGaussian(p=[sigma,phase]))
            self.norms.append(ampl)
            norms=np.asarray(self.norms)
            if norms.sum() >1:
                norms *= 1./norms.sum()
            self.do_fit()
            template = LCTemplate(self.primitives,norms=norms)
            pl.clf()
            light_curve(self.phases,weights=self.weights,nbins=self.nbins,template=template)
            pl.draw()


    def write_template(self,outfile):
        if not hasattr(self,'fitter'):
            print('Must do fit first!'); return
        self.fitter.write_template(outfile)

    def write_profile(self,outfile,nbin,integral=True,suppress_bg=True):
        if not hasattr(self,'fitter'):
            print('Must do fit first!'); return
        self.fitter.template.write_profile(outfile,nbin,integral=integral,
                suppress_bg=suppress_bg)

#equivalent to main class of java : First code evaluated
if __name__ == '__main__':

    desc="Read an FT1 file containing PULSE_PHASE info and interactively fit a template."""
    parser=OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
    parser.add_option('--nhistbins',type='int',default=50,help="Number of bins to use in phase histogram.")
    parser.add_option('--nprofbins',type='int',default=256,help="Number of bins to use output tabular profile.")
    parser.add_option('-u','--unbinned',action='store_true',default=False,help="Perform fit with unbinned likelihood.")
    parser.add_option('-a','--align',action='store_true',default=False,help="Align template such that peak falls at phase 0.")
    parser.add_option('-w','--weights',action='store_true',default=False,help='Use weighted light curve')
    parser.add_option('-c','--weightcol',type='string',default='WEIGHT',help='Column in FITS file that holds the weight')
    parser.add_option('-p','--prof',type='string',default=None,help='Output name for products (default itemplate.*')
    parser.add_option('-m','--min_weight',type='float',default=1e-2,help='Minimum weight to include in fit.')
    parser.add_option('-T','--tmax',type='float',default=999999999,help='Maximum time to include in fit.')
    parser.add_option('-e','--errors',action='store_true',default=False,help='Compute errors on components.')
    parser.add_option('-P','--pulsar',action='store',default = 'unidentified pulsar',type='string',dest='pulsar_name')

    ## Parse arguments
    (options,args) = parser.parse_args()
    if len(args) < 1:
        raise ValueError('Must provide an input FITS file!')

    phases,weights = get_phases(args[0],get_weights=options.weights,weightcol=options.weightcol)

    if options.weights:
        phases = phases[weights > options.min_weight]
        print('%d of %d photons survive weight cut'%(len(phases),len(weights)))
        weights = weights[weights > options.min_weight]

    if False:
        pass
    else:
        #main()
        intf = BiasedFitter(phases,nbins=options.nhistbins,weights=weights,errors=options.errors,unbinned=options.unbinned)
        #pl.close()
        intf.do_fit(doshow=True)
        #intf.do_fit(doshow=False)

        if options.align:
            intf.fitter.template.align_peak()

        if options.prof is not None:
            # check that specified directory exists
            out = options.prof
            dname = os.path.dirname(out)
            if len(dname) > 0 and not os.path.exists(dname):
                raise IOError('Specified directory %s does not exist!'%(os.path.dirname(out)))
            prof = options.prof
        else:
            # risk of overriding exisiting file
            prof = 'itemplate'
        intf.write_template(prof+'.gauss')
        intf.write_profile(prof+'.prof',options.nprofbins)
        pickle.dump(intf.fitter.template,open(prof+'.pickle','wb'))

print('\n \n \n WARNING \n \n \n The fit is based on the user \'s initial estimates. \n╭( ✖_✖ )╮If the fitting failed, please re-execute the script without using the < -P > flag in order to input the right initial estimates. \n \n \nIf the fitting was successful, fetch phase/ampl/fwhm data by calling < cat itemplate.gauss> \n \nGoodbye ~~')
