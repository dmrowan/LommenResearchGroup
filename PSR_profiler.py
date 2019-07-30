import scipy
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import math
import scipy.optimize as fitter
from astropy.io import fits
import matplotlib
#change global variable = # of histogram bins
BINS = 300
#read evt file and get photon-get_phases
def getPhases(evtFile):
    f=fits.open(evtFile)
    phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
    mask = f['EVENTS'].data.field('TIME')<999999999
    phases = phases[mask]
    f.close()
    return phases
#Define all possible FITS :
def asymmetricalgaussian(x, amp1, cen1, sigma0, sigma1, amp2, cen2, sigma2, offs):
   lowind = np.where(x < cen1)
   highind = np.where(x >= cen1)
   lowy = amp1 * np.exp((-(x[lowind] - cen1)**2)/(2 * sigma0**2))
   highy= amp1 * np.exp((-(x[highind] - cen1)**2)/(2 * sigma1**2))
   return np.concatenate((lowy,highy), axis=None) + amp2 * np.exp((-(x - cen2)**2)/(2 * sigma2**2)) + offs
def asymmetricallorentzian(x, amp1, cen1, width0, width1, amp2, cen2, width2, offs):
   lowind = np.where(x < cen1)
   highind = np.where(x >= cen1)
   lowy = amp1 / (1+((x[lowind]-cen1)/width0)**2)
   highy= amp1 / (1+((x[highind]-cen1)/width1)**2)
   return np.concatenate((lowy,highy), axis=None) + amp2 / (1+((x-cen2)/width2)**2) + offs
   '''
   Define Moffat function (here for documentation )==> see model for composite function
   * alpha determines the overall shape of the profile
   * gamma is a scale factor such that FWHM = 2*R for Moffat_func(R)= 0.5 amp
   * thus gamma = FWHM/(2*((2**(1/alpha)-1)**1/2))
   '''
def Moffats(x,amp1,x01,fwhm1,alpha1,amp2,x02,fwhm2,alpha2,offs):
    hwhm1= fwhm1/2
    gamma1= hwhm1*((2**(1/alpha1)-1)**(-0.5))
    hwhm2= fwhm2/2
    gamma2= hwhm2*((2**(1/alpha2)-1)**(-0.5))
    return amp1*(1+(np.abs(x-x01)/gamma1)**2)**(-alpha1)+amp2*(1+(np.abs(x-x02)/gamma2)**2)**(-alpha2)+offs
def Gaussians(x,amp1,cen1,sigma1,amp2,cen2,sigma2,offs):
    return amp1*np.exp((-(x-cen1)**2)/(2*sigma1**2))+amp2*np.exp((-(x-cen2)**2)/(2*sigma2**2))+offs
def triple_Gaussian(x,amp1,cen1,sigma1,amp2,cen2,sigma2,amp3,cen3,sigma3,offs):
    return amp1*np.exp((-(x-cen1)**2)/(2*sigma1**2))+amp2*np.exp((-(x-cen2)**2)/(2*sigma2**2))+amp3*np.exp((-(x-cen3)**2)/(2*sigma3**2))+offs
def Lorentzians(x,amp1,cen1,fwhm1,amp2,cen2,fwhm2,offs):
    return amp1/(1+((x-cen1)/fwhm1)**2)+amp2/(1+((x-cen2)/fwhm2)**2)+offs
#Define template-writing function => kept itemplate-nomenclature for use with other scripts
def writeTemplate(opt,cov,output):
    file=[]
    dashes = '-'*25
    nprims = len(opt)//3
    i=0
    k=0
    while i< nprims:
        print(nprims)
        i+=1
        amp =[]
        amp.append(opt[0+k])
        amp.append(np.sqrt(np.diag(cov)[0+k]))
        print(amp)
        phas = []
        phas.append(opt[1+k])
        phas.append(np.sqrt(np.diag(cov)[1+k]))
        print(phas)
        fwhm=[]
        fwhm.append(opt[2+k])
        fwhm.append(np.sqrt(np.diag(cov)[2+k]))
        print(fwhm)
        k+=3
        for st,va in zip(['phas','fwhm','ampl'],[phas,fwhm,amp]):
            file += ['%s%d = %.5f +/- %.5f'%(st,i,va[0],va[1])]
    const = 'const = %.5f +/- %.5f'%(opt[len(opt)-1],np.sqrt(np.diag(cov)[len(opt)-1]))
    rstring = [dashes] + [const] + file + [dashes]
    f = open('itemplate.'+output,'w')
    f.write('#'+output+'\n')
    for s in rstring: f.write(s+'\n')
    f.close()
#Define function to compute X2-goodness of fit estimator
def X2_compute(fit,nbins):
    err= (fit*float(nbins)/len(phases))**0.5
    chi_num = (yarr-fit)**2
    chi_denom=err**2
    chi_unit=chi_num/chi_denom
    i=0
    chi=0
    while i< len(chi_unit):
        chi+=chi_unit[i]
        i+=1
    return chi
#Plot residuals
def add_Residuals(model):
    axes_resid = ax.twinx()
    resids = yarr/model-1
    err= (model*float(BINS)/len(phases))**0.5
    axes_resid.errorbar(xarr,resids,yerr=err/model,color='green',ls=' ',marker='o',alpha=0.5)
    axes_resid.axhline(0)
    axes_resid.set_ylabel('Fractional Residuals')
    axes_resid.axis([0,1,axes_resid.axis()[2],axes_resid.axis()[3]])
#main-method
if __name__== '__main__':
    desc = 'Fit a pulsar profile to one of four functions'
    parser=OptionParser(usage=" %prog [options] [filename]",description=desc)
    parser.add_option("-P","--pulsar",action='store',type='string',default='unrecognized psr',help='name of pulsar ',dest='pulsar_name')
    parser.add_option("-F","--fit",action='store',type='string',default='none',help='fitting_method : lorentzian .. alorentz ...gaussian ..agauss..tripleGauss... moffat',dest='fit')
    parser.add_option("-R","--residuals",action='store',type='string',default='true',help='compute and plot fitting residuals',dest='res')
    (options,args) = parser.parse_args()
    print(options.fit)
    #raise error
    if len(args) < 1:
        raise ValueError('Must provide an input evt/fits file!')
    phases = getPhases(args[0])
    #Plot histogram = light_curve
    fig,ax = plt.subplots(1,1,figsize=(8,5))
    yarr,bins,patches = plt.hist(phases,bins = BINS, color ='b', alpha=0.5, density = True)
    #create x-axis
    xarr = []
    for i in range(len(bins)-1):
        xarr.append((bins[i]+bins[i+1])/2)
    xarr=np.array(xarr)
    #perform Gaussiann fit
    if options.fit == 'gaussian':
        gopt,gcov=fitter.curve_fit(Gaussians,xarr,yarr,p0=[2,0.72843,0.023,1.5,0.28,0.04,0.9])
        ax.plot(xarr,Gaussians(xarr,*gopt),'r-')
        print('The value of the CHI2 test for the Gaussian fit is ',X2_compute(Gaussians(xarr,*gopt),BINS) )
        if options.res == 'true':
            add_Residuals(Gaussians(xarr,*gopt))
    #perform Moffat fit
    if options.fit == 'moffat':
        popt,pcov = fitter.curve_fit(Moffats,xarr,yarr,p0=[3.17,0.027,0.03,4.765,1.35737,0.557907,0.02,4.765,0.9],bounds=(0,10))
        ax.plot(xarr,Moffats(xarr,*popt),'r-')
        #correct in case of null convergence
        i=0
        if (popt[i] == 0) :
            popt[i]=0.0000001
            i+=1
        print('The value of the CHI2 test for the Moffat fit is ',X2_compute(Moffats(xarr,*popt),BINS) )
        if options.res == 'true':
            add_Residuals(Moffats(xarr,*popt))
        writeTemplate(popt,pcov,'moffat')
    #perform Lorentzian fitting
    if options.fit == 'lorentzian':
        lopt,lcov= fitter.curve_fit(Lorentzians,xarr,yarr,p0=[2,0.72843,0.023,1.5,0.28,0.04,0.9])
        plt.plot(xarr,Lorentzians(xarr,*lopt),'r-')
        print('The value of the CHI2 test for the Lorentzian fit is ',X2_compute(Lorentzians(xarr,*lopt),BINS))
        if options.res == 'true':
            add_Residuals(Lorentzians(xarr,*lopt))
        writeTemplate(lopt,lcov,'lorentz')
    #add triple-Gaussian
    if options.fit == 'tripleGauss'
        topt,tcov= fitter.curve_fit(triple_Gaussian,xarr,yarr,p0=[3.17,0.027,0.03,1.35737,0.557907,0.02,1.2,0.04,0.05,0.9])
        plt.plot(xarr,triple_Gaussian(xarr,*topt),'b-')
        print('The value of the CHI2 test for the triple Gaussian fit is ',X2_compute(triple_Gaussian(xarr,*topt),BINS))
        if options.res == 'true':
            add_Residuals(triple_Gaussian(xarr,*topt))
        writeTemplate(topt,tcov,'gauss')
    #fit to asymmetrical Gaussian
    if options.fit == 'agauss':
        aopt,acov = fitter.curve_fit(asymmetricalgaussian,xarr,yarr,p0=[2,0.72843,0.023,0.1,1.5,0.28,0.04,0.9])
        plt.plot(xarr,asymmetricalgaussian(xarr,*aopt),'b--')
        print('The value of the CHI2 test for the assymetrical Gaussian fit is ',X2_compute(asymmetricalgaussian(xarr,*aopt),BINS))
        if options.res == 'true':
            add_Residuals(asymmetricalgaussian(xarr,*aopt))
    #fit to asymmetrical Lorentzians
    if options.fit == 'alorentz':
        fopt,fcov = fitter.curve_fit(asymmetricallorentzian,xarr,yarr,p0=[2,0.72843,0.023,0.1,1.5,0.28,0.04,0.9])
        plt.plot(xarr,asymmetricallorentzian(xarr,*fopt),'b--')
        print('The value of the CHI2 test for the assymetrical Lorentzian fit is ',X2_compute(asymmetricallorentzian(xarr,*fopt),BINS))
        if options.res == 'true':
            add_Residuals(asymmetricallorentzian(xarr,*fopt)


    plt.show()
