#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import scipy.integrate as integrate

import niutils
import xspeclog

desc="""
Integrate a XSPEC modeled spectrum
"""

def trapezoid_sum(xd, i):
    dA = ((xd.data['energy'][i+1] - xd.data['energy'][i]) * 
           (xd.data['model'][i]+xd.data['model'][i+1])/2)
    return dA

def numerical_integrate(fname, lower_energy, upper_energy, plot=True):

    #Load data in with xspecdata object
    xd = xspeclog.xspecdata(fname)

    if plot:
        #Plotting Routine
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax = niutils.plotparams(ax)
        ax.scatter(xd.data['energy'], xd.data['model'])
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_ylim(top = 1e-3, bottom=1e-8)
        ax.set_xlabel('Energy keV', fontsize=20)
        ax.set_ylabel(r'Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$', fontsize=20)

        #ax.set_xlim(left=4.5, right=8.5)
        ax.set_xlim(left=.2, right=10)
        ax.set_ylim(bottom=10e-6, top=.0001)

    idx_lower = np.where(xd.data['energy'] == min(xd.data['energy'], 
                         key=lambda x: abs(x-lower_energy)) )[0][0]

    idx_upper = np.where(xd.data['energy'] == min(xd.data['energy'], 
                         key=lambda x:abs(x-upper_energy)) )[0][0]

    #First calculate the total integration over the energy range
    total_integration = 0
    for i in range(idx_lower, idx_upper):
        #Use trapezoidal riemann sum
        total_integration += trapezoid_sum(xd, i)

    #We want to find the energy where half of the integration occurs
    ratios = []
    idx_list = np.arange(idx_lower+1, idx_upper)
    #Iterate by increasing the lower energy bound
    for i in idx_list:
        left_sum = 0
        for left_i in range(idx_lower, i):
            left_sum += trapezoid_sum(xd, left_i)
        #At each index, compute the ratio between the calculated sum and the total integration/2
        ratios.append(left_sum/(total_integration/2))

    #The middle will be where the ratio is closest to 1
    ratio_idx_middle = np.where(np.array(ratios) == min(ratios, key=lambda x: abs(x-1)))[0][0]
    #We want the idx corresponding to the original dataframe
    idx_middle = idx_list[ratio_idx_middle]

    print(f"The center of the integrated energy range is {xd.data['energy'][idx_middle]} keV")

    if plot:
        ax.axvline(xd.data['energy'][idx_lower], color='gray', ls=':')
        ax.axvline(xd.data['energy'][idx_upper], color='gray', ls=':')
        ax.axvline(xd.data['energy'][idx_middle], color='xkcd:azure', lw=4)
        
        plt.show()

    return xd.data['energy'][idx_middle]

def polynomial_integrate(fname, lower_energy, upper_energy, plot=True, degree=8):

    #Load data in with xspecdata object
    xd = xspeclog.xspecdata(fname)

    if plot:
        #Plotting Routine
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax = niutils.plotparams(ax)
        ax.scatter(xd.data['energy'], xd.data['model'])
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_ylim(top = 1e-3, bottom=1e-8)
        ax.set_xlabel('Energy keV', fontsize=20)
        ax.set_ylabel(r'Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$', fontsize=20)

        ax.set_xlim(left=.2, right=10)
        ax.set_ylim(bottom=10e-6, top=.0001)

    coefs = poly.polyfit(xd.data['energy'], xd.data['model'], degree)
    ffit = poly.polyval(xd.data['energy'], coefs)

    def mypoly(x):
        return poly.polyval(x, coefs)

    full_integration = integrate.quad(mypoly, lower_energy, upper_energy)[0]
    ratios = []
    nspace = 100000

    ratios = [ integrate.quad(mypoly, lower_energy, upper)[0]/(full_integration/2)
               for upper in np.linspace(lower_energy, upper_energy, nspace) ]

    idx = np.where(np.array(ratios) == min(ratios, key=lambda x: abs(x-1)))[0][0]
    center_energy = np.linspace(lower_energy, upper_energy, nspace)[idx]
    
    print(f"The center of the integrated energy range is {center_energy} kev")

    if plot:
        ax.plot(xd.data['energy'], ffit, color='xkcd:violet')
        ax.axvline(center_energy, color='red')
        plt.show()


    return center_energy

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fname", help='data file', type=str)
    parser.add_argument("emin", help='lower energy in keV', type=float)
    parser.add_argument("emax", help='upper energy in keV', type=float)
    parser.add_argument("--degree", help='degree for polyfit', type=int, 
                        default=8)

    args = parser.parse_args()


    numerical_integrate(args.fname, args.emin, args.emax, plot=False)
    polynomial_integrate(args.fname, args.emin, args.emax, degree=args.degree, plot=True)





