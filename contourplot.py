#!/usr/bin/env python
import argparse
import collections
from fuzzywuzzy import process
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import scipy

from spectraplots import plotparams
from xspeclog import logfile

rc('text', usetex=True)

#Dom Rowan 2019
desc="""

Parse xspec generated contour plots and generate in Python

Example command line calls:

   $ contourplot.py -c leading_contour.txt trailing_contour.txt -l 
     "Leading Edge" "Trailing Edge" -o ContourEdges.pdf -s "PSR B1821-24"

   $ contourplot.py -c nicer_contour.txt nustar_contour.txt xmm_contour.txt -l
     "NICER" "NuSTAR" "XMM-Newton" -o ContourPlot.pdf

Example function calls:

    >>> plotcontour(['file1.txt', 'file2.txt', 'file3.txt'],
                    [' label1', 'label2', label3'],
                    outputpdf.pdf
                    source="PSR B1937+21)

No function return -- plot is saved in current directory. 


Things to do:
    Add arrows for left and right boundary errors
"""

#Plot contours from xspec log
def plotcontour(fname_list, labels, output, source=None):

    #Plot first contour
    log = logfile(fname_list[0])
    contour_params = log.get_contour()

    #Set extent of contour plots
    extent = [min(contour_params.param_x_bins), 
              max(contour_params.param_x_bins),
              min(contour_params.param_y_bins), 
              max(contour_params.param_y_bins)]

    #Color groups for contours
    contour_colors = [ ["#4c167c", "#351394", "#7B1394" ],
                      ["#0135d2", "#0C7DE8", "#0BAFDE" ], 
                      ["#c91c00", "#E00B14", "#E0470B" ],
                      ['#123611', '#287825', '#3DB839' ] ]

    #Set labels as rich text formatting
    labels = [ r'$'+l+r'$' for l in labels ]

    #Set the figure size
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
    plt.subplots_adjust(top=.95)

    #Set the aspect ratio
    aspect=1

    #Plot the 2d array
    im = ax.imshow(contour_params.array, origin='lower', extent=extent, aspect='auto',
                   cmap='Purples')

    #Isolate the min chi2 and mark it with +
    minvals = np.where(contour_params.array == np.min(contour_params.array))
    min_chisq = contour_params.array[ minvals[0][0], minvals[1][0] ]
    ax.scatter([contour_params.param_x_bins[minvals[1][0]]], 
               [contour_params.param_y_bins[minvals[0][0]]], 
               color=contour_colors[0][0], marker='+', s=200, 
               label=labels[0],zorder=5)

    #chisq_contours = scipy.stats.chi2.isf([.1, .05], params.dof)
    #Significance levels from: https://ned.ipac.caltech.edu/level5/Wall2/Wal3_4.html
    chisq_contours = [ min_chisq+val for val in [4.61, 9.21] ]
    #Plot contours
    ax.contour(contour_params.array, chisq_contours, extent=extent, 
               origin='lower', colors=contour_colors[0], alpha=.8)
    ax = plotparams(ax)

    plt.setp(ax.get_yticklabels()[0], visible=False)    
    plt.setp(ax.get_yticklabels()[-1], visible=False)
    plt.setp(ax.get_xticklabels()[0], visible=False)    
    plt.setp(ax.get_xticklabels()[-1], visible=False)
    ax.set_xlim(left=extent[0], right=extent[1])
    ax.set_ylim(bottom=extent[2], top=extent[3])

    ax.set_xlabel(contour_params.xlabel, fontsize=20)
    ax.set_ylabel(contour_params.ylabel, fontsize=20)


    #Add second axis for colorbar to make it easier to move and resize
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05, aspect=30*aspect)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(r'$NICER$ $\chi^2$', fontsize=15, rotation=270, labelpad=20)
    
    #plt.subplots_adjust(left=.08, bottom=.08, right=.94, top=.98)

    #Iterate through the other files to overplot the contours (but not the images)
    for i in range(1, len(fname_list)):

        l = logfile(fname_list[i])
        cp =l.get_contour()
        cc = contour_colors[i]

        mv = np.where(cp.array == np.min(cp.array))
        mchsq = cp.array[ mv[0][0], mv[1][0] ]
        chisqc = [ mchsq + val for val in [4.61, 9.21] ]
        ax.contour(cp.array, chisqc, extent=extent,
                   origin='lower', colors=cc, alpha=.8)
        #If the point is below the ylim bottom, plot at lim with down arrow
        if cp.param_y_bins[mv[0][0]] < extent[2] or cp.param_y_bins[mv[0][0]] == extent[2]:
            ax.scatter([cp.param_x_bins[mv[1][0]]],
                       [extent[2]+.1], marker="$\downarrow$",
                       color=cc[0], s=200, label=labels[i])
        elif cp.param_y_bins[mv[0][0]] > extent[3] or cp.param_y_bins[mv[0][0]] == extent[3]:
            ax.scatter([cp.param_x_bins[mv[1][0]]],
                       [extent[3]-.1], marker=r'$\uparrow$',
                       color=cc[0], s=200, label=labels[i])
        else:
            ax.scatter([cp.param_x_bins[mv[1][0]]], 
                       [cp.param_y_bins[mv[0][0]]], 
                       color=cc[0], marker='+', s=200, zorder=5, label=labels[i])
            print(cp.param_x_bins[mv[1][0]], cp.param_y_bins[mv[0][0]])


    ax.legend(fontsize=15, edgecolor='black', loc=2)

    if source is not None:

       source= process.extract(source, [r'PSR B1821$-$24', 
                                         r'PSR B1937$+$21', r'PSR J0218$+$4232'], 
                                limit=1)[0][0]
       fig.text(.95, .95, source, ha='right', va='top',
                 transform=ax.transAxes, fontsize=20)

    fig.savefig(output, dpi=1000)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-c", help="Contour files", nargs='+', type=str)
    parser.add_argument("-l", help="Labels (in order of contours", nargs='+', type=str)
    parser.add_argument("-o", help="output", default='ContourPlot.pdf', type=str)
    parser.add_argument("-s", help="source label", default=None, type=str)
    args = parser.parse_args()

    plotcontour(args.c, args.l, args.o, args.s)

