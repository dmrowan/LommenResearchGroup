#!/usr/bin/env python

import argparse
import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import Table
from astropy import log
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patheffects import withStroke
import pandas as pd
from tqdm import tqdm

from background_LC import convert_time

desc="""
Generate a skymap showing the BKGD_RXTE regions
"""

rc('text', usetex=True)

def find_exposure():

    #Open list of event files
    with open('/students/pipeline/heasoft6.26/bkgd_all_evt_list', 'r') as f:
        paths = f.readlines()

    #Grab BKGD_RXTE number
    bkgd_number = [ p[10] for p in paths ]
    paths = [ '/students/pipeline/heasoft6.26/'+p.strip('\n') for p in paths ]

    #Make dictionary
    exp_totals = {'1':0, '2':0, '3':0, '4':0, '5':0, '6':0, '8':0}

    full_exp_list = []
    #Pull exp val from each table
    for i in tqdm(range(len(paths))):
        """
        tab = Table.read(paths[i], hdu=1)
        exp = tab.meta['EXPOSURE']
        if exp > 1e8: print(paths[i])
        full_exp_list.append(exp)
        """
        hdul = fits.open(paths[i])
        tstart = convert_time(hdul[0].header['TSTART'])
        tstop = convert_time(hdul[0].header['TSTOP'])
        exp = (tstop - tstart).total_seconds()
        exp_totals[bkgd_number[i]] += exp

    return exp_totals

def count_obs():

    #Open list of event files
    with open('/students/pipeline/heasoft6.26/bkgd_all_evt_list', 'r') as f:
        paths = f.readlines()

    bkgd_number = [ p[10] for p in paths ]

    bkgd_number.sort()

    obsID_totals = dict(
            zip(set(bkgd_number), 
                [bkgd_number.count(k) for k in set(bkgd_number)]))

    return obsID_totals

def make_table(output):

    #Background numbers
    bkgd_n = [ str(i) for i in range(1, 7) ] + ['8']
    #Get number of observations in each region
    obs_n_dic = count_obs()
    #Find exposure
    exp_dic = find_exposure()

    #Coordinates from nasa data access site
    coord_list = [ (5, -67), (60, 2), (138, 15), (30, 10), (345-360, -18), 
               (160, 72.57), (183.7-360, 53.3) ]

    #Create the dataframe
    df = pd.DataFrame({
        'Region': ['BKGD\_RXTE\_1']+bkgd_n[1:],
        'Number of ObsIDs': [ obs_n_dic[k] for k in obs_n_dic.keys() ],
        'Total Exposure (ks)': [ int(round(exp_dic[k]/(1e3))) 
                                for k in exp_dic.keys() ],
        'Right Ascension': [ r'${0}^\circ$'.format(str(c[0])) 
                             for c in coord_list ],
        'Declination': [ r'${0}^\circ$'.format(str(c[1]))
                         for c in coord_list ]})

    #First write to latex table format
    with open(output, 'w') as f:
        f.write(df.to_latex(index=False, escape=False, column_format='rccrr'))
        f.close()

    #Read in the file again
    with open(output, 'r') as f:
        lines = f.readlines()
        f.close()

    #Insert a mid line and total at the bottom of the table
    i = 11
    lines.insert(i, '\midrule\n')
    lines.insert(i+1, 'Total: & {0} & {1} & &\\\ \n'.format(
        sum(df['Number of ObsIDs']), 
        int(round(sum([ exp_dic[k]/(1e3) for k in exp_dic.keys() ])))))

    #Save the file again
    with open(output, 'w') as f:
        for l in lines:
            f.write(l)


def plot_skymap(savefig=None):
    labels = ['1', '2', '3', '4', '5', '6', '8' ]
    #Coordinates from nasa target list site
    coord_list = [ (5, -67), (60, 2), (138, 15), (30, 10), (345-360, -18), 
               (160, 72.57), (183.7-360, 53.3) ]

    #Use sky cood to grab info
    ra = coord.Angle(np.array([c[0] for c in coord_list]), unit=u.degree)
    ra.wrap_at(180*u.degree)
    dec = coord.Angle(np.array([c[1] for c in coord_list]), unit=u.degree)

    #Setup plot with projection
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), 
                           subplot_kw={'projection': "mollweide"})
    ax.grid(True, zorder=0)
    ax.set_axisbelow(True)
    plt.subplots_adjust(top=.98, right=.98, bottom=.02, left=.04)

    log.info("Calculating Exposure Totals")
    #Pull exposure values
    exp_totals = find_exposure()

    #Rescale exposure to point size
    key_max = max(exp_totals.keys(), key=(lambda k: exp_totals[k]))
    scalar = 600
    exp_scaled = [ (exp_totals[k]/exp_totals[key_max])*scalar
                   for k in exp_totals.keys() ]

    #Plot the points
    ax.scatter(ra.radian, dec.radian, s=exp_scaled, color='#003f5c', zorder=1)
    #Add text labels
    myeffectw = withStroke(foreground="black", linewidth=1)
    txtkwargsw = dict(path_effects=[myeffectw])
    afont = {'fontname':'monospace'}
    for i in range(len(labels)):
        if labels[i] == '8':
            rai = ra[i].radian + 10*np.pi/180
        elif labels[i] != '6':
            rai = ra[i].radian + 7*np.pi/180
            va='bottom'
        else:
            rai = ra[i].radian - 26*np.pi/180
            va='top'
        ax.text(rai, dec[i].radian, labels[i],
                ha='left', va=va, weight='normal', fontsize=12,
                zorder=2, color='black', **txtkwargsw, **afont)



    #Plot mikly way plane
    lvals = np.linspace(0, 360, 10000)
    gal = coord.SkyCoord(l=lvals, b=0, 
                         frame='galactic', unit=u.degree)
    mw = gal.transform_to('icrs')

    ax.scatter(mw.ra.wrap_at(180*u.degree).radian, mw.dec.radian, color='gray', alpha=.6, s=2)


    #This is literally the worst
    fig.canvas.draw()
    ylabels = [ item.get_text() for item in ax.get_yticklabels() ]
    new_y_labels = [ s.replace('°', r'$^{\circ}$') for s in ylabels ]
    ax.set_yticklabels(new_y_labels)

    xlabels = [ item.get_text() for item in ax.get_xticklabels() ]
    new_x_labels = [ s.replace('°', r'$^{\circ}$') for s in xlabels ]
    ax.set_xticklabels(new_x_labels)
    
    #Save or show
    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig)
        #if savefig.endswith('pdf'):
        #    fig.savefig(savefig.replace('pdf', 'png'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--skymap", 
                        help='Create skymap of BKGD RXTE observations', 
                        default=False, action='store_true')
    parser.add_argument("--table",
                        help='Create table with BKGD RXTE info',
                        default=False, action='store_true')
    args = parser.parse_args()

    if args.skymap:
        plot_skymap('skymap.pdf')

    if args.table:
        make_table('bkgd_rxte_table.tex')
