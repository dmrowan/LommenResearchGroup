#!/usr/bin/env python

from astropy import log
from astropy.table import Table
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

from niutils import plotparams

#Dom Rowan 2020

#Show and plot columns from attitude file
def attitude_example(att_table):
    log.info('Print Attitude Table')
    print(att_table)
    
    print('ST_VALID values: ', list(set(att_table['ST_VALID'])))
    print('SUBMODE_AZ values', list(set(att_table['SUBMODE_AZ'])))
    print('SUBMODE_EL values', list(set(att_table['SUBMODE_EL'])))
    
    log.info("Generating Plot")
    
    colors=["#cb6a49", "#a46cb7", "#7aa457"]

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax = plotparams(ax)
    ax.plot(att_table['TIME'], att_table['SUBMODE_EL'], color=colors[0], label='SUBMODE_EL')
    ax.plot(att_table['TIME'], att_table['SUBMODE_AZ'], color=colors[1], label='SUBMODE_AZ')
    ax.plot(att_table['TIME'], att_table['ST_VALID'], color=colors[2], label='ST_VALID')
    ax.legend(fontsize=12, edgecolor='black', loc='upper right')
    ax.set_xlabel('Time (s)', fontsize=15)
    ax.set_ylabel('Parameter Flag', fontsize=15)
    plt.show()
   

#Show catalog table
def catalog_example(cat_table):
    log.info("Print Catalog File")
    print(cat_table)
    
    
#Show and plot columns from filter file
def filter_example(mkf_table):
    print("List of all the columns:")
    print(mkf_table.colnames)
    
    log.info("Print Filter File")
    print(mkf_table)

    
    log.info("Generating Plot")
    fig, ax = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
    plt.subplots_adjust(hspace=0)
    colors=["#cb6a49", "#a46cb7", "#7aa457"]

    ax[0].plot(mkf_table['TIME'], mkf_table['SUN_ANGLE'], color=colors[0], label='Sun Angle')
    ax[0].plot(mkf_table['TIME'], mkf_table['MOON_ANGLE'], color=colors[1], label='Moon Angle')
    ax[0].plot(mkf_table['TIME'], mkf_table['ELV'], color=colors[2], label='Earth Elevation')

    ax[0].set_ylabel('Angle (Degrees)', fontsize=15)

    ax[1].plot(mkf_table['TIME'], mkf_table['SAT_LAT'], color=colors[0], label='Satellite Lattitude')
    ax[1].plot(mkf_table['TIME'], mkf_table['SAT_LON'], color=colors[1], label='Satellite Longitude')
    ax[1].set_ylabel('Degrees', fontsize=15)

    ax[2].plot(mkf_table['TIME'], mkf_table['COR_SAX'], color=colors[0])
    ax[2].set_ylabel('Cutoff Rigidity (GeV/c)', fontsize=15)


    for a in ax:
        a = plotparams(a)
        if a.get_legend_handles_labels() != ([], []):
            a.legend(loc='upper right', fontsize=12, edgecolor='black')    
    plt.show()
    
#Show and plot columns from orbit file
def orbit_example(orb_table):
    log.info("Print Orbit File")
    print(orb_table)
    
    log.info("Plot orbit position")
    fig, ax = plt.subplots(1, 1, figsize=(6,6), subplot_kw={'projection':'3d'})
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)

    #Generating 3d plot of position, colored by time
    ax.scatter3D(orb_table['X'], orb_table['Y'], orb_table['Z'], c=np.arange(len(orb_table)), 
                 cmap=plt.get_cmap('jet'), s=5)
    plt.show()

#Show housekeeping table
def housekeeping_example(hk_table):
    log.info("Printing housekeeping file")
    print(hk_table)
    
    
if __name__ == '__main__':
    att_file = 'ni1070020426.att'
    cat_file = 'ni1070020426.cat'
    mkf_file = 'ni1070020426.mkf'
    orb_file = 'ni1070020426.orb'
    hk_file = 'ni1070020426_0mpu0.hk'
    table_list = [ Table.read(fname, hdu=1) for fname in [att_file, cat_file, mkf_file, orb_file, hk_file] ]
    att_table, cat_table, mkf_table, orb_table, hk_table = table_list
    
    
    attitude_example(att_table)
    filter_example(mkf_table)
    orbit_example(orb_table)
    housekeeping_example(hk_table)
