#!/usr/bin/env python

from astropy import log
import numpy as np
import os
import subprocess

# Dom Rowan 2020

desc="""
A wrapper for calling the 3C50 background estimator on BKGD observations

While it would be possible to merge the UFAs from all observations together, it is less computationally intensive to generate bkgd for each observation and merge those
"""

def main(dir_head='/students/pipeline/heasoft6.26/',
         bkgd_list='bkgd_all_evt_list',
         outdir='./'):

    bkgd_list = os.path.join(dir_head, bkgd_list)
    #Get list of observations
    with open(bkgd_list) as f:
        paths = f.readlines()

    #Find paths to BKGD directories
    paths = [p.strip("_pipe/cleanfilt.evt\n") for p in paths]
    rootdir = [os.path.join(dir_head, f'BKGD_RXTE_{p[10]}')
               for p in paths]

    #Iterate through paths
    #for i in range(len(paths)):
    for i in range(10):
        bkgd_n = paths[i][10]
        obsID = paths[i].lstrip(f"BKGD_RXTE_{bkgd_n}").replace('/', '')

        #See if background already exists
        if os.path.isfile(f'3C50_{obsID}_tot.pha'):
            log.info("Background already generated for {obsID}")
            continue
        else:
            log.info(f"Running nibackgen on {obsID}")
            cmd = ['nibackgen3C50', f"rootdir={rootdir[i]}", f'obsid={obsID}', 
                   'bkgdir=/packages/heasoft-6.26.1/backdir/bg_model_3C50_RGv5/',
                   f"totspec=3C50_{obsID}_tot.pha", 
                   f"bkgspec=3C50_{obsID}_bkg.pha"]
            subprocess.call(cmd)

def merge():

    all_files = os.listdir()
    bkg_paths = []
    for p in all_files:
        if p.endswith('_bkg.pha'):
            bkg_paths.append(p)
        else:
            continue

    groups = np.array_split(np.array(bkg_paths), 30)
    for i in range(len(groups)):
        with open(f'bkg_list_{i}', 'w') as f1:
            for item in groups[i]:
                if os.path.isfile(item):
                    f1.write(item+'\n')
                else:
                    print(item, " not found")
        
        cmd1 = ['addspec', f'bkg_list_{i}', f'merged_bkg_{i}', 
                'qaddrmf=FALSE', 'qsubback=FALSE']
        subprocess.run(cmd1)

    with open('bkg_big_list_for_merge.txt', 'w') as f2:
        for i in range(len(groups)):
            f2.write(f'merged_bkg_{i}.pha\n')

    cmd2 = ['addspec', 'bkg_big_list_for_merge.txt', 'merged_3C50_bkg',
            'qaddrmf=FALSE', 'qsubback=FALSE']
    subprocess.call(cmd2)

if __name__ == "__main__":
    #main()
    merge()
