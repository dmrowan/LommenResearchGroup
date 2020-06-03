#!/usr/bin/env python

from astropy import log
import os
from nicergof.bkg import bkg_estimator as be

def gen_br_slices():
    pha_list = ['br_earth_--40.pha', 'br_earth_40--60.pha',
                'br_earth_60--80.pha', 'br_earth_80--180.pha', 
                'br_earth_180--.pha']

    mkf_list = [ p.replace('pha', 'mkf') for p in pha_list ]

    for mkf in mkf_list:
        if os.path.isfile(mkf+'3'):
            continue
        else:
            be.add_kp(mkf)

    for i in range(len(pha_list)):
        if os.path.isfile(pha_list[i].replace('.pha', '_bkg.pha')):
            continue
        log.info(f"Generating BKG for {pha_list[i]}")
        bkg_chan, bkgspectot, btotexpo = be.mk_bkg_spec_evt(
                pha_list[i], mkf3file=mkf_list[i]+'3')


def gen_all():
    mkf_all = 'bkgd_merged.mkf3'
    pha_all = 'all_events_spectra.pha'
    bkg_chan, bkcspectot, btotexpo = be.mk_bkg_spec_evt(pha_all, mkf_all)

if __name__ == '__main__':
    gen_all()
