#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
import os
import pexpect
import sys
import subprocess
import time
from astropy.table import Table
import spectraplots
#Dom Rowan 2019

desc = """
Procedure to generate spectra for pulsar in phase region
"""

#Make phase selections using fselect
#Input event file .evt, output .fits file
def fselect_phase(evt, output, phase_lower, phase_upper, clobber=False):
    if os.path.isfile(output) and (not clobber):
        print("File already exists")
        return
    else:
        cmd = [ 'fselect', evt, output]
        #If the phase region overlaps 0 (or 1), use an or statement
        if not (phase_lower <= 1 <= phase_upper):
            cmd.append(
                    f"PULSE_PHASE >= {phase_lower} &&" \
                    f"PULSE_PHASE <= {phase_upper}")
        else:
            cmd.append(
                    f"PULSE_PHASE >= {phase_lower} ||" \
                    f"PULSE_PHASE <= {phase_upper-1}")
        cmd.append('clobber=yes')
        subprocess.run(cmd)

#Update exposure keyword to reflect phase selection
def update_exp(fits, output, new_exp, hdu=1):
    t = Table.read(fits, hdu=hdu)
    t.meta['EXPOSURE'] = new_exp
    t.write(output, overwrite=True)
        
#Find exposure of fits/evt/pha file
def get_exposure(f):
    t = Table.read(f, hdu=1)
    return t.meta['EXPOSURE']

#Double check fselect phase with plot and min/max
def test_fselect(fits, plot=False):
    t = Table.read(fits, hdu=1)
    print(t['PULSE_PHASE'].min(), t['PULSE_PHASE'].max())
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.hist(t['PULSE_PHASE'], edgecolor='black', color='xkcd:violet')
        plt.show()

#Spawn xselect child 
def xselect(lower_e=0, upper_e=1200, lower_phase=0, upper_phase=1, 
            epoch=50000, period=1, datadir=None, eventfile=None,
            session="autopython"):

    assert(all([type(val) in [int, float] for val in [
                lower_phase, upper_phase, 
                epoch, period]]))
    assert(type(lower_e) == int)
    assert(type(upper_e) == int)
    assert(lower_phase >= 0)
    assert(upper_phase <= 1)
    assert(all([type(val) == str for val in [datadir, eventfile, session]]))
    assert(os.path.isfile(f"{datadir}{eventfile}"))

    phase = f"{lower_phase}-{upper_phase}"

    #Spawn xsel child
    xsel = pexpect.spawn("xselect")
    xsel.expect("Enter session name >")
    xsel.sendline(session)
    xsel.expect(f"{session}:SUZAKU")
    xsel.sendline(f"set datadir {datadir}")
    xsel.expect(f"{session}:SUZAKU")
    xsel.sendline(f"read event {eventfile}")
    xsel.expect("Reset the mission")
    xsel.sendline("yes")
    xsel.expect(f"{session}:NICER-XTI-PHOTON")
    #Perform energy filter
    if not (lower_e==0 and upper_e==1200):
        print("Performing energy filter in xselect")
        xsel.sendline(f"filter PHA_CUTOFF {lower_e} {upper_e}")
        xsel.expect(f"{session}:NICER-XTI-PHOTON")
    #Perform phase filter
    if not (lower_phase==0 and upper_phase==1):
        print("Performing phase filter in xselect")
        xsel.sendline(f"filter phase {epoch} {period} {phase}")
        xsel.expect(f"{session}:NICER-XTI-PHOTON")
    xsel.sendline(f"extract spec")
    xsel.expect(f"{session}:NICER-XTI-PHOTON")
    xsel.sendline(f"save spec {session}_spec")
    #Save pha file using session name
    if os.path.isfile(f"{session}_spec.pha"):
        xsel.expect("overwrite")
        xsel.sendline("yes")
    xsel.expect(f"{session}:NICER-XTI-PHOTON")
    xsel.sendline("exit no")
    xsel.close()

#Pexpect wrapper for xspec to generate xspec plot
def xspec_wrapper(phafile, channel_lower, channel_upper, save=None):
    #Spawn xspec child with pexpect
    xspec = pexpect.spawn("xspec")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"data 1:1 {phafile}")
    xspec.expect("XSPEC12>")
    xspec.sendline("resp 1 /packages/caldb/data/nicer/xti/cpf/" \
            "rmf/nixtiref20170601v001.rmf")
    xspec.expect("XSPEC12>")
    xspec.sendline("arf 1 /packages/caldb/data/nicer/xti/cpf/" \
            "arf/nixtiaveonaxis20170601v002.arf")
    xspec.expect("XSPEC12>")
    xspec.sendline("ig **-0.3, 10.-**")
    xspec.expect("XSPEC12>")
    if save is None:
        xspec.sendline("cpd /xs")
    else:
        xspec.sendline(f"cpd {save}.ps/cps")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"ig {channel_upper}-**")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"ig **-{channel_lower}")
    xspec.expect("XSPEC12>")
    xspec.sendline("plot data")
    time.sleep(1)
    

#Pexpect wrapper for grppha to bin spectra
def grppha_wrapper(pha_in, pha_out, nchan):
    #Spawn grppha child with pexpect
    grppha = pexpect.spawn("grppha")
    grppha.expect("Please enter PHA filename")
    grppha.sendline(f"{pha_in}")
    grppha.expect("Please enter output filename")
    grppha.sendline(f"{pha_out}")
    grppha.expect("GRPPHA")
    grppha.sendline(f"group min {nchan}")
    grppha.expect("GRPPHA")
    grppha.sendline(f"exit !{pha_out}")
    grppha.wait()
    grppha.close()
    

#Conver ps to pdf
def convertPDF(psfile, display=False):
    cmd = ['ps2pdf', psfile, psfile.replace('.ps', '.pdf')]
    subprocess.run(cmd)
    if display:
        cmd2 = subprocess.run(['xdg-open', psfile.replace('.ps', '.pdf')])
                                                                
#Wrap heasarc calls to generate spectra with energy/phase selections
def gen_spectra(evt, phase_lower, phase_upper, 
                channel_lower, channel_upper, nchan,
                save_pha=None, save_plot=None, display=False):

    t = Table.read(evt, hdu=1)
    original_exp = t.meta['EXPOSURE']
    new_exp = original_exp * (phase_upper-phase_lower)

    print("-----Running fselect-----")
    fselect_phase(evt, "autofits.fits", 
                  phase_lower, phase_upper, clobber=True)

    #test_fselect("autofits.fits", plot=True)

    print("-----Running xselect-----")
    xselect(datadir='./', eventfile="autofits.fits", 
            session='autoxselect')

    print("-----Updating exposure-----")
    update_exp("autoxselect_spec.pha", "autoxselect_spec_2.fits", 
               new_exp)
    subprocess.run(['mv', 'autoxselect_spec_2.fits',
                    'autoxselect_spec_2.pha'])

    print("-----Running grppha-----")
    grppha_wrapper("autoxselect_spec_2.pha", "autoxselect_spec_grppha.pha", nchan)

    if save_pha is not None:
        subprocess.run(['cp', 'autoxselect_spec_grppha.pha', save_pha])
    print("-----Running xspec-----")
    xspec_wrapper('autoxselect_spec_grppha.pha', 
                  channel_lower, channel_upper, save=save_plot)

    if save_plot is not None:
        print("-----Converting to PDF-----")
        convertPDF(f"{save}.ps", display=display)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--pha", 
                        help="Spectra file path for spec", 
                        type=str, default=None)
    parser.add_argument("--evt", 
                        help="Event file path for lc", 
                        type=str, default=None)
    parser.add_argument("--le", 
                        help="Lower PHA_CUTOFF for xselect filtering",
                        type=int, default=0)
    parser.add_argument("--ue", 
                        help="Upper PHA_CUTOFF for xselect filtering",
                        type=int, default=1200)
    parser.add_argument("--lp", 
                        help="Lower phase for xselect filtering", 
                        type=float, default=0)
    parser.add_argument("--up", 
                        help="Upper phase for xselect filtering", 
                        type=float, default=1)
    parser.add_argument("--nchan", 
                        help="Number of channels for grppha",
                        default=1000, type=int)
    parser.add_argument("--save_plot", 
                        help="Save plot", 
                        default=None, type=str)
    parser.add_argument("--save_pha", 
                        help="Save pha file", 
                        default=None, type=str)
    parser.add_argument("--display", 
                        help="Display saved plot after converting to pdf", 
                        default=False, action='store_true')

    args = parser.parse_args()

    gen_spectra(args.evt, args.lp, args.up, args.le, args.ue, 
                args.nchan, save_pha=args.save_pha, 
                save_plot=args.save_plot, display=args.display)

