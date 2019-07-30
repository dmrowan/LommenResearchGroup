#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import argparse
from astropy import log
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import pexpect
import subprocess
import time
from astropy.table import Table

import niutils

#Dom Rowan 2019

desc = """
Procedure to generate phase resolve spectra

Example command line calls:

    $ genspectra.py --evt eventfile --lp lower_phase --up upper_phase 
                    --save_pha outout.pha --nchan photonbins

    $ genspectra.py --evt J0218.evt --lp .2 --up .4 --save_pha mypha.pha

    #Can use one or two phase regions
    $ genspectra.py --evt B1821.evt --lp .2 .4 --up .7 .9 
                    --save_pha mypha.pha --nchan 750

Example function calls:

    >>> genspectra(eventfile, lowerphase, upperphase, 
                   grphha_channels, save_pha=output.pha
                   save_plot=output.plot)

    >>> genspectra("1821data.evt", .97, 1.06, 1000,
                   save_pha="output.pha")

    #Can input one or two phase ranges
    >>> genspectra("1821data.evt", (.2, .4), (.7, .9),
                   700, save_pha="backgroundspec.pha"
                  
TODO:
    (1) Make rmf and arf import paths work for generating pha file outside current directory


"""

#Make phase selections using fselect
#Input event file .evt, output .fits file
def fselect_phase(evt, output, phase_lower, phase_upper, 
                  clobber=False, verbose=True):
    if verbose:
        log.info("Selecting phase with fselect")
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

#version of phase selection that has two ranges
def fselect_two_phases(evt, output, phase_1, phase_2, clobber=False,
                       verbose=True):

    if verbose:
        log.info("Selecting two phase regions with fselect")
    if os.path.isfile(output) and (not clobber):
        print("File already exists")
        return
    else:
        #Build fselect command
        cmd = [ 'fselect', evt, output ]
        if not (phase_1[0] <= 1 <= phase_1[1]):
            first_phase = f"PULSE_PHASE >= {phase_1[0]} &&" \
                          f"PULSE_PHASE <= {phase_1[1]}"
        else:
            first_phase = f"PULSE_PHASE >= {phase_1[0]} ||" \
                          f"PULSE_PHASE <= {phase_1[1]-1}"

        second_phase = f"PULSE_PHASE >= {phase_2[0]} &&"\
                       f"PULSE_PHASE <= {phase_2[1]}"

        #Connect both phases with or statement
        both_phases = f"({second_phase}) || ({first_phase})"
        cmd.append(both_phases) 
        cmd.append('clobber=yes')
        if verbose:
            print(cmd)
        subprocess.run(cmd)

#Update exposure keyword with fmodhead
def fmodhead_exp(input_name, output, new_exp, verbose=True):
    if verbose:
        log.info("Updating exposure with fmodhead")

    with open('htemp.dat', 'w') as f:
        f.write(f"EXPOSURE {float(new_exp)} / Value of Exposure Changed")

    #update keyword in each header
    cmds = [['fmodhead', f"{input_name}[0]", 'htemp.dat'], 
            ['fmodhead', f"{input_name}[1]", 'htemp.dat'],
            ['fmodhead', f"{input_name}[2]", 'htemp.dat']]
    for cmd in cmds:
        subprocess.run(cmd)

    #Change file name if necessary
    if input_name != output:
        subprocess.run(['cp', input_name, output])


#Set the header keyword for RESP and ANCR
def fparkey_resp_arf(input_name, verbose=True):
    if verbose:
        log.info("Updating resp and arf with fparkey")

    #Check if paths exist
    # (we can't use absolute paths because of fits header restrictions)
    check_resp_arf_path()

    #Run fparkey twice
    cmd = ['fparkey', 'nixtiref20170601v001.rmf', 
           f"{input_name}[1]", 'RESPFILE']

    cmd2 = ['fparkey', 'nixtiaveonaxis20170601v002.arf',
            f"{input_name}[1]", 'ANCRFILE']

    subprocess.run(cmd)
    subprocess.run(cmd2)


#Make copies of resp and arf in current directory
def check_resp_arf_path():
    if not os.path.isfile('nixtiref20170601v001.rmf'):
        cmd = ['cp', '/packages/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf',
               '.']
        subprocess.run(cmd)

    if not os.path.isfile('nixtiaveonaxis20170601v002.arf'):
        cmd2 = ['cp', '/packages/caldb/data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf',
                '.']
        subprocess.run(cmd2)
        
#Find exposure of fits/evt/pha file
def get_exposure(f):
    t = Table.read(f, hdu=1)
    return t.meta['EXPOSURE']

#Double check fselect phase with plot and min/max
def test_fselect(f, plot=False):
    t = Table.read(f, hdu=1)
    print(t['PULSE_PHASE'].min(), t['PULSE_PHASE'].max())
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.hist(t['PULSE_PHASE'], edgecolor='black', color='xkcd:violet')
        ax = niutils.plotparams(ax)
        plt.show()

#Spawn xselect child 
def xselect(lower_e=0, upper_e=1200, lower_phase=0, upper_phase=1, 
            epoch=50000, period=1, datadir=None, eventfile=None,
            session="autopython", verbose=True):

    if verbose:
        log.info("Producing pha in xselect")
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
def xspec_wrapper(phafile, 
                  save=None, verbose=True):
    #Spawn xspec child with pexpect
    if verbose:
        log.info("Running xspec")
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
    xspec.sendline("ig **-0.5, 10.-**")
    xspec.expect("XSPEC12>")
    if save is None:
        xspec.sendline("cpd /xs")
    else:
        xspec.sendline(f"cpd {save}.ps/cps")
    xspec.expect("XSPEC12>")
    xspec.sendline("plot data")
    time.sleep(1)
    

#Pexpect wrapper for grppha to bin spectra
def grppha_wrapper(pha_in, pha_out, nchan, verbose=True):
    #Spawn grppha child with pexpect
    if nchan <= 0:
        log.info("No grppha perfomed")
        cmd = ['cp', pha_in, pha_out]
    elif verbose:
        log.info("Grouping energy bins with grppha")

    grppha = pexpect.spawn("grppha")
    grppha.expect("Please enter PHA filename")
    grppha.sendline(f"{pha_in}")
    grppha.expect("Please enter output filename")
    grppha.sendline(f"{pha_out}")
    grppha.expect("GRPPHA")
    #grppha.sendline(f"group 0 {1499-nchan} {nchan}")
    grppha.sendline(f"group min {nchan}")
    grppha.expect("GRPPHA")
    grppha.sendline(f"exit !{pha_out}")
    grppha.wait()
    grppha.close()

    cmd = ['fparkey', str(nchan), f"{pha_out}[1]", 'GROUPING', 'add=yes']
    subprocess.run(cmd)
    

#Conver ps to pdf
def convertPDF(psfile, display=False, verbose=True):
    if verbose:
        log.info("Converting to pdf")
    cmd = ['ps2pdf', psfile, psfile.replace('.ps', '.pdf')]
    subprocess.run(cmd)
    if display:
        cmd2 = subprocess.run(['xdg-open', psfile.replace('.ps', '.pdf')])
                                                                
#Wrap heasarc calls to generate spectra with energy/phase selections
def gen_spectra(evt, phase_lower, phase_upper, 
                nchan, save_pha=None, 
                save_plot=None, display=False,
                run_xspec=True, verbose=True):

    if all( [type(p) not in [list, tuple, np.ndarray] 
             for p in [phase_lower, phase_upper]]):
        using_two_regions = False
    elif len(phase_lower) != 1 and len(phase_upper) != 1:
        using_two_regions = True
    else:
        phase_lower = phase_lower[0]
        phase_upper = phase_upper[0]
        using_two_regions = False

    #Fiind exposure from event file
    t = Table.read(evt, hdu=1)
    original_exp = t.meta['EXPOSURE']
    #Calculate exposure corresponding to phase selection
    if using_two_regions:
        new_exp = original_exp * ( (phase_upper[1]-phase_upper[0]) 
                                  +(phase_lower[1]-phase_lower[0]) )
    else:
        new_exp = original_exp * (phase_upper-phase_lower)

    #Apply phase selection to evt file
    if using_two_regions:
        fselect_two_phases(evt, "autofits.fits", 
                           phase_lower, phase_upper, 
                           clobber=True, verbose=verbose)
    else:
        fselect_phase(evt, "autofits.fits", 
                      phase_lower, phase_upper, clobber=True, 
                      verbose=verbose)

    #test_fselect("autofits.fits", plot=False)

    #Run xselect to create pha file
    xselect(datadir='./', eventfile="autofits.fits", 
            session='autoxselect', verbose=verbose)


    #Modify the exposure keyword to correspond to the phase selection
    fmodhead_exp("autoxselect_spec.pha", "autoxselect_spec_2.pha",
             new_exp, verbose=verbose)

    #Modify the resp and arf keywords
    fparkey_resp_arf("autoxselect_spec_2.pha", verbose=verbose)

    #Apply grppha to set the grouping keyword
    grppha_wrapper("autoxselect_spec_2.pha", "autoxselect_spec_grppha.pha", 
                   nchan, verbose=verbose)

    #Save pha to output name
    if save_pha is not None:
        subprocess.run(['cp', 'autoxselect_spec_grppha.pha', save_pha])

    #Run xspec
    if run_xspec:
        xspec_wrapper('autoxselect_spec_grppha.pha', 
                      save=save_plot,
                      verbose=verbose)

        if save_plot is not None:
            convertPDF(f"{save}.ps", display=display, verbose=verbose)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--pha", 
                        help="Spectra file path for spec", 
                        type=str, default=None)
    parser.add_argument("--evt", 
                        help="Event file path for lc", 
                        type=str, default=None)
    parser.add_argument("--lp", 
                        help="Lower phase for xselect filtering", 
                        nargs='+',
                        type=float, default=0)
    parser.add_argument("--up", 
                        help="Upper phase for xselect filtering", 
                        nargs='+',
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


    gen_spectra(args.evt, args.lp, args.up, 
                args.nchan, save_pha=args.save_pha, 
                save_plot=args.save_plot, display=args.display)



