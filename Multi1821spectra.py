#!/usr/bin/env python
import spectraplots
import genspectra
import pexpect
import time
import os

#Dom Rowan 2019

def autofitting(fname_head):
    xspec = pexpect.spawn("xspec")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"data 1:1 {fname_head}.pha")
    xspec.expect("XSPEC12>")
    #This is an xspec script that loads data and fits model
    xspec.sendline("@autofitting.xcm")
    xspec.expect("XSPEC12>")
    #Save the model params with this command
    xspec.sendline(f"save model model_{fname_head}")
    if os.path.isfile(f"model_{fname_head}.xcm") and clobber:
        xspec.sendline("y")
    xspec.expect("XSPEC12>")
    #Set the xaxis to be keV for our plots (and txt files)
    xspec.sendline("setplot energy")
    xspec.expect("XSPEC12>")
    #Plot both at same time then save txt file
    xspec.sendline(f"plot ufspec delchi")
    xspec.expect("XSPEC12>")
    xspec.sendline("ipl")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"wdata data_{fname_head}.txt")
    if os.path.isfile(f"data_{fname_head}.txt") and clobber:
        xspec.expect("XSPEC12>")
        xspec.sendline("yes")
    xspec.expect("XSPEC12>")
    xspec.sendline("exit")


def main(clobber=True):
    offpeak_range = (0.1, 0.4)
    mincounts_offpeak = 1200 #Not sure how much this matters
    onpeak_ranges = [(.97, 1.06), (.98, 1.05), (.99, 1.04), (0, .02)]
    mincounts_onpeak = [1200, 1000, 800, 600]
    mincounts_interpulse = [1200, 1000, 800, 600]
    print("*****Generating Off-Peak Spectra*****")
    genspectra.gen_spectra("../PSR_B1821-24_combined.evt", 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           mincounts_offpeak, 
                           save_pha="offpeak.pha")

    print("*****Generating On-peak spectra*****")
    #Define phasebins kinda manually, iterate through them
    for i, tup in enumerate(onpeak_ranges):
        #This generates the pha file we normally open in xspec manually
        genspectra.gen_spectra("../PSR_B1821-24_combined.evt", tup[0], tup[1],
                               0, 1200, mincounts_onpeak[i], save_pha=f"onpeak_{i}.pha")

        #This is opening xspec in python and doing basic fitting 
        autofitting(f"onpeak_{i}")


    print("*****Generating interpulse spectra*****")
    for i, tup in enumerate([(.50, .62), (.51, .61), (.52, .60), (.53, .59)]):
        genspectra.gen_spectra("../PSR_B1821-24_combined.evt", tup[0], tup[1],
                               0, 1200, mincounts_interpulse[i], save_pha=f"interpulse_{i}.pha")

        autofitting(f"interpulse_{i}")


if __name__ == '__main__':
    main()
