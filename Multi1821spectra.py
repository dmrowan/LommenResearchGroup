#!/usr/bin/env python
import spectraplots
import genspectra
import pexpect
import time

#Dom Rowan 2019

print("*****Generating Off-Peak Spectra*****")
genspectra.gen_spectra("../PSR_B1821-24_combined.evt", 0.1, 0.4, 0, 1200, 
                       1200, save_pha="offpeak.pha")

print("*****Generating On-peak spectra*****")
#Min counts should probably get smaller with narrower phase selections
cts = [1200, 1000, 800, 600]
#Define phasebins kinda manually, iterate through them
for i, tup in enumerate([(.97, 1.06), (.98, 1.05), (.99, 1.04), (0, .02)]):
    #This generates the pha file we normally open in xspec manually
    genspectra.gen_spectra("../PSR_B1821-24_combined.evt", tup[0], tup[1],
                           0, 1200, cts[i], save_pha=f"onpeak_{i}.pha")

    #This is opening xspec in python and doing basic fitting 
    xspec = pexpect.spawn("xspec")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"data 1:1 onpeak_{i}.pha")
    xspec.expect("XSPEC12>")
    #This is an xspec script that loads data and fits model
    xspec.sendline("@autofitting.xcm")
    xspec.expect("XSPEC12>")
    #Save the model params with this command
    xspec.sendline(f"save model onpeak_model_{i}")
    xspec.expect("XSPEC12>")
    #Set the xaxis to be keV for our plots (and txt files)
    xspec.sendline("setplot energy")
    xspec.expect("XSPEC12>")
    #Plot both at same time then save txt file
    xspec.sendline(f"plot ufspec delchi")
    xspec.expect("XSPEC12>")
    xspec.sendline("ipl")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"wdata data_onpeak_{i}.txt")
    xspec.expect("XSPEC12>")
    xspec.sendline("exit")



print("*****Generating interpulse spectra*****")
for i, tup in enumerate([(.50, .62), (.51, .61), (.52, .60), (.53, .59)]):
    genspectra.gen_spectra("../PSR_B1821-24_combined.evt", tup[0], tup[1],
                           0, 1200, cts[i], save_pha=f"interpulse_{i}.pha")

    xspec = pexpect.spawn("xspec")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"data 1:1 interpulse_{i}.pha")
    xspec.expect("XSPEC12>")
    xspec.sendline("@autofitting.xcm")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"save model interpulse_model_{i}")
    xspec.expect("XSPEC12>")
    xspec.sendline("setplot energy")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"plot ufspec delchi")
    xspec.expect("XSPEC12>")
    xspec.sendline("ipl")
    xspec.expect("XSPEC12>")
    xspec.sendline(f"wdata data_interpulse_{i}.txt")
    time.sleep(3)
    xspec.expect("XSPEC12>")
    xspec.sendline("exit")


