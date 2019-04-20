#!/usr/bin/env python
import spectraplots
import genspectra
import pexpect
import time
import os
import pandas as pd

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


def main(evt, clobber=True):
    offpeak_range = (0.1, 0.4)
    mincounts_offpeak = 1200 #Not sure how much this matters
    onpeak_ranges = [(.97, 1.06), (.98, 1.05), (.99, 1.04), (0, .02)]
    mincounts_onpeak = [1200, 1000, 800, 600]
    mincounts_interpulse = [1200, 1000, 800, 600]
    print("*****Generating Off-Peak Spectra*****")
    genspectra.gen_spectra(evt, 
                           offpeak_range[0], offpeak_range[1], 
                           0, 1200, 
                           mincounts_offpeak, 
                           save_pha="offpeak.pha")

    print("*****Generating On-peak spectra*****")
    #Define phasebins kinda manually, iterate through them
    for i, tup in enumerate(onpeak_ranges):
        #This generates the pha file we normally open in xspec manually
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts_onpeak[i], 
                               save_pha=f"onpeak_{i}.pha")

        #This is opening xspec in python and doing basic fitting 
        autofitting(f"onpeak_{i}")


    print("*****Generating interpulse spectra*****")
    for i, tup in enumerate([(.50, .62), (.51, .61), (.52, .60), (.53, .59)]):
        genspectra.gen_spectra(evt, tup[0], tup[1],
                               0, 1200, mincounts_interpulse[i], 
                               save_pha=f"interpulse_{i}.pha")

        autofitting(f"interpulse_{i}")


#This class reads in the txt file output
class xspecdata:
    def __init__(self, filename):
        assert(os.path.isfile(filename))
        with open(filename, 'r') as h:
            lines = h.readlines()
        breakidx = np.where(np.array(lines) == 'NO NO NO NO NO\n')[0][0]
        df0 = pd.read_csv(filename, skiprows=3, delimiter=" ", header=None,
                          nrows=breakidx-3)
        df1 = pd.read_csv(filename, skiprows=breakidx+1, 
                          delimiter=" ", header=None)
        df0.columns = ['energy', 'energy_err', 
                       'counts', 'counts_err', 'model']
        df1.columns = ['energy', 'energy_err', 
                       'delchi', 'delchi_err', 'model']
        self.data = df0
        self.residuals = df1

        self.width = None
        self.lower = None
        self.upper = None
        self.component = None
    def set_component(self, component):
        self.component = process.extract(component, ['primary', 'interpulse'],
                                         limit=1)[0][0]
    def set_phaserange(self, p1, p2):
        self.width=p2-p1
        self.lower = p1
        self.upper = p2
    def get_label(self):
        if self.lower is None or self.upper is None:
            print("No phase region specified")
            return -1
        else:
            return f"Phase: {self.lower} -- {self.upper}"

def multi_ufspec(sourcename, primarytxts, interpulsetxts, p_ranges, i_ranges):

    #Init figure
    fig = plt.figure(figsize=(10, 11))
    plt.subplots_adjust(top=.98, right=.98, hspace=.15, left=.15)
    outer = gridspec.GridSpec(2, 1, height_ratios=[1,1])

    inner_p = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                               hspace=0, height_ratios=[3, 1])
    inner_i = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1],
                                               hspace=0, height_ratios=[3,1])
    axp1 = plt.Subplot(fig, inner_p[1])
    axp0 = plt.Subplot(fig, inner_p[0], sharex=axp1)
    axi1 = plt.Subplot(fig, inner_i[1])
    axi0 = plt.Subplot(fig, inner_i[0], sharex=axi1)
    

    primary_data = []
    interpulse_data = []
    for i in range(len(primarytxts)):
        xd = xspecdata(primarytxts[i])
        xd.set_phaserange(p_ranges[i][0], p_ranges[i][1])
        primary_data.append(xd)
    for i in range(len(interpulsetxts)):
        xd = xspecdata(interpulsetxts[i])
        xd.set_phaserange(i_ranges[i][0], i_ranges[i][1])
        interpulse_data.append(xd)
    alldata = [primary_data, interpulse_data]

    #Match sourcename
    sourcename = process.extract(sourcename, 
                                 ['PSR_B1821-24', 'PSR_B1937+21'],
                                 limit=1)[0][0]
    assert(sourcename in ['PSR_B1821-24', 'PSR_B1937+21'])


    labels=["Primary Pulse", "Interpulse"]


    #Plot data
    colors = ["#0096cf",
            "#a7cb20",
            "#6800a2",
            "#71b07b",
            "#ff355f"]
    labels=["Primary Pulse", "Interpulse"]
    for i, ax in enumerate([axp0, axi0]):
        for j in range(len(alldata[i])):
            ax.errorbar(alldata[i][j].data['energy'], 
                        alldata[i][j].data['counts'],
                        xerr = alldata[i][j].data['energy_err'],
                        yerr = alldata[i][j].data['counts_err'],
                        ls=' ', marker='o', color=colors[j],
                        label=alldata[i][j].get_label(),
                        zorder=i)
            ax.plot(alldata[i][j].data['energy'],
                    alldata[i][j].data['model'],
                    ls='-', lw=3, color=colors[j],
                    zorder=len(alldata[i])+i, 
                    label='_nolegend_')
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.text(.95, .95, labels[i], transform=ax.transAxes, 
                fontsize=20, ha='right', va='top')
        ax.legend(loc=(.40, 0.05), fontsize=13, edgecolor='black')
        ax.set_xlim(right=10)
        fig.add_subplot(ax)

    #Plot residuals
    for i, ax in enumerate([axp1, axi1]):
        for j in range(len(alldata[i])):
            ax.errorbar(
                    alldata[i][j].residuals['energy'].astype(float), 
                    alldata[i][j].residuals['delchi'].astype(float),
                    xerr=alldata[i][j].residuals['energy_err'].astype(float), 
                    yerr=alldata[i][j].residuals['delchi_err'].astype(float),
                    ls=' ', marker='.', color=colors[j], alpha=0.8,
                    zorder=i)
        ax = plotparams(ax)
        ax.set_xscale('log')
        ax.set_xlim(right=10)
        fig.add_subplot(ax)


    plt.setp(axi0.get_xticklabels(), visible=False)
    plt.setp(axp0.get_xticklabels(), visible=False)

    fig.text(.03, .55, "Normalized Flux", ha='center', va='center', 
             rotation='vertical', fontsize=30)
    axi1.set_xlabel("Energy (keV)", fontsize=30)
    fig.savefig("plotunfolded.pdf", dpi=300)

if __name__ == '__main__':
    main()


    #Example use of multi_ufspec:
testufspec("PSR_B1821-24", [f"data_onpeak_{i}.txt" for i in range(4)], 
            [f"data_interpulse_{i}.txt" for i in range(4)], 
            [(0.97, 1.06), (0.98, 1.05), (0.99, 1.04), (0.0, 0.02)],
            [(0.50, .62), (.51, .61), (.52, .60), (.53, 59)])
