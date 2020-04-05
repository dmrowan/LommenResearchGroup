#!/usr/bin/env python

from LCClass import LightCurve
import matplotlib.pyplot as plt
import niutils

def main():
    lc1821 = LightCurve("PSR_B1821-24/PSR_B1821-24_combined.evt")
    lc0218 = LightCurve("PSR_J0218+4232/PSR_J0218+4232_combined.evt")

    fig, ax = plt.subplots(2, 1, figsize=(8, 8))

    ax[0], _ = lc1821.fit_two('lorentzian', ax=ax[0], label=False, annotate=False)
    ax[1], _ = lc0218.fit_two('gaussian', ax=ax[1], label=False, annotate=False)

    ax[1].set_xlabel("Pulse Phase", fontsize=25)
    ax[0].text(.08, .95, r'PSR B1821$-$24', ha='left', va='top', 
               fontsize=20, transform=ax[0].transAxes,
               bbox=dict(facecolor='white', edgecolor='none', alpha=0.6))
    ax[1].text(.08, .95, r'PSR J0218$+$4232', ha='left', va='top', 
               fontsize=20, transform=ax[1].transAxes,
               bbox=dict(facecolor='white', edgecolor='none', alpha=0.6))

    ax[0].tick_params(labelbottom=False)
    #plt.setp(ax[0].get_yticklabels()[0], visible=False)
    
    fig.text(.04, .5, r'Photon Counts', ha='center', va='center',
             rotation='vertical', fontsize=25)

    plt.subplots_adjust(hspace=0, bottom=.08, top=.94, right=.98, left=.15)

    fig.savefig("poster_plot.svg")



if __name__ == '__main__':
    main()
