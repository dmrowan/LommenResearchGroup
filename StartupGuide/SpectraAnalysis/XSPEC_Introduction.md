# XSPEC Introduction

XSPEC is a command line software packaged used for spectral analysis of X-ray observations. It can be used with a variety of X-ray observatories, including NICER.

There is also a Python implementation, PyXSPEC. Here, we focus only on the command line program. We also include discussion on how to use XSPEC in Python with `pexpect`.

[Papers by Keith Arnaud on XSPEC development](https://ui.adsabs.harvard.edu/search/filter_author_facet_hier_fq_author=AND&filter_author_facet_hier_fq_author=author_facet_hier%3A%220%2FArnaud%2C%20K%22&fq=%7B!type%3Daqp%20v%3D%24fq_author%7D&fq_author=(author_facet_hier%3A%220%2FArnaud%2C%20K%22)&q=title%3A%22XSPEC%22&sort=date%20desc%2C%20bibcode%20desc&p_=0)

## Requirements
XSPEC is included in HEAsoft. Therefore, if the relevant HEAsoft information is included in the .bashrc, XSPEC should be able to be called from the command line. 

The CALDB database is also required. Again, this is part of the HEASARC software, and the CALDB, CALDBCONFIG, and CALDBALIAS environmental variables should be defined. 

Information on HEAsoft installation can be found [here](https://heasarc.gsfc.nasa.gov/lheasoft/install.html).

CALDB documentation can be found [here](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_doc.html). 

For those working on the Haverford Cluster, HEAsoft and CALDB are already installed. As of June 2020, the current HEAsoft version is 6.27.2. The following lines can be added to the .bashrc file:

```
export HEADAS=/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23
alias heainit=". $HEADAS/headas-init.sh"
heainit

CALDB=/packages/heasoft-6.27.2/caldb; export CALDB
CALDBCONFIG=$CALDB/software/tools/caldb.config; export CALDBCONFIG
CALDBALIAS=$CALDB/software/tools/alias_config.fits; export CALDBALIAS
```

XSPEC can be launched on the command line by simply typing `xspec`. This will bring up a command promt `XSPEC12>`. 

To exit XSPEC at any time, simply type `exit`.

## Example Data Set

In this guide we will use an emission spectra from millisecond pulsar PSR B1937+21. The spectrum file is 1937spectrum.pha. 

We will also need three other files:
* nixtiref20170601v001.rmf
* nixtiaveonaxis20170601v002.arf
* background.pha

The usage of these files is explained below in the 'Instrument Response' and 'Background Spectra' sections. 

## Loading spectra into XSPEC

Launch XSPEC on the command line by typing `xspec`. This will bring up version information and a prompt `XSPEC12>`
![xspec0](assets/xspec_0.png)

Before analyzing a spectra, it is typically good to create a log file. This will save all the commands and output run during the interactive XSPEC session to a file. Later in this guide, we discuss how this can be parsed with Python. To create a log file, type `log` and the name of the logfile.
```
XSPEC12> log 1937_logfile.txt
```
We are now ready to load in our pulsar emission spectrum. XSPEC allows for multiple data sets or observations to be modeled independently, so spectra are numbered when loaded in. For our purposes, this simply means that we are assigning our spectra to the first spectrum position. To load the data, type
```
XSPEC12> data 1 1937spectra.pha
```
The output will show some information about the spectrum:
![xspec2](assets/xspec_2.png)

XSPEC gave us a warning here -- we're missing a response file. Depending on how you've generated a spectrum, this warning may or may not show up. In the next section, we go over what a response file does and how to load it into XSPEC. 

### Instrument Response

Here is some text from Dom's thesis about the instrument response:

> When we observe an X-ray source, the recorded data is the number of photons detected in each instrument channel $n$ (Arnaud, 1999). This is related to the intrinsic spectrum of the source, ![f(E)](https://render.githubusercontent.com/render/math?math=f(E)), where ![E](https://render.githubusercontent.com/render/math?math=E) is energy, and the response function of the instrument, ![R(n, E)](https://render.githubusercontent.com/render/math?math=R(n%2C%20E)), by
> 
> ![C(n) = \int_0^{\infty} f(E) R(n,E)\ dE](https://render.githubusercontent.com/render/math?math=C(n)%20%3D%20%5Cint_0%5E%7B%5Cinfty%7D%20f(E)%20R(n%2CE)%5C%20dE)
> 
> In reality, we don't have a continuous function for ![R(n, E)](https://render.githubusercontent.com/render/math?math=R(n%2C%20E)), but rather a matrix that defines a discrete function defined as
> 
> ![R_D(n, m) = \frac{\int_{E_{m-1}}^{E_m} R(n,E)\ dE}{E_m - E_{m-1}}](https://render.githubusercontent.com/render/math?math=R_D(n%2C%20m)%20%3D%20%5Cfrac%7B%5Cint_%7BE_%7Bm-1%7D%7D%5E%7BE_m%7D%20R(n%2CE)%5C%20dE%7D%7BE_m%20-%20E_%7Bm-1%7D%7D)
> 
> where ![m](https://render.githubusercontent.com/render/math?math=m) is the index of the energy range. 

In other words, a response matrix allows us to relate the energy channel of the photon with the energy in units of keV. [Page, M. (2015)](https://arxiv.org/abs/1506.07015) has a good description of the relation. 

Where can we find a response matrix? The NICER instrument response can be found in the calibration database (part of the HEAsoft installation). It can be found on the Haverford cluster at /packages/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf. It is also included in this github directory. 

To apply the instrument response, we ue the command `resp`. We want to apply our response file to the first spectrum in our XSPEC session. Therefore we write:
``` 
XSPEC12> resp 1 nixtiref20170601v001.rmf
```
There is also an ancillary response file. This file contains information about the effective area of the detector. This is also found in CALDB, at /packages/caldb/data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf. Again, this file is also included in the github directory. We load the ARF in a way similar to the response file:
 ```
 XSPEC12> arf 1 nixtiaveonaxis20170601v002.arf
 ```
 If both of these commands were successful, we see
 
 ![xspec_3](assets/xspec_3.png)
 

## Plotting Data

XSPEC allows us to plot the spectra and in a variety of forms, both with and without a fit model. First, we have to choose a plot device. The following command will bring up a plotting window

```
XSPEC12> cpd /xs
```
The window that appears should look like this:
![window](assets/xspec_plot_0.png)

The first thing we can do is simply plot the data:
```
XSPEC12> plot data
```
This will update the plot window and should look something like this
![xspec_plot_1](assets/xspec_plot_1.png)

Even though we've added the instrument response and ARF, we still have the channel axis on the bottom. If we want to plot the energy, we can adjust the axis with `setplot`
```
XSPEC12> setplot energy
```
We can then refresh the plot
```
XSPEC12> plot
```
This should now show an energy axis. 
![xspec_plot_2](assets/xspec_plot_2.png)

Looks like somethings not right. We have a log scale on the x-axis and it appears to have struggled with the lowest energy bin. The strategy here will be to use ignore/notice to select energies channels we wish to include

## Ignoring Channels / Energy Ranges

In most cases, we are only interested in a range of energies that the detector collects. For NICER, we know that events with energies less than 0.2 and greater than 18keV are detected, yet this is not the optimal range for the calibrated instrument. 

We can choose what energies or channels to include in our spectral analysis using the `ignore` and `notice` commands. If we want to ignore the energies less than 1keV, we can write
```
XSPEC12> ignore **-1.0
```
Alternatively, if we were working with the channel axis (you can use `setplot channel` to get back to that plot), we could say something like
```
XSPEC12> ignore **-24
```
There is a very subtle difference in ignoring energy bins and channel bins. If we want to ignore energy ranges, we must include a decimal point. If not, it will be interpreted as a channel range. 

We can also ignore energies up to a certain value and above a certain value. For example:
```
XSPEC12> ignore **-1.0, 10.0-**
```
If the command is successful, XSPEC will show you how many channels were ignored:
![xspec_5](assets/xspec_5.png)

After refreshing the plot, we see that the x-axis range has changed. Let's say that we changed our mind and we actually want to include some energies less than 1.0 keV. We can use the notice command to include channels that were previously ignored
```
XSPEC12> notice 0.3-1.0
```
Again, the use of decimal points indicates the selection is in units of keV, rather than by channel.

Note that if multiple spectra are used, ignoring/noticing can be done separately for each spectrum. See the [ignore documentation](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node74.html) for more information. 

## Background spectrum

We will probably need a background spectrum to be able to distinguish the source from the background. Here, the background.pha file can be loaded in with
```
XSPEC12> back 1 background.pha
```
This analysis is an example of a difference spectra. Check out the other section of the StartupGuide for information on background spectrum generation.

## Fitting a Model

We are now ready to fit a model to our data. There are a ton of models available in XSPEC, and they can be combined additively and multiplicatively. To see a list of models (and a short description) go to the [XSPEC alphabetical summary page](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node130.html)

For PSR B1937+21, we want to use an absorbed power law model. Looking on the list of models, we see that there isn't a single model that suits our needs, so we will have to use a combination of the 'tbabs' and 'powerlaw' models.

Before using the tbabs models, we have to [set the abundances](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node117.html)
```
XSPEC12> abund wilm
```
Now we are ready to define the model
```
XSPEC12> model tbabs*powerlaw
```
This will pull up a slightly different prompt. In this case it should look like `1:TBabs:nH>`. This is where you can enter initial values before fitting. You also have the option to leave them at default by hitting enter. For now, just hit enter
 three times. The output in XSPEC should look something like this:
 ![xspec_7](assets/xspec_7.png)
 
 As we can see, the fit is really terrible using the default values. Type `renorm` to fix the normalization. Then fit the model with
 ```
 XSPEC12> fit <iterations>
 ```
 The computation time is pretty short, so I normally go for a larger number of iterations (~100,000). You should see the fit statistics printed again after the fitting is complete. 

### Computing errors

The best fit parameters will have errors listed next to them. __These errors should not be used for reporting the significance__. Instead, use the `error` command:
```
XSPEC12> error 1
```
This prints the error on the first parameter. You can change the desired confidence level and also select multiple parameters. Additional examples are given on the bottom of the [documentation page](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node80.html)

### Plotting the model
Now that we've fit a model we can see how it compares against the data. There are multiple axis we can plot, but I usually like to do
```
XSPEC12> plot ufspec delchi
```
This will make two panels, with the top showing the unfolded spectrum and the bottom showing the residuals. 

## Plotting Spectra in Python
The plot looks great, but it's probably not as nice as your other matplotlib plots. To export the data shown in the current plot, use the following commands
```
XSPEC12> ipl
PLT> wdata mydata.txt
PLT> exit
```
This will save the current plot information for all panels. It might take a bit of investigating to figure out what the columns are. There's a class `xspecdata` in `xspeclog.py` that is able to parse the pulsar data from ufspec delchi. Mileage may vary for other sources. 

## Using Log Files

The log file is also super useful for keeping track of many spectra. We might generate many spectra with different data selections (see multispectra.py for an example) and want to save the results in a table. We can parse the log file as a text file to find the fit parameter values, errors, and chi2. See `xspeclog.py` for an example implementation. 

## Scripting in XSPEC

Doing XSPEC many times really sucks. Luckily, it's not too hard to script. Write a file with extension xcm that lists all the commands:
```
data 1:1 onpeak_1.pha 2:2 nu30101031002_on_off_bk_sr.pha.grp 3:3 b1937_xmm_mos_src.pha.grp
resp 1:1 /packages/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf
arf 1:1 /packages/caldb/data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf
back 1:1 offpeak.pha

cpd /xs
abund wilm

ignore 1:**-1., 9.-**
ignore 2:80.-**
ignore 3:**-1.,10.-** 

model powerlaw*tbabs
```
If we name the file 'xspecscript.xcm' (and make it executable), this can be called in XSPEC
```
XSPEC12>@xspecscript.xcm
```

## XSPEC in Python

There's a python module for XSPEC. The documentation can be found [here](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/PyXspec.pdf). I don't have much experience using this, but I know Lauren had some success with it. 

You can also call command line XSPEC in python using pexpect. Here's an example of how I used it in `genspectra.py`:
```
import pexpect

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
```
(this method could also be used to script other command line ftools like XSELECT).

To see a more involved example where we use this in cojunction with xcm scripts, check out `multispectra.py`. 
