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

### Instrument Response
*Mention fparkey here

### Background Spectra

### Note about Exposure Keyword

## Plotting Data

## Ignoring Channels / Energy Ranges

## Fitting a Model

### Choosing a model

### Computing errors

### Plotting the model

## Plotting Spectra in Python

## Using Log Files
