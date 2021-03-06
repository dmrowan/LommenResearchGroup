# Lommen Research Group

A collection of procedures designed for analysis of NICER observations. Developed by students in Prof. Lommen's research group at Haverford College. 

<a href="https://gameon.nasa.gov/projects/deep-space-x-ray-navigation-and-communication/"><img src="https://gameon.nasa.gov/files/2017/03/nicer_logo.png" title="NICERlogo" alt="NICERlogo" width="400"></a>

## Getting Started
1. Ensure all prerequisites are met
2. Use git clone to clone the directory
```
git clone https://github.com/dmrowan/LommenResearchGroup
```
3. Modify PATH and PYTHONPATH in the bashrc file. If the directory you cloned LommenResarchGroup into is `<basedir>`, this would look like:
```
  export PATH=<basedir>/LommenResearchGroup:$PATH
  export PYTHONPATH=<basedir>/LommenResearchGroup/:$PYTHONPATH
```
Code that is in subdirectories can be used in multiple ways. Code can be imported in python:
```
from pipeline import pulsar_pipe
pulsar_pipe.allprocedures(*args, **kwargs)
```
For code that can be executed on the command line, the corresponding git subdirectory can be added to the .bashrc PATH variable. For example, to add LommenResearchGroup/pipeline to the PATH, include:
```
export PATH=<basedir>LommenResearchGroup/pipeline:$PATH
```


### Prerequisites
* Python 3
* PINT
* Heasoft
* [NICERsoft](https://github.com/paulray/nicersoft)
* See requirements.txt for python package requirements (dom still has to do this)

### PINT Path Requirements
We need the TEMPO2 variable to be defined. On the Haverford cluster, this can be added to the bashrc with:
```
export TEMPO2=/packages/tempo2/T2runtime
```

### HEASoft Path Requirements
Code that uses ftools, XSPEC, or other HEASoft commands require a few environmental variables. On the Haverford cluster, these can be added to the bashrc with:
```
export HEADAS=/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23
alias heainit=". $HEADAS/headas-init.sh"
heainit

CALDB=/packages/heasoft-6.27.2/caldb; export CALDB
CALDBCONFIG=$CALDB/software/tools/caldb.config; export CALDBCONFIG
CALDBALIAS=$CALDB/software/tools/alias_config.fits; export CALDBALIAS
```

## Project Directories
* bright_earth -- investigation into soft X-ray background dependence on bright earth angle
* pipeline -- collection of pipeline routines for various sources


## Authors

* **Dominick Rowan**
* **Liam Lynch**
* **Andrea Lommen**
* **Zaynab Ghazi**
* **Lauren Lugo**
* **Nate Ruhl**
* **Noah Schwab**
* **Mackenzie Tygh**
* **Sasha Levina**

See also the list of [contributors](https://github.com/dmrowan/LommenResearchGroup/contributors) who participated in this project.

