# Lommen Research Group

A collection of procedures designed for analysis of NICER observations. Developed by students in Prof. Lommen's research group at Haverford College. 

<a href="https://gameon.nasa.gov/projects/deep-space-x-ray-navigation-and-communication/"><img src="https://gameon.nasa.gov/files/2017/03/nicer_logo.png" title="NICERlogo" alt="NICERlogo" width="400"></a>

## Getting Started
1. Ensure all prerequisites are met
2. Use git clone to clone the directory
```
git clone https://github.com/dmrowan/LommenResearchGroup
```
3. Modify PATH and PYTHONPATH in the bashrc file. If the directory you cloned LommenResarchGroup into is <basedir>, this would look like:
```
  export PATH=<basedir>/LommenResearchGroup:$PATH
  export PYTHONPATH=<basedir>/LommenResearchGroup/:$PYTHONPATH
```
### Prerequisites
* Python 3
* PINT
* Heasoft
* [NICERsoft](https://github.com/paulray/nicersoft)
* See requirements.txt for python package requirements (dom still has to do this)

## Project Directories
* bright_earth -- investigation into soft X-ray background dependence on bright earth angle
* pipeline -- collection of pipeline routines for various sources


## Authors

* **Dominick Rowan**
* **Liam Lynch**
* **Zaynab Ghazi**
* **Lauren Lugo**

See also the list of [contributors](https://github.com/dmrowan/LommenResearchGroup/contributors) who participated in this project.

