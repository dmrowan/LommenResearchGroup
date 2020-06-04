# Overview of data reduction procedures

NICER observations can be downloaded from the [NASA archive](https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html) using a curl command. The NICERsoft script [ni_data_download.py](https://github.com/paulray/NICERsoft/blob/master/scripts/ni_data_download.py) can be used for downloading observations in bulk*.

After downloading and unzipping files, the observation directory will be structured as described in the [ObsID section of the startup guide](https://github.com/dmrowan/LommenResearchGroup/tree/master/StartupGuide/ObsID). In the majority of cases, the next step of processing is to run `nicerl2`. This provides the standard level 2 processing procedures. 

After nicerl2 is run, cleaned event files can be found in the xti/even_cl directory of the observation. Depending on your source and the type of analysis one is trying to perform, additional filtering might be needed. For pulsars, this typically includes filtering on spaceweather criteria, hot detectors, and count rate. When performing additional filtering the three Heasoft/Ftools commands most typically used are `nimaketime`, `niextract-events`, and `fltime`. 

Here, we give an overview of these routines and discuss how they might be used for different analysis procedures. To see a bulk implementation of these techniques, see the [LommenResearchGroup pipeline procedures](https://github.com/dmrowan/LommenResearchGroup/tree/master/pipeline).

# nicerl2

[HEADAS Help File](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerl2.html5)

`nicerl2` combines multiple nicer calibration routines into one step to allow for easy handling of observations. In it's simplest form, nicerl2 can be run as
```
$ nicerl2 obsID
```
where obsID is the number corresponding to the observation directory. 

nicerl2 combines four steps:
1. [nicercal](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicercal.html)
2. [niprefilter2](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/niprefilter2.html)
3. [nimaketime](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nimaketime.html)
4. [nicermergeclean](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicermergeclean.html)

In the nicercal step, gain calibration and clock corrections are applied for each MPU separately. This calculates the PI columns.

In the niprefilter2 step, the filter file is transformed into a 'level 2' filter file. This means that additional columns are added to the mkf, notably X-ray rates, overshoot and undershoot rates, deadtime information, and Soyuz vehicle docking information. 

In the nimaketime step, the standard filtering criteria are used to generate good time intervals (GTI). We discuss these criteria in more detail below. 

Finally, nicermergeclean generates the ni\*_0mpu7_cl.evt file, which applies the GTI screening criteria to the (now calibrated) MPU events. This is also where the trumpet cut is made (through a call to [nicerclean](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerclean.html)). 

# What are the default filtering criteria? 

There are an abundance of optional arguments listed on the nicerl2 [help page](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerl2.html). Here, we briefly discuss each in the order they are listed under 'Parameters'. 

As mentioned above, the __`indir`__ input is the NICER observation ID (which is a directory). This directory contains the relevant event files and filter file, so the paths should almost always be kept as the defaults for the `ufdir`, `cldir`, `ufafile`, `clfile`, and `mkfile` arguments.

The `mpulist` arguments gives the MPU selections. The NICER detectors are grouped into 7 MPUs (decribed on the [Mission Guide](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/). The default is to use all the MPUs.


To apply the calibrations in nicercal, the NICER calibration database, stored in CALDB, should be used. The `picalfile`, `pifastcalfile`, and `timebiascalfile` arguments use CALDB as the default value. 

Corrections based on leap seconds can be modified with the `leapint` argument. This corection is done automatically by default. 

In addition to X-ray photon events, the detectors also record overshoot events, undershoot events, and forced triggers. The `nicercal_filtexpr` argument extracts only X-ray events and forced triggers by default. While the event flag can be changed during the call to nicerl2, it is often simpler to do so in a latter call to niextract-events (discussed below). Each NICER event is given a 8-bit event flag, as described in the [Mission Guide](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/). The default expression is 'bxxxx00'. 

The `calstation` argument simply says where the observation was performed. This is set to "FLIGHT" by default for in-flight observations. 

The use of niprefilter (not niprefilter2) can be controlled with the `niprefilter` flag. This is set to "YES" by default. 

The next two arguments, `orbfile` and `attfile` also give paths to the auxiliary files. In the [ObsID](https://github.com/dmrowan/LommenResearchGroup/tree/master/StartupGuide/ObsID) section of the startup guide we discuss how to view these files in python. 

When the niprefilter (not niprefilter2) task is called, specific columns can be computed. This is controlled using the `prefilter_columns` in nicerl2. By default, all relevant columns are computed. 



# nimaketime

# niextract-events

# fltime


### Notes:
\* on the pipeline account, ni_data_download has been modified to download selected obsIDs, but this change does not exist in the main NICERsoft branch. 

