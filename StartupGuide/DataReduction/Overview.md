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

As mentioned above, the __`indir`__ input is the NICER observation ID (which is a directory). This directory contains the relevant event files and filter file, so the paths should almost always be kept as the defaults for the __`ufdir`__, __`cldir`__, __`ufafile`__, __`clfile`__, and __`mkfile`__ arguments.

The __`mpulist`__ arguments gives the MPU selections. The NICER detectors are grouped into 7 MPUs (decribed on the [Mission Guide](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/). The default is to use all the MPUs.


To apply the calibrations in nicercal, the NICER calibration database, stored in CALDB, should be used. The __`picalfile`__, __`pifastcalfile`__, and __`timebiascalfile`__ arguments use CALDB as the default value. 

Corrections based on leap seconds can be modified with the __`leapint`__ argument. This corection is done automatically by default. 

In addition to X-ray photon events, the detectors also record overshoot events, undershoot events, and forced triggers. The __`nicercal_filtexpr`__ argument extracts only X-ray events and forced triggers by default. While the event flag can be changed during the call to nicerl2, it is often simpler to do so in a latter call to niextract-events (discussed below). Each NICER event is given a 8-bit event flag, as described in the [Mission Guide](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/). The default expression is 'bxxxx00'. 

The __`calstation`__ argument simply says where the observation was performed. This is set to "FLIGHT" by default for in-flight observations. 

The use of niprefilter (not niprefilter2) can be controlled with the __`niprefilter`__ flag. This is set to "YES" by default. 

The next two arguments, __`orbfile`__ and __`attfile`__ also give paths to the auxiliary files. In the [ObsID](https://github.com/dmrowan/LommenResearchGroup/tree/master/StartupGuide/ObsID) section of the startup guide we discuss how to view these files in python. 

When the niprefilter (not niprefilter2) task is called, specific columns can be computed. This is controlled using the __`prefilter_columns`__ in nicerl2. By default, all relevant columns are computed. 

The South Atlantic Anomaly is a region where satellites experience high radiation levels. We typically want to filter out events detected while the ISS is flying over the SAA. This is done with the NICER region file, defined in CALDB. This is set using __`saaregfile="CALDB"`__. 

Similarly, information on ISS mission docking, ISS manuevering, and ISS robotics are stored in CALDB. The __`vehiclefile`__ argument can be used to disable retrieval of docking information or select on specific spacecraft. The __`issmanfile`__ and __`robofile`__ can be used to disable retrieval of maneuver and robotics information. 

Filtering of the South Atlantic Anomaly is turned on by default. The filtering of the NICER defined contour for the SAA is controled with __`nicersaafilt=YES`__. There is also a more general SAA region. Since the default is to use the NICER contour, __`saafilt`__ is set to zero. 

There are three flags corresponding to the instrument attitude that can be used to identify times of good tracking. These are added to the MKF file during the niprefilter call. `ATT_MODE` gives the pointing contorl mode. It is set to 1 for science mode. ` ATT_SUBMODE_AZ` gives the azimuthal tracking flag. If the flag is 2, tracking is good. Similarly, `ATT_SUBMODE_EL==2` indicates good tracking in elevation. These criteria are combined in a boolean expression to indicate tracking. `(ATT_MODE==1 && ATT_SUBMODE_AZ==2 && ATT_SUBMODE_EL==2)`. In nicerl2, using the default  __`trackfilt=YES`__ selects times of good tracking where the boolean expression is True. While we will usually want to keep this flag set to YES, for some science purposes, such as horizon crossings, we might want to allow times of bad tracking. 

The NICER startracker is used to achieve good pointing. A valid star tracker solution is saved in the flag `st_valid`. nicerl2 filters for times where __`ST_VALID=YES`__ by default. As was the case for trackfilt, this condition is typically left at its default unless performing observations like a horizon crossing. 

The angular distance gives the minimum offset between the pointing and the target. The default is __`ang_dist=0.015`__. This can be modified to change acceptable pointing offsets. 

The ELV parameter gives the elevation angle from the Earth limb to the target. The minimum value set in nicerl2 is __`ELV=15`__. This can be changed to a more or less restrictive criteria.

When the Earth is illuminated, we typically want a larger minimum value for the elevation. This is called the bright earth angle. The default value in nicerl2 is __`br_earth=30`__. 

The cutoff rigidity can be used to filter out cosmic ray events. nicerl2 gives the option to choose a range of cutoff rigidities to accept using the __`cor_range`__ parameter. There is no cutoff rigidity range defined defaulty. 

Undershoot events occur when the detector experiences optical loading or dark current. The maximum underhoot count range per MPU is typically set to 200. The user can specify the undershoot range rate using __`underonly_range`__

Overshoot events occur when high-energy particles saturate the detector. These events are due to the particle background, and should be filtered out. This is done using an expression that incorporates the cutoff rigidity `COR_SAX`. The default expression is
that the overshoot rate is less than 1.52\*COR_SAX\*\*(-0.633). This expression is controlled with the __`overonly_expr`__ argument. Liek the undershoot event rate, the overshoot event rate is also subject to an upper and lower limit. This is set to 0 to 1 by default with the __`overonly_range`__ argument.

We can restrict our events to times where a minimum number FPMs were on. By default, the __`min_fpm=7`__. 

Not all possile filtering options are given as arguments of nicerl2. Using the __`gtifiles`__ argument, additional GTI filters can be applied. 

The range of energies kept by nicerl2 is given with the __`pirange`__ argument. The default 0s 0.2 keV to 15.0 keV. This is written in the 10eV units of PI, so the default range is written as 20:1500. If you only wish to have energies up to 6keV, for example, the argument would be given as `pirange="20:600"`.

The trumpet filter is used to reject background events using the PI_RATIO calculation from nicercal. This is turned on by default. To turn the trumpet filter off, nicerl2 can be called with the argument __`trumpetfilt=NO`__. 

Not all nicerclean arguments are given as options for nicerl2. Additional arguments one would typically give nicerclean can instead be passed through nicerl2 using __`nicerclean_args`__. 

During the various stages of nicerl2, various temporary files are generated. The unfiltered events files are kept defaulty, but can be removed using __`cleanup_ufa_files`__. The other temporary files (like the GTIs used for the event extraction) are removed defaultly. This can be changed with the __`cleanup`__ flag. 

If nicerl2 has already been called on an observation, any subsequent calls will not overwrite unless the __`clobber`__ argument is set to yes. 

In case the task is producing unexpected results, the __`chatter`__ argument can control what information is printed. This is given on a range of zero to 5. The standard value of 2 gives typically logging of what step in the process is being performed. chatter=5 prints full debugging output. 

The nicerl2 parameters are saved into the output heaader if __`history`__ is set to YES. 



# nimaketime

# niextract-events

# fltime


### Notes:
\* on the pipeline account, ni_data_download has been modified to download selected obsIDs, but this change does not exist in the main NICERsoft branch. 

