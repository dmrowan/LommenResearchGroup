# What is in a NICER observation directory?

If you head to the [NICER observations page](https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html) you'll find a table listing all the observations for each source. The leftmost column gives the observation ID (obsID). When you download the observation, it will be in the form of a directory with that file name. What's in that directory?

Here is the directory structure for obsID 1200040101:

![Image of DIR structure](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/figures/nicer_dirs.png)

Depending on if you've unzipped the files, you might not see the .gz extension.

The auxil subdirectory contains three four files:
* The attitude file, with extension .att, gives information about the telescopes slewing and position. It is in the form of a FITS file. It is used in niprefilter2
* The orbit file, with extension .orb, gives information about the orbital parameters at the observation time. 
* The filter file, with extension .mkf, is a FITS file containing over 90 different columns with information about the observation as a function of time. Columns include earth elevation, number of active FPMs, and cutoff rigidity
* The catalog file with extension .cat gives ...

The log directory contains two html files, an error log and a job log. I don't think the job log keeps track of calls to nicerl2, but it does tell you about what reduction procedures were applied to the data downloaded. 

The final directory, xti, contains detection information from the x-ray telescope instrument. There are three subdirectories
* event_cl gives the 'cleaned events' following the reduction procedures described in the log
    * ni*_0mpu7_cl.evt is the merged event file from all the MPUs 0-6
    * ni*_0mpu7_ufa.evt is the unfiltered event file from all the MPUs 0-6
* event_uf gives the unfiltered events for each MPU
* The houskeeping directory (hk) contains a FITS file for eac MPU. I think these give relavent parameters about the detectors and is used for niprefilter2

We can read each of these files with Python. See the Jupyter Notebook for an example. 
