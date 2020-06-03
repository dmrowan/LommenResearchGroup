# Overview of data reduction procedures

NICER observations can be downloaded from the [Nasa archive](https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html) using a curl command. The NICERsoft script [ni_data_download.py](https://github.com/paulray/NICERsoft/blob/master/scripts/ni_data_download.py) can be used for downloading observations in bulk*.

After downloading and unzipping files, the observation directory will be structured as described in the [ObsID section of the startup guide](https://github.com/dmrowan/LommenResearchGroup/tree/master/StartupGuide/ObsID). In the majority of cases, the next step of processing is to run `nicerl2`. This provides the standard level 2 processing procedures. 

After nicerl2 is run, cleaned event files can be found in the xti/even_cl directory of the observation. Depending on your source and the type of analysis one is trying to perform, additional filtering might be needed. For pulsars, this typically includes filtering on spaceweather criteria, hot detectors, and count rate. When performing additional filtering the three Heasoft/Ftools commands most typically used are `nimaketime`, `niextract-events`, and `fltime`. 

Here, we give an overview of these routines and discuss how they might be used for different analysis procedures. To see a bulk implementation of these techniques, see the [LommenResearchGroup pipeline procedures](https://github.com/dmrowan/LommenResearchGroup/tree/master/pipeline).

# nicerl2

# nimaketime

# niextract-events

# fltime


### Notes:
\* on the pipeline account, ni_data_download has been modified to download selected obsIDs, but this change does not exist in the main NICERsoft branch. 

