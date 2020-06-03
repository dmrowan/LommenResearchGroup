# Overview of data reduction procedures

NICER observations can be downloaded from the [Nasa archive](https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html) using a curl command. The NICERsoft script [ni_data_download.py](https://github.com/paulray/NICERsoft/blob/master/scripts/ni_data_download.py) can be used for downloading observations in bulk*.

After downloading and unzipping files, the observation directory will be structured as described in the [ObsID section of the startup guide](https://github.com/dmrowan/LommenResearchGroup/tree/master/StartupGuide/ObsID). In the majority of cases, the next step of processing is to run `nicerl2`. 



### Notes:
\* on the pipeline account, ni_data_download has been modified to download selected obsIDs, but this change does not exist in the main NICERsoft branch. 

