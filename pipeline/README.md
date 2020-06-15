# NICER Pipeline Guide

The code in this directory is used for downloading & processing NICER observations for a handful of different sources. It is currently used for the following sources:
* PSR_B1937+21
* PSR_B1821-24
* PSR_J0218+4232
* PSR_B0531+21
* BKGD_RXTE (all 7 regions)
* V4641_Sgr

For each pipeline script the procedure generally follows these steps:
1. Download, unarchive, and decrypt all observations for source
2. Run nicerl2
3. Add geomagnetic index information
4. Perform additional filtering with nimaketime
5. Merge evt (and mkf, for BKGD)
6. Perform count rate cut (for pulsars)

The pipeline also automatically manages backups of evt and mkf files.

All the pipeline code is currently designed for multiprocessing. This means that the command line output of nicerl2, for example, looks like garbage because there are multiple sources at once. Someone might wanna add a multiprocessing flag in the future... 

__If the pipeline broke__ ... yeah sorry about that ... check out the section at the bottom of the file for things that have broken in the past and ideas of where to start. 

## Directory Structure

The pipeline is designed to work with the directory structure used on the Haverford cluster. 

In the /students/pipeline directory, there are separate directories for different versions of heasoft. We keep multiple versions to help track down problems that arise after version updates. At the time of writing, the current HEAsoft version is 6.27.2. Whenever heasoft 6.28 releases, a new directory should be created and the pipeline re-run. 

In the heasoft6.27 directory, there is a separate directory for each source. The sources are named to match the names given on the NASA site. For the background sources, there are separate subdirectories (e.g. BKGD_RXTE/BKGD_RXTE_1). 

When the pulsar pipe is run on PSR_B1821-24, it must be called in the PSR_B1821-24 directory. The pipeline code will catch most path errors, but still, best be careful. 

Within each source directory, a 'tmp' subdirectory is required. This is used for managing heasoft environmental variables in batch processing (discussed later). We also need a 'evt_backups' to store the backup files/log. For BKGD_RXTE sources, an additional 'mkf_backups' directory is required. If these directories don't exist, the pipeline will prompt the user before creating them. 

## Pulsar parameter files

We will need parameter files to add phases to the pulsar events. These are kept in the /students/pipeline/parfiles directory. Note that the crab pulsar requires multiple par files, and thus is a subdirectory. 

The path to the par file is given as an input to the pulsar pipeline.

## Pulsar Pipeline

The pulsar pipe is found in the pulsar_pipe.py file. It takes advantage of psrpipe.py in step #4 of the steps listed above. Therefore, environmental variables must be set as described in the NICERsoft readme.md. 

The primary command line invocation of the script looks something like this:
```
$ pulsar_pipe.py --update PSR_B1937+21 --par <parfile> --user <username> --passwd <passwd> -k
```
The `--update` argument calls the full pulsar pipeline for the source (B1937+21, in this case). We supply the par file with `--par`. 

We need to give the username and password for the NASA archive with the `--user` and `--passwd` arguments. The `-k` argument is used for the decryptkey. This will pullup a password prompt for the NICER decryption key. 

An alternative way to enter the key is by storing it in a file. The key should just be on the first line of the file. Then, the prompt can be bypassed using the argument `--k_file`:
```
$ pulsar_pipe.py --update PSR_B1937+21 --par <parfile> --k_file <keyfile>
```

Additional arguments are available:

`--emin`: Minimum energy to include in filtering. Default is 0.25 keV.

`--emax`: Maximum energy to include in filtering. Default is 12.0 keV.

`--cut`: Count rate cut in cts/sec. Default is 2.0 cts/s.

`--filterbinsize`: bin size for count rate cut. Default is 16s. 

`--trumpet`: Use the trumpet cut. Default is True.

`--keith`: Use the standard Keith Gendreau pipeline filters. Default is True


If you want to call the pipeline in python, it can be imported and run as:
```
>>> from pipeline import pulsar_pipe
>>> pulsar_pipe.update('PSR_B1937+21', heasarc_user, heasarc_passwd, 
                       decrypt_key, <parfile>, emin=0.25, emax=12.0,
                       cut=2.0, filterbinsize=16, trumpet=True,
                       keith=True)
```
(where heasarc_user and heasarc_passwd are the username and password for the data download. Decrypt_key is the decryption key string.

However you call the pipeline, after running a merged event file is created. If the sourcename is PSR_B1937+21, it will be 'PSR_B1937+21_combined_cut.evt'. It will automatically be entered into the backup directory and logged. 

The pipeline can also run each step separately if necessary:
1. run the pipeline on a single observation
2. download all the data, but don't do any filtering
3. merge all the event files in the current directory
4. create a backup

To run the pipeline on a single observation, specify the observation ID number with the obsID flag:
```
pulsar_pipe.py --obsID 3070020518 --par <parfile>
```
Since a count rate cut is normally done on a merged event file, this will prompt the user before performing to perform the cut. If it is selected, a new event file, '3070020518_pipe/cleanfilt_cut.evt', in this case, will be created. 

The clobber argument can be included to re-run the pipeline on an observation:
```
pulsar_pipe.py --obsID 3070020518 --par <parfile> --clobber
```

To download all the files for a given source, without doing the filtering, use the download argument:
```
$ pulsar_pipe.py --download PSR_B1821-24 --user <username> --passwd <passwd> -k
```
By default, existing observations are not overwritten. This can be changed using the clobber argument:
```
$ pulsar_pipe.py --download PSR_B1821-24 --user <username> --passwd <passwd> -k --clobber
```

The data download can also be called in Python using pipeline_utils.run_datadownload. This is a wrapper of custom_data_download.py, which in turn is a modification to NICERsoft ni_data_download.py. It's confusing, I know. Here's an example of calling a data download in python:
```
>>> from pipeline import pipeline_utils
>>> pipeline_utils.run_datadownload('PSR_B1821-24', <username>, <passwd>, './', <decryptkey>, clobber=False)
```
Note that the output directory for the downloaded observations can be changed by modifiying the outdir argument, which is set to './' in the example above. 

To merge all the event files currently piped, we can use the `--merge`. The most simple call would be:
```
$ pulsar_pipe.py --merge PSR_B1821-24
```
We can also choose to only merge events from observations up to a certain date. Say, for example we want our paper to only include observations taken before July 2019, we would call
```
$ pulsar_pipe.py --merge PSR_B1821-24 --max_date 2019-07-01
```
Note that the date must be formatted %Y-%m-%d. We can specify the output file with `--output`:
```
$ pulsar_pipe.py --merge PSR_B1821-24 --max_date 2019-07-01 --output merged_events.evt
```
The count rate cut is performed by default, so the --cut and --filterbinsize arguments can also be used with the merge argument. 

To call a merge in python, we can use the `run_niextract` function:
```
>>> import datetime
>>> from pipeline import pulsar_pipe
>>> dt = datetime.datetime.strptime('2019-07-01', '%Y-%m-%d')
>>> pulsar_pipe.run_niextract('PSR_B1821-24', max_date=dt, cut=2.0, filterbinsize=16.0, output='merged_events.evt')
```
Note that when the output file is given as an argument, the actual output will be merged_events_cut.evt because of how the count rate cut changes filenames. Someone should probably fix that, but eh, maybe later. 

Finally, we can make a backup of our event file using the `--backup` flag. We also have the option to enter a message with `--message`
```
$ pulsar_pipe.py --backup merged_events_cut.evt --message "I made this backup because I can"
```
Of course, we can also do this in python
```
>>> from pipeline import pipeline_utils
>>> pipeline_utils.product_backup('merged_events_cut.evt', message='remember to take out the trash before the next backup')
```

Again, even though we've shown how each of these pipeline functionalities can be called independently, either from the command line or a python session, they are all done together using the `--update` flag described at the top of the section.

## Background Pipeline

## V4641 Pipeline

## Cron job setup

## Backup system


## Uhoh it broke -- start here
So, the pipeline broke. Don't panic, [click here](https://www.youtube.com/watch?v=dQw4w9WgXcQ) for a guide on how to fix it. 
