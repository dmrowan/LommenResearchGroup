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
There is another argument `--silent_curl`, that removes the progress bar from the download. This is useful for when the output is redirected to a file, like we do in the crontab (see that section towards the bottom). 


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

Our pipeline code allows us to easily update a dataset when new data is available. If we want to automatically update a dataset, we can use the linux utility `cron` to schedule updates at given frequencies

We choose what scripts we want to run, and when to run them, using the crontab file. To edit, type `crontab -e`. The basic format of a cronjob is:

![cron_format](https://www.ostechnix.com/wp-content/uploads/2018/05/cron-job-format-1.png)

So, if we wanted to run 'myscript.py' every day at 1am, we would add the following line:
```
0 1 * * * <path_to_script>/myscript.py
```
(this assumes that the she-bang is at the top of the python script)

[Crontab.guru](https://crontab.guru/) is a great tool to help with the time information. 

Unfortunately, we can't add our pipeline scripts to the crontab the same way. Our pipeline scripts take advantage of HEAsoft, CALDB, PINT. In our bashrc, we have the following lines to specify all that info:
```
export HEADAS=/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23
alias heainit=". $HEADAS/headas-init.sh"
heainit

CALDB=/packages/heasoft-6.27.2/caldb; export CALDB
CALDBCONFIG=$CALDB/software/tools/caldb.config; export CALDBCONFIG
CALDBALIAS=$CALDB/software/tools/alias_config.fits; export CALDBALIAS
export TEMPO2=/packages/tempo2/T2runtime
```
A cron job doesn't use the environmental variables of the user that defined it. We must specify these again in the crontab file. There are a handful of ways to do this. What worked best for me was to be as explicit as possible and re-define the environmental variables. At the top of the file we can define the environmental variables:
```
HEADAS=/packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23
CALDB=/packages/heasoft-6.27.2/caldb
CALDBCONFIG=/packages/heasoft-6.27.2/caldb/software/tools/caldb.config
CALDBALIAS=/packages/heasoft-6.27.2/caldb/software/tools/alias_config.fits
TEMPO2=/packages/tempo2/T2runtime
```
I had some trouble getting `export` to work, so I just did everything explicitly here.

Our pipeline code also takes advantage of other scripts in the LommenResearchGroup and nicersoft repositories. We need to add that relevant information to the PATH and PYTHON_PATH. We also need to define LD_LIBRARY_PATH to define the library information. I found it was easier just to print each on the command line (e.g. `echo $PATH`) and copy them into the crontab. 

We also need to define two new environmental variables to make HEAsoft work with non-interactive shells. Add the following two lines to the crontab file:
```
HEADASNOQUERY=""
HEADASPROMPT=/dev/null
```
We should also specify that we are using bash:
```
SHELL=\bin\bash
```

At this point we've defined all the environmental variables to be the same as our bashrc. You might have noticed that we're missing one thing from the bashrc snippet. We use an alias `heainit` to intialize all the HEASoft environmental variables. We can't use the alias defined in the bashrc, so we have to remember to do this explicitly when we call our pipeline script. 

Let's take the pulsar PSR_B1821-24 as an example. We want to do a dataset update every Monday at midnight. In cron language, this is `0 0 * * 1`. 

We need to do a few things when this job is called:
1. Initialize heasoft (like we do in the bashrc)
2. cd to the PSR_B1821-24 directory
3. call pulsar_pipeline.py
4. (optional) write a log that the update was completed

We can either do all these steps explicitly in the crontab file (joining them with &&), or write a shell script and just call that. I think the second method is a bit easier. In the LommenResearchGroup/pipeline/cron_scripts directory there is a short script 'update1821':
```
#!/bin/bash

. /packages/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.23/headas-init.sh

cd /students/pipeline/heasoft6.27/PSR_B1821-24

/packages/python3.6.8/bin/python /homes/pipeline/scripts/LommenResearchGroup/pipeline/pulsar_pipe.py --update PSR_B1821-24 --user nicer_team --passwd sextant --k_file ~/decrypt_key --par /students/pipeline/parfiles/J1824-2452-frozen-dmx.par.nicerrefit --silent_curl

echo -n 'Completed PSR_B1821-24 ' >> ~/cron_log.txt
date >> ~/cron_log.txt
```

This is equivalent to the four steps listed above. A few notes: 
* All the paths are full absolute paths (not relative)
* We explicitly state which python we want to call for calling pulsar_pipe.py
* We use the --silent_curl argument to get rid of the curl progress bar. There would be a mess in /var/mail/pipeline if we don't

We can now add the following line in our crontab file:
```
0 0 * * 1 /homes/pipeline/scripts/LommenResearchGroup/pipeline/cron_scripts/update1821
```

### How to add another cron job
The PSR_B1821-24 example above can be replicated for any source. Here are the steps:
1. Identify the python command
2. Create a shell script that includes the four steps listed in the previous section
3. Define your cron time field
4. Add the shell script to the crontab file

How do we know if it works? On the Haverford cluster cron is setup so that every time a job is run, whether it fails or completes, it is entered into /var/mail/<user> file. Going all the way to the bottom of that will show information about the last job that ran. 
  
__Note that cron jobs are machine specific. If you create a job on dave.astro.haverford.edu it won't show up on frank.astro.haverford.edu, even for the same user__ (all the cron jobs in summer 2020 will be on Dave)

## Backup system


## Uhoh it broke -- start here
So, the pipeline broke. Don't panic, [click here](https://www.youtube.com/watch?v=dQw4w9WgXcQ) for a guide on how to fix it. 
