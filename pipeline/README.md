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

We will need parameter files to add phases to the pulsar events. These are kept in the /students/pipeline/parfiles directory. Note that the crab pulsar requires multiple par files, and thus is a subdirectory. See the section on the Crab for more info on working with multiple par files and matching par files to observations.

The path to the par file is given as an input to the pulsar pipeline.

## Pulsar Pipeline

The pulsar pipe is found in the pulsar_pipe.py file. It takes advantage of psrpipe.py in step #4 of the steps listed above. Therefore, environmental variables must be set as described in the NICERsoft readme.md. 

For downloading & processing bulk data, the primary command line invocation of the script looks something like this:
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
Since a count rate cut is normally done on a merged event file, this will prompt the user before performing to perform the cut. If it is selected, a new event file, '3070020518_pipe/cleanfilt_cut.evt', in this case, will be created. The count rate cut and filterbin size can be changed with `--cut` and `--filterbinsize`, as mentioned above.

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
To download a single observation use `--obsID`. Note that this __does not__ process the observation, just downloads it. You would have to do a call with `--obsID`, like above, to process it after downloading. (working with single observations this way isn't really the intended use of pulsar_pipe, so the command line calls are a bit convoluted)

```
$ pulsar_pipe.py --download PSR_B0531+21 --obsID 1013010127 --user <username> --passwd <passwd> -k
```

There is another argument `--silent_curl`, that removes the progress bar from the download. This is useful for when the output is redirected to a file, like we do in the crontab (see that section towards the bottom). 


The data download can also be called in Python using pipeline_utils.run_datadownload. This is a wrapper of custom_data_download.py, which in turn is a modification to NICERsoft ni_data_download.py. It's confusing, I know. Here's an example of calling a data download in python:
```
>>> from pipeline import pipeline_utils
>>> pipeline_utils.run_datadownload('PSR_B1821-24', <username>, <passwd>, './', <decryptkey>, clobber=False)
```
Note that the output directory for the downloaded observations can be changed by modifiying the outdir argument, which is set to './' in the example above. You can specify which observations to download using the optional argument `obsIDs=['1013010121', '1013010127']`, for example. 

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

### Crab Pipeline

The Crab is run as part of the pulsar_pipeline but there are some differences in file handling worth mentioning. 

There are multiple parameter files generated from the M&M Crab paper table. These are stored in /students/pipeline/parfiles/crab. We have to match observations to par files. We do this with the `crab_utils.crab_par_match` function. This function identifies which observations have a par file and saves the info in a CSV.

```
from pipeline import crab_utils

par_dir='/students/pipeline/parfiles/crab'
username = <nasa site username>
passwd = <nasa site passwd>

crab_utils.crab_par_match(username, passwd, par_dir=par_dir)
```
This will write a csv and tex file (in the current directory) that shows which observations correspond to which par files. There is also an option to have some leway in the date range -- for example if you have an observation on May 15th but your par file only goes to May 10th, you can use `par_date_clearance=5` in `crab_par_match` to include it. By default this is turned off, so the clearance column in the CSV is all zeros. 

We can use this CSV to see what par files correspond to which obsID. For example, the Crab observation 1013010127 corresponds to 'march_2018_1.par'. We can call pulsar_pipe for a single obsID (that has already been downloaded): 

```
pulsar_pipe.py --obsID 1013010127 --par /students/pipeline/parfiles/march_2018_1.par --crab
```
Notice that there is an additional `--crab` argument that we didn't have above. The Crab observations have way more events so this argument changes the way multiprocessing is handled. For pulsars like PSR_B1937+21 pulsar_pipe is setup to run multiple observations at a time using multiple CPUs. When we use the `--crab` argument we use multiple CPUs for the photonphase processing step.

(in theory we could use --crab for any pulsar, it just wouldn't be useful unless the individual observations contained large (>100,000) events after filtering)


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


## Miscellaneous Notes

__Why did I get an error about a missing 'tmp' directory?__
When we call nicerl2 with multiprocessing, we need to be explicit about environmental variables used within HEASoft. If we don't do this, some processes won't be able to use the heasoft tools like nicerl2 correctly. The pipeline creates a path in the tmp directory for each observation. The PFILES variable is set before running nicerl2. To see explicitly where this is done in the python code, go to [LommenResearchGroup/pipeline/pipeline_utils.py](https://github.com/dmrowan/LommenResearchGroup/blob/master/pipeline/pipeline_utils.py) under the `run_nicerl2` function. To read more about why we need to do this, go to the [batch processing help page](https://heasarc.gsfc.nasa.gov/lheasoft/scripting.html)


## Uhoh it broke -- start here
So, the pipeline broke. Don't panic, [click here](https://www.youtube.com/watch?v=dQw4w9WgXcQ) for a guide on how to fix it. 

### How can I tell if the pipeline is working?
There are a few things you can do to see if the cron jobs are working:
* read the last entry in mail at /var/mail/pipeline. This shows the full output of the pipeline and is the best way to check for errors and warnings. It also shows the timestamp. See the 'check_cron.py' section below for more ways to use the mail file
* see if there are any recent changes in the data set directory. Example: `ls -ltr /students/pipeline/heasoft6.27/PSR_B1821-24`
* read the backup log for event files. This can be found within each dataset directory. Example: `vim /students/pipeline/heasoft6.27/PSR_B1821-24/evt_backups/log.txt`. The final step of the pipeline is creating a backup which will show up in the log with a time stamp. 

Even if the directory is updated and there is a new event file in the log, there might still be problems hiding. If something breaks new observations might not be downloaded/processed correctly. The best way to catch this is to use the pipeline.pipeline_utils.check_recent_obs function. Here's an example of using that function:

```
from pipeline import pipeline_utils

#This example is being called in the PSR_B1821-24 directory
# (/students/pipeline/heasoft6.27/PSR_B1821-24)

df_out = pipeline_utils.check_recent_obs('PSR_B1821-24', <username>, <passwd>, ncheck=5)
```
where <username> and <passwd> are the login credentials for the NASA target summary page. This will return a dataframe that says if the observation directory, pipe directory, event file exist. It also gives the length of the event file:
```
           obsdir pipedir evt  len(evt)
3070010331      ✓       ✓   ✓         0
3070010329      ✓       ✓   ✓      2720
3070010332      ✓       ✓   ✓      8685
3070010330      ✓       ✓   ✓      9890
3070010328      ✓       ✓   ✓      4619
```
(It will also print out this information)
  
 
So, if this showed that none of the recent observations had pipedir or event files, that would be cause for concern. It would also be worriesome if all the observations have zero length evts. 

#### Using check_cron.py
The output of the cron jobs is redirected to a mail file. On the Haverford cluster this is at /var/mail/pipeline. There is a lot of output from the various HEASoft and pipeline code and it can be difficult to identify errors. The check_cron.py script allows the user to search for errors/warnings on specific jobs quickly. To check for errors on the latest call of the 'update0218' routine, use this command:

```
$ check_cron.py --mail /var/mail pipeline --job update0218
```
The `--mail` argument gives the path to the mail file and the `--job` argument gives the job(s) you want to search for. Multiple jobs can be included in the search:
```
$ check_cron.py --mail /var/mail pipeline --job update0218 update1821
```

You might find that certain errors/warnings aren't super relevant. For example some nicersoft scripts have outdated xlim/ylim matplotlib args so there's a warning there. We don't necessarily care about this, so we can hide all warnings that have the word 'matplotlib' with `--skip`
```
$ check_cron.py --mail /var/mail/pipeline --job update0218 --skip "matplotlib"
```
If we also wanted to skip hot detector warnings, we'd say:
```
$ check_cron.py --mail /var/mail/pipeline --job update0218 --skip "matplotlib" "hot detectors"
```
It gets kinda annoying to have to say these every time, so __check_cron.py remembers what errors/warnings you skipped previously__. In other words, next time we call check_cron.py we won't see any matplotlib errors even if we don't use the --skip argument. 

If we had previously skipped an error/warning string but now want to see it, use the `--include` argument:
```
$ check_cron.py --mail /var/mail/pipeline --job update0218 --include "hot detectors"
```
Again, this will be remembered for future calls to check_cron. Everytime this script is called the errors on the 'denylist' are printed to the terminal so you'll know what you're missing. 


### EOFError: EOF when reading a line

### gzip unexpected end of file

### invalid columns

### curl response reading failed

### no event cl dir for obsID
