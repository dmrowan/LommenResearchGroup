#!/bin/bash python

from bs4 import BeautifulSoup as bs
import argparse
from astropy import log
import os, sys
from subprocess import check_call, check_output
import xml
import datetime
import lxml
import fnmatch
import gnupg
import lxml.html as LH

# Liam's attempt at automatically scraping the NASA website every day to check for new uploaded data.
current_time = datetime.datetime.utcnow()
yesterday = current_time - datetime.timedelta(days=1)
log.info("current UTC time = {0}".format(current_time))
gpg = gnupg.GPG()
gpg.encoding = 'utf-8'
print(gpg)

desc=""" Scrapes the NASA website to find the observations requested either by target ID, processing date, or (eventually) OBSID. Downloads desired files in a new directory based on the download_path it receives. Currently the directory name is the date of the previous day, as the goal of the original implementation was to get all data uploaded in the last day from when the script was run. At some point i'll most likely change that, so that this script can be run actively with any of the below options, (currently dependent on the time argument). 


"""

parser=argparse.ArgumentParser(description=desc)
parser.add_argument("--target_id", help= "Only download files from a specific target", action= 'store', default= None)
parser.add_argument("--download_path", help= "Path to where you want to download new files, can be absolute or relative (example: '/student/ltlynch/NICER') Default is current working directory", default= os.getcwd())
parser.add_argument("--start_time", help= "Starting date you would like to look for files. Must be in UTC format (example: '2018-06-24T16:22:30') Default is 1 day ago.", action='store', default= yesterday.strftime('%Y-%m-%dT%H:%M:%S'))
parser.add_argument("--end_time", help= "End date for file search. Must be in UTC format. Default is current UTC time", default= current_time.strftime('%Y-%m-%dT%H:%M:%S'))
# Shouldn't have both a pipeline and par argument, temporary fix.
parser.add_argument("--pipeline", help= "Whether or not you want to run the pipeline (default= False)", default= False, action= 'store_true')
parser.add_argument("--par", help="Until I get something better to work you'll have to input an absolute path to a par file. Just so I can test the rest of the script.", action='store', default=None)
#parser.add_argument("--exposure_plot",help="Output figure file name (default=None)", default=None)

## Parse arguments
args = parser.parse_args()

# Let's you run shell commannds from this python script.

def runcmd(cmd):
	# CMD should be a list of strings since it is not processed by a shell
	log.info('CMD: '+" ".join(cmd))
	log.info(cmd)
	check_call(cmd,env=os.environ)

# Slightly altered the standard runcmd command in NICERsoft in order to let me store the output in a variable. Using check_output instead of check_call so that it actually waits for the function to finish running, and can be returned to a variable outside.
def returncmd(cmd):
	# CMD should be a list of stargs.rings args.since it is not processed by a shell
	log.info('CMD: '+" ".join(cmd))
	log.info(cmd)
	return check_output(cmd,env=os.environ)

# The main component of this script. Hopefully when it's done this will be able to scrape all of the NASA html data off of their website. As it's implemented now, I'll need to use wget in order to get passed the authentication form that pops up, and then store the html in a file name defined below. As of now scrape_HTML then goes through and seperates out the columns that I want: data "processing date" for the purpose of comparing it to current_time. Eventually it, or another seperate function, will make the connection between the time I found to be in my desired range, and it's associated file path, at which point i'll just have to curl that particular file. Loop through all files that were uploaded in the desired range. woo

def compare_times(date_dict):
     good_index_dict = {}
     UTC_convert = lambda x: datetime.datetime.strptime(x,'%Y-%m-%dT%H:%M:%S')
     for key in date_dict:
          UTC_datetime = UTC_convert(date_dict[key])
          if ((UTC_datetime.strftime("%s") < end_time.strftime("%s") and (UTC_datetime.strftime("%s")) > start_time.strftime("%s"))):
               #print("Index found in time interval: {0}, Date: {1}, Seconds: {2}".format(key, UTC_datetime, UTC_datetime.strftime("%s")))
               good_index_dict[key] = date_dict[key]
     return good_index_dict
     #print(map(lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'), date_list))

if args.start_time:
	start_time = datetime.datetime.strptime(args.start_time,'%Y-%m-%dT%H:%M:%S')
if args.end_time:
	end_time = datetime.datetime.strptime(args.end_time,'%Y-%m-%dT%H:%M:%S')
	log.info("download start_time: {0}".format(start_time))
	log.info("download end_time: {0}".format(end_time))

new_dir = start_time.strftime('%Y-%m-%d').replace("-", "_") + "-" + end_time.strftime('%Y-%m-%d').replace("-", "_")
log.info("Path where the data will be downloaded: {}".format(args.download_path))

#TODO Haven't really finished yet, so I'm gonna delete these passwords. If you want to use it you gotta put in the credentials for getting onto the website
cmd = ["wget", "--user", "Username here", "--password", "Password here", "--output-document", "NICER_data.html", "https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html"]
print(cmd)
runcmd(cmd)
with open('NICER_data.html', 'r') as f:	
	soup = bs(f, 'lxml') # Switching to lxml now, from html.parser
table = soup.body.tbody.children

# TODO All of this should eventually be changed to using pandas and csv. Right now it loops through rtd's and tr's and makes a ton of directories and it's awful. Rewrite eventually.

tr_list = {}
index_date_dict = {}
index_id_dict = {}
index_file_dict = {}
obsid_index_dict = {}
for tr_index, tr in enumerate(table):
	#print("Index: {}, Value: {}".format(tr_index,tr))
	for td_index, td in enumerate(tr):
		#if td_index == 0:
			#if td == '':
				#print("Index: {0} is empty".format(tr_index))
			#else:
				#print(type(td.encode('ascii')))
				#print(td.encode('ascii'))
				#print(str(td))
				#print(LH.fromstring(td.encode('ascii')).text_content())
		if td_index == 18:
			if td.string == None:
				log.info("Removing index: {}, because no processing date was found".format(tr_index))
				#index_date_dict.append[tr_index] = None
			else:
				#print(type(td))
				index_date_dict[tr_index] = str(unicode(td.string))
		elif td_index == 26:
			if not td.string == "|":
				#print("Index: {0}, td is: {1}".format(td_index,td.a.next_sibling.next_sibling["href"]))
				file_path = td.a.next_sibling.next_sibling["href"]
				index_file_dict[tr_index] = file_path
		elif td_index == 6:
			index_id_dict[tr_index] = str(unicode(td.string))
print("----------------------------")
print("CHECK: Index = 229, Date = {0}, Filepath = {1}".format(index_date_dict[229], index_file_dict[229]))
print("---------------------------------------------------------------------")
good_index_dict = compare_times(index_date_dict)
if args.target_id:
	target = args.target_id
	for key in good_index_dict.keys():
		if (not index_id_dict[key] == target):
			log.info("Removing index: {0}, because target: {1} was not specified".format(key, index_id_dict[key]))	
			del good_index_dict[key]
	print("-----------------------------------------------------------------")
	print("Indexes for on target observations uploaded in the last day: {}".format(good_index_dict))

# TODO Maybe this close statement should be made higher up, once soup variable is save, but troubleshoot when you move it.
f.close()
os.remove('NICER_data.html')

# This is really inconvenient, but I can't figure out how to get --create-dirs option to work in curl, so I'm going to just manually make a directory for these files. 
#print(download_path)
os.chdir(args.download_path)
# As of now, problematic, because if you run it once with one target_id and then again with another, you'd be using the same directory, but the par files would be different, psrpipe will then error.
if not os.path.isdir(new_dir):
	log.info("Creating directory '{0}' at path '{1}' ".format(new_dir, args.download_path))
	os.mkdir(new_dir)
os.chdir(new_dir)
#print(os.path.isdir(download_path))
for count, key in enumerate(good_index_dict):
	# Had to add this if statement check because NASA started uploading rows with all information filled in, but no file path.
	if key in index_file_dict.keys():
		cmd = ["wget", "--no-clobber", index_file_dict[key]]
		runcmd(cmd)
	else:
		log.info("Skipping download on index: {}, because no file was found on HEASARC".format(key))

print("------------------------------------------------------")
log.info("All downloading complete, data processing beginning")
#First go through and decrypt it all
log.info("Files have been downloaded, starting to process")
# One thing to note is that if you put a backslash in a cmd command it will error out. For using the find command, the semicolon doesn't need to be escaped.
print(os.getcwd())
log.info("Unzipping tar files")	

cmd = ["find", ".", "-name", '*tar.gz*', "-exec", "tar", "xvfp", "{}", ";", "-exec", "rm", "-f", "{}", ";"]	  
runcmd(cmd)

log.info("Decrypting gpg files")
# TODO Right now passphrase is just hardcoded in, should eventually switch to an argument like pasphrase-fd that allows the passphrase to be read from a file or something. More secure and don't have to have the passphrase be in the file if this ever ends up on github.	

# Also, if we ever upgrade to python 3, should find different implementation for this, (pathlib)?
for root, dirnames, filenames in os.walk(os.path.join(args.download_path, new_dir)):
  for filename in fnmatch.filter(filenames, '*.gpg'):
	  print(os.path.join(root,filename))
	  with open(os.path.join(root,filename),'r+') as f:
		  print(filename[:-4])
		  #TODO Same thing as before, until I do this a different way, i'm just deleting password
		  status = gpg.decrypt_file(f, passphrase="Decyprtion password here", output= os.path.join(root, filename[:-4]))
		  print("ok: {}".format(status.ok))
		  print("status: {}".format(status.status))
		  print("stderr: {}".format(status.stderr))
		  f.close()

log.info("Deleting gpg files")
cmd = ["find", ".", "-name", '*.gpg', "-exec", "rm", "-f", "{}", ";"]
runcmd(cmd)

log.info("Unzipping gz files")
cmd = ["find", ".", "-name", '*.gz', "-exec", "gunzip", "-f", "{}", ";"]
runcmd(cmd)

if args.par:
	log.info("Running the pipeline")
	for obsdir in os.listdir(os.getcwd()):
		if os.path.isdir(obsdir):
			print(obsdir)
			cmd = ["psrpipe.py", obsdir, "--mask", "-1", "--cormin", "5", "--filtpolar", "--par", args.par]
			runcmd(cmd)
# Could just implement the functionality for a list of files now, rather than looping through obsids, but will only do that if I have time. Wrote a quick little code that creates such lists elsewhere, just reference that, create the lists, then call psrpipe with @ (could also get a free merge option from that which would be nice)
	
