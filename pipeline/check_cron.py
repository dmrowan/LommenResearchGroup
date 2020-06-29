#!/usr/bin/env python

import argparse
from astropy import log
import numpy as np
import os
import pickle

#LommenResearchGroup imports
import niutils

# Dom Rowan 2020

desc="""
Parse the recent pipeline cron jobs to identify errors quickly

Read in file
Find lines corresponding to cron starts (get pipeline@dave and datetime)
If the top lines aren't the "head" of a cron job, get rid of them
    (we cant do anything about a partial log that happens if someone removes p      art of the mail file)

get list of lines for each cron job
search for errors
blacklist/whitelist errors

"""

class DenyList:
    
    def __init__(self):
        self.picklefile = os.path.join(os.environ.get('HOME'), 
                                       '.cron_denylist.pickle')
        if not os.path.isfile(self.picklefile):
            log.info("Creating cron denylist pickle")
            self.denylist = []
            self.write_pickle()
        else:
            self.read_pickle()

    def write_pickle(self):
        with open(self.picklefile, 'wb') as p:
            pickle.dump(self.denylist, p)

    def read_pickle(self):
        with open(self.picklefile, 'rb') as p:
            self.denylist = pickle.load(p)

    def add(self, key):
        if key not in self.denylist:
            self.denylist.append(key)
            self.write_pickle()

    def remove(self, key):
        if key in self.denylist:
            self.denylist.remove(key)
            self.write_pickle()

    def get(self):
        return self.denylist
        

def read_log(fname):
    if not os.path.isfile(fname):
        raise FileNotFoundError
        return -1

    with open(fname, 'r') as f:
        lines = f.readlines()

    criteria = ['From', 'Return-path:', 'Delivery-date:', 
                'To:', 'Subject:', 'X-Cron-Env:', 'Date:']

    first_parts = [ lines[i].split()[0] if (len(lines[i].split())>0)
                    else "" for i in range(len(lines)) ]

    if not all([ c in first_parts for c in criteria]):
        raise TypeError("Input file does not meet criteria for cron mail")
        return -1

    #Now want to find the lines where the cron jobs start
    start_idx = [ i for i in range(len(first_parts)-2) 
                  if first_parts[i:][:3] == ['From', 
                                             'Return-path:', 
                                             'Envelope-to:'] ] 

    job_groups = [ lines[start_idx[i]:start_idx[i+1]] 
                   if i != len(start_idx)-1
                   else lines[start_idx[i]:]
                   for i in range(len(start_idx)) ]

    return job_groups

def identify_job(job_lines):
    for l in job_lines:
        if l.startswith('Subject: Cron'):
            job = l.split('>',1)[-1].strip('\n')
            return job

def check_errors(lines):

    dl = DenyList()
    skip = dl.get()

    log.info("Skipping the following errors: {}".format(skip))

    bad_lines = []
    for l in lines:
        if any([ s in l.lower() for s in ['error', 'warning', 'traceback']]):
            if any([s in l.lower() for s in skip]):
                continue
            else:
                bad_lines.append(l)

    return bad_lines

def main(mail_file, jobs_to_check):

    if not niutils.check_iter(jobs_to_check):
        jobs_to_check = [jobs_to_check]

    #Get all cron log output 
    job_groups = np.array(read_log(mail_file))

    job_names = np.array([ identify_job(g) for g in job_groups ])
    
    idx_to_check = []
    #For each job we're interested in
    for job in jobs_to_check:

        #Find which cron lines match
        idx = [ i for i in range(len(job_names))
                if job in job_names[i] ]

        if len(idx) == 0:
            log.warning("No cron jobs found matching {0}".format(job))

        #only want to do error search on most recent
        idx_to_check.append(idx[-1])

    bad_lines = [ check_errors(job_groups[i]) for i in idx_to_check ]

    for i,j in zip(idx_to_check, jobs_to_check):
        log.info("Errors/Warnings for {}".format(j))
        bad_lines = check_errors(job_groups[i])
        for l in bad_lines:
            print(l.strip('\n'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("--job", 
                        help="job(s) to check most recent errors for",
                        nargs='+', type=str, default=None)
    parser.add_argument("--mail",
                        help="path to mail file",
                        type=str, default=None)

    parser.add_argument("--skip", type=str, nargs='+', default=None,
                        help="strings to skip in error output")
    parser.add_argument("--include", type=str, nargs='+', default=None,
                        help="strings to include in error output")

    args = parser.parse_args()

    if args.include is not None:
        dl = DenyList()
        for item in args.include:
            dl.remove(item)

    if args.skip is not None:
        dl = DenyList()
        for item in args.skip:
            dl.add(item)

    main(args.mail, args.job)
                    
