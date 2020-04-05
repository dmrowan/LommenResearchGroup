#!/usr/bin/env python

import argparse
import os
import numpy as np
import xspeclog
import tempfile
import subprocess

desc="""
Script to run an xspeclog file in an interactive XSPEC window
"""

def main(fname, default_params=True):
    logfile = xspeclog.logfile(fname)

    commands = logfile.get_commands()

    print(commands)
    commands = list(filter(None, commands))
    commands = [ c for c in commands if c[0] != '@' ]
    if 'exit' in commands[-1]:
        commands = commands[:-1]

    if default_params:
        for i in range(len(commands)):
            if commands[i].startswith('mo'):
                commands.insert(i+1, '/*\n')

    with open('tmpfile.xcm', 'w') as f:
        for c in commands:
            f.write(c)


    with open("tmpexpect", 'w') as f:
        f.write("#!/usr/bin/env expect\n")
        f.write("spawn xspec\n")
        f.write('expect "XSPEC12>"\n')
        f.write("send @tmpfile.xcm\n")
        f.write("interact\n")
        f.close()

    subprocess.run(['chmod', 'u+x', 'tmpexpect'])
    subprocess.run(['./tmpexpect'])



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fname", help='logfile', type=str)

    args = parser.parse_args()

    main(args.fname)
