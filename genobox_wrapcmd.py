#!/panvol1/simon/bin/python2.7

import argparse
import subprocess
import logging
import random
import string
import os
import sys

parser = argparse.ArgumentParser(prog='genobox_wrapcmd.py', description='''Take input as command and submit to msub (this way pipes and redirects can be done)''', usage='%(prog)s [cmd]')

parser.add_argument('cmd', help='command to submit')
parser.add_argument('--name', help='run name', default='run_genobox_wrapcmd')
parser.add_argument('--l', help='extension string', default='nodes=1:ppn=1,mem=1gb,walltime=86400')
parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('time'.split())

# start logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set home
home = os.getcwd()

# create pbsjob
N = 10
rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))
filename = 'pbsjob.tmp%s' % rand

fh = open(filename, 'w')
fh.write('#!/bin/sh\n\n')
fh.write('%s\n' % args.cmd)
fh.close()

# submit job to msub
call = 'msub -d %s -l %s -q %s -r y -N %s %s' % (home, args.l, args.queue, args.name, filename)

logger.info(call)
id = subprocess.check_output(call, shell=True)
id = id.split('\n')[1]

# remove pbsfile
genobox_modules.rm_files([filename])

# output job-id
sys.stdout.write('%s\n' % id)
