#!/panvol1/simon/bin/python2.7

import argparse
import subprocess
import logging
import genobox_modules

def samFilterSort(i, q, m, o):
   '''Filters bam on quality and sort'''
   
   paths = genobox_modules.setSystem()
   sam_cmd = paths['samtools_svn_home'] + 'samtools'
   if q == 0:
      call = '%s sort -m %s %s %s' % (sam_cmd, m, i, o)
   else:
      call = '%s view -u -q %i %s |  %s sort -m %s - %s' % (sam_cmd, q, i, sam_cmd, m, o)
   logger.info(call)
   subprocess.check_call(call, shell=True)
   

parser = argparse.ArgumentParser(description='''
   Filters bam for quality and sort
   ''')

# add the arguments
parser.add_argument('--i', help='input bam')
parser.add_argument('--q', help='quality cutoff [30]', default=30, type=int)
parser.add_argument('--m', help='memory for sort', default=500000000, type=int)
parser.add_argument('--o', help='prefix output bam')
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--i test.bam --q 0 --m 5000000 --o test.sort '.split())

# set logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# filter
samFilterSort(args.i, args.q, args.m, args.o)
