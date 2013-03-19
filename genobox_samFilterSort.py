#!/panvol1/simon/bin/python2.7

import argparse
import subprocess
import logging
import genobox_modules

def samFilterSort(i, q, m, o, F):
   '''Filters bam on quality and sort'''
   
   paths = genobox_modules.setSystem()
   sam_cmd = paths['samtools_svn_home'] + 'samtools'
   if q == 0:
      if F:
         call = '%s view -u -F %i %s |  %s sort -m %s - %s' % (sam_cmd, F, i, sam_cmd, m, o)
      else:
         call = '%s sort -m %s %s %s' % (sam_cmd, m, i, o)
   else:
      if F:
         call = '%s view -u -F %i -q %i %s |  %s sort -m %s - %s' % (sam_cmd, F, q, i, sam_cmd, m, o)
      else:
         call = '%s view -u -q %i %s |  %s sort -m %s - %s' % (sam_cmd, q, i, sam_cmd, m, o)
   #logger.info(call)
   subprocess.check_call(call, shell=True)

def picardFilterSort(i, q, o):
   '''Filters bam on quality and sort using picard'''
   
   paths = genobox_modules.setSystem()
   call = '''%sjava -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms1500m -Xmx1500m -jar %s/ViewSam.jar INPUT=%s ALIGNMENT_STATUS=Aligned VALIDATION_STRINGENCY=LENIENT | perl -ane 'if ($_ =~ m/^@/) {print $_;} else {if ($F[4] >= %i) { print $_ }}' | %sjava -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms4500m -Xmx4500m -jar %s/SortSam.jar INPUT=/dev/stdin OUTPUT=%s SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/panvol1/simon/tmp MAX_RECORDS_IN_RAM=1000000''' % (paths['java_home'], paths['picard_home'], i, q, paths['java_home'], paths['picard_home'], o)
   subprocess.call(call, shell=True)
   


parser = argparse.ArgumentParser(description='''
   Filters bam for quality and sort
   ''')

# add the arguments
parser.add_argument('--i', help='input bam')
parser.add_argument('--q', help='quality cutoff [30]', default=30, type=int)
parser.add_argument('--F', help='bit-field to exclude [None]', default=None, type=int)
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
#samFilterSort(args.i, args.q, args.m, args.o, args.F)
picardFilterSort(args.i, args.q, args.o)