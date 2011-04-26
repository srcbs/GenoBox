#!/usr/bin/python

import argparse
import logging
import subprocess
import pipelinemod
from ruffus import *

def read_libs(lib_file):
   '''Read library file and return dict'''
   
   from collections import defaultdict
   
   fh = open(lib_file, 'r')
   bam2lib = {}
   lib2bam = defaultdict(list)
   for line in fh:
      line = line.rstrip()
      fields = line.split('\t')
      bam2lib[fields[0]] = (fields[1], fields[2])
      
      lib2bam[fields[2]].append(fields[0])
      
   return (bam2lib, lib2bam)

parser = argparse.ArgumentParser(description='''
   Process BAM alignment files:
   Reads lib file that should contain (one line pr. bamfile)
   bam_file\tqual_threshold\tlib_name\n
   
   sort
   filter for mapping quality
   merge to libraries
   remove pcr duplicates
   merge to final bam
   ''')

# add the arguments
parser.add_argument('--lib', help='input file with libraries')
parser.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/')
parser.add_argument('--n', help='number of cpus for parallel runs [1]', default=1, type=int, choices=range(1, 17))
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--i Kleb-10-213361.bam --lib kleb.lib'.split())
#args = parser.parse_args('--lib NA12891.libs --n 5'.split())

# set logging
logger = logging.getLogger('genobox_bwa.py')
hdlr = logging.FileHandler('genobox_bwa.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set paths
paths = pipelinemod.setSystem()

# read library file
(bam2lib,lib2bam) = read_libs(args.lib)


# Sort bam #

task_params = []
for bam in bam2lib.keys():
   task_params.append([bam, bam+'.sort.bam'])

@files(task_params)
def sort_bam(in_bam, out_bams):
   '''Sort bam files'''
   
   paths = pipelinemod.setSystem()
   sam_cmd = paths['samtools_svn_home'] + 'samtools sort'
   
   out_bam_prefix = bam + '.sort'
   sam_arg = ' %s %s' % (bam, out_bam_prefix)
   sam_call = sam_cmd + sam_arg
   logger.info(sam_call)
   subprocess.check_call(sam_call, shell=True)


# Filter bam #

@follows(sort_bam)
@transform(sort_bam, suffix('.sort.bam'), '.sort.flt.bam', bam2lib)
def filter_bam(sort_bam, out_bam, bam2lib):
   '''Filter bam on qualities'''
   
   paths = pipelinemod.setSystem()
   
   in_bam = sort_bam.replace('.sort.bam', '')
   sam_cmd = paths['samtools_svn_home'] + 'samtools view'
   sam_arg = ' -b -q %s -o %s %s' % (bam2lib[in_bam][0], out_bam, sort_bam)
   sam_call = sam_cmd + sam_arg
   
   logger.info(sam_call)
   subprocess.check_call(sam_call, shell=True)


# Merge to libraries #

task_params = []
for lib in lib2bam.keys():
   sort_flt_bams = []
   for bam in lib2bam[lib]:
      sort_flt_bams.append(bam + '.sort.flt.bam')
   task_params.append([sort_flt_bams, lib+'.sort.flt.bam'])

@files(task_params)
@follows(filter_bam)
def merge_bam(list_bams, out_bam):
   '''Merge bam files'''
   
   if len(list_bams) == 1:
      mv_call = 'cp %s %s' % (' '.join(list_bams), out_bam)
      logger.info(mv_call)
      subprocess.check_call(mv_call,shell=True)
   else:
      paths = pipelinemod.setSystem()
      
      sam_cmd = paths['samtools_svn_home'] + 'samtools merge'
      sam_arg = ' %s %s' % (out_bam, ' '.join(list_bams))
      sam_call = sam_cmd+sam_arg
      
      logger.info(sam_call)
      subprocess.check_call(sam_call, shell=True)


# Remove duplicates #

@posttask(touch_file("alignment.complete.flag"))
@jobs_limit(5)
@follows(merge_bam)
@transform(merge_bam, suffix('.sort.flt.bam'), '.sort.flt.rmdup.bam', args.tmpdir)
def rmdup(in_bam, out_bam, tmpdir):
   '''Run rmdup on bam-file'''
   
   paths = pipelinemod.setSystem()
   
   java_call = paths['java_home']+'java -jar'
   
   metrics_file = in_bam + '.metrics.txt'
   picard_cmd = paths['picard_home'] + 'MarkDuplicates.jar'
   picard_arg = ' INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=%s VALIDATION_STRINGENCY=LENIENT' % (in_bam, out_bam, metrics_file, tmpdir)
   picard_call = picard_cmd + picard_arg
   
   call = '%s %s' % (java_call, picard_call)
   logger.info(call)
   subprocess.check_call(call, shell=True)


# Run pipeline #
pipeline_run([rmdup], multiprocess = args.n)
