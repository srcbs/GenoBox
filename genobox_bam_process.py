#!/usr/bin/python

import argparse
import logging
import subprocess
import pipelinemod
import os

# read library file 
def read_libs(lib_file):
   '''Read library file and return dict.
   Reads lib file should contain (one line pr. bamfile):
   bam_file\tqual_threshold\tlib_name\n
   '''
   
   from collections import defaultdict
   
   fh = open(lib_file, 'r')
   bam2lib = {}
   lib2bam = defaultdict(list)
   for line in fh:
      if line.startswith('#'):
         continue
      
      line = line.rstrip()
      fields = line.split('\t')
      bam2lib[fields[0]] = (fields[1], fields[2])
      
      lib2bam[fields[2]].append(fields[0])
   return (bam2lib, lib2bam)

# filter bam and sort
def bam_filter_sort(lib2bam, bam2lib):
   '''Filter bam on quality and sort on stream'''
   
   paths = pipelinemod.setSystem()
   cmd = 'python2.7 ' + paths['pyscripts_home'] + 'samFilterSort.py'
   calls = []
   
   # set infiles and outfiles
   infiles = lib2bam.values()
   
   calls = []
   outfiles = []
   for lib_files in infiles:
      for f in lib_files:
         sort_prefix = f + '.flt.sort'
         out_bam = f + '.flt.sort.bam'
         outfiles.append(out_bam)
         arg = ' --i %s --q %s --o %s' % (f, bam2lib[f][0], sort_prefix)
         calls.append(cmd+arg)
   return (calls, outfiles)

# merge to libraries
def merge_bam(libs, lib_infiles, add_suffix=False, final_suffix=''):
   '''Merge bam files to libraries'''
   
   calls = []
   outfiles = []
   for i in range(len(libs)):
      lib = libs[i]
      
      # set input and output files
      # add suffix to files (this is if they are given as original filenames, before filter+sort)
      if add_suffix:
         list_bams = []
         for infile in lib_infiles[i]:
            list_bams.append(infile + '.flt.sort.bam')
      else:
         list_bams = lib_infiles[i]
      
      # add suffix to outfile if set else do not
      out_bam = lib + final_suffix
      outfiles.append(out_bam)
      
      if len(list_bams) == 1:
         call = 'cp %s %s' % (' '.join(list_bams), out_bam)
      else:
         sam_cmd = paths['samtools_svn_home'] + 'samtools merge'
         sam_arg = ' %s %s' % (out_bam, ' '.join(list_bams))
         call = sam_cmd+sam_arg
      calls.append(call)
   return (calls, outfiles)

# rmdup
def rmdup(infiles, tmpdir):
   '''Run rmdup on bam-files'''
   
   calls = []
   outfiles = []
   java_call = paths['java_home']+'java -jar '
   picard_cmd = paths['picard_home'] + 'MarkDuplicates.jar'
   
   for i,f in enumerate(infiles):
      out_bam = f + '.rmdup.bam'
      outfiles.append(out_bam)
      metrics_file = f + '.metrics.txt'
      picard_arg = ' INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=%s VALIDATION_STRINGENCY=LENIENT' % (f, out_bam, metrics_file, tmpdir)
      call = java_call + picard_cmd + picard_arg
      calls.append(call)
   return (calls, outfiles)
   

parser = argparse.ArgumentParser(description='''
   Process BAM alignment files:
   Reads lib file that should contain (one line pr. bamfile)
   bam_file\tqual_threshold\tlib_name\n
   
   filter for mapping quality
   sort
   merge to libraries
   remove pcr duplicates
   merge to final bam
   ''')

# add the arguments
parser.add_argument('--lib', help='input file with libraries')
parser.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/')
parser.add_argument('--q', help='queue to submit jobs to (cbs, urgent) [cbs]', default='cbs')
parser.add_argument('--o', help='output file [final.flt.sort.rmdup.bam]', default='final.flt.sort.rmdup.bam')
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--i Kleb-10-213361.bam --lib kleb.lib'.split())
#args = parser.parse_args('--lib NA12891.libs'.split())

# set logging
logger = logging.getLogger('genobox_bam_process.py')
hdlr = logging.FileHandler('genobox_bam_process.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set queueing
paths = pipelinemod.setSystem()
home = os.getcwd()
cpuA = 'nodes=1:ppn=1,mem=512mb'
cpuC = 'nodes=1:ppn=1,mem=2gb'
cpuE = 'nodes=1:ppn=1,mem=5gb'
cpuF = 'nodes=1:ppn=2,mem=2gb'
cpuB = 'nodes=1:ppn=16,mem=10gb'

# read library file
(bam2lib,lib2bam) = read_libs(args.lib)


## CREATE CALLS ##

# filter bam and sort
(filter_sort_calls, filter_sort_files) = bam_filter_sort(lib2bam, bam2lib)

# merge to libs
(merge_lib_calls, lib_files) = merge_bam(lib2bam.keys(), lib2bam.values(), add_suffix=True, final_suffix='.flt.sort.bam')

# rmdup on libs
(rmdup_calls, rmdup_files) = rmdup(lib_files, args.tmpdir)

# merge to final file
(merge_final_call, final_file) = merge_bam([args.o], [rmdup_files], add_suffix=False)



## SUBMIT JOBS ##

filter_sort_ids = pipelinemod.submitjob(filter_sort_calls, home, paths, logger, 'run_filtersort', args.q, cpuC, False)
merge_lib_ids = pipelinemod.submitjob(merge_lib_calls, home, paths, logger, 'run_lib_merge', args.q, cpuC, True, 'complex', map(len, lib2bam.values()), True, *filter_sort_ids)
rmdup_ids = pipelinemod.submitjob(rmdup_calls, home, paths, logger, 'run_rmdup', args.q, cpuC, True, 'one2one', 1, True, *merge_lib_ids)
merge_final_ids = pipelinemod.submitjob(merge_final_call, home, paths, logger, 'run_final_merge', args.q, cpuC, True, 'conc', len(rmdup_ids), True, *rmdup_ids)


# release jobs #
allids = []
allids.extend(filter_sort_ids) ; allids.extend(merge_lib_ids) ; allids.extend(rmdup_ids) ; allids.extend(merge_final_ids)
releasemsg = pipelinemod.releasejobs(allids)

