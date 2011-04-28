#!/usr/bin/python

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

#def make_lib_file(se, pe1, pe2, bams, q, lib):
#   '''Create lib_file from input paramters'''
#   
#   import os
#   
#   fh = open('libs', 'w')
#   if se:
#      se = map(os.path.abspath, se)
#      for i,fq in enumerate(se):
#         fh.write('%s\t%i\t%s\n' % ('alignment/' + os.path.split(se[i])[1] +'.bam', q[i], 'alignment/' + lib[i]))
#   if pe1:
#      pe1 = map(os.path.abspath, pe1)
#      for i,fq in enumerate(pe1):
#         fh.write('%s\t%i\t%s\n' % ('alignment/' + os.path.split(pe1[i])[1] +'.bam', q[i], 'alignment/' + lib[i]))
#   if bams:
#      bams = map(os.path.abspath, bams)
#      for i, bam in enumerate(bams):
#         fh.write('%s\t%i\t%s\n' % ('alignment/' + os.path.split(bam)[1], q[i], 'alignment/' + lib[i]))
#   fh.close()
#   return 'libs'

def make_lib_file(bams, q, lib):
   '''Create lib_file from input paramters'''
   
   import os
   import sys
   
   fh = open('libs', 'w')
   bams = map(os.path.abspath, bams)
   
   # check if shorter arguments is passed to q or lib
   # if so the first argument is reused
   if len(q) < len(bams):
      sys.stderr.write("warning: length of mapq is shorter than length of bamfiles, reusing %i\n" % q[0])   
      while len(q) < len(bams):
         q.append(q[0])
   if len(lib) < len(bams):
      sys.stderr.write("warning: length of lib names is shorter than length of bamfiles, reusing %s\n" % lib[0])   
      while len(lib) < len(bams):
         lib.append(lib[0])
   # create library file
   for i, bam in enumerate(bams):
      fh.write('%s\t%i\t%s\n' % ('alignment/' + os.path.split(bam)[1], q[i], 'alignment/' + lib[i]))
   fh.close()
   return 'libs'


# filter bam and sort
def bam_filter_sort(lib2bam, bam2lib, m=500000000):
   '''Filter bam on quality and sort on stream'''
   
   import pipelinemod
   paths = pipelinemod.setSystem()
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_samFilterSort.py'
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
         arg = ' --i %s --q %s --m %i --o %s' % (f, bam2lib[f][0], m, sort_prefix)
         calls.append(cmd+arg)
   return (calls, outfiles)

# merge to libraries
def merge_bam(libs, lib_infiles, add_suffix=False, final_suffix=''):
   '''Merge bam files to libraries'''
   
   import pipelinemod
   paths = pipelinemod.setSystem()
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
   
   import pipelinemod   
   paths = pipelinemod.setSystem()
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
   

def start_bamprocess(lib_file, bams, mapq, libs, tmpdir, queue, final_bam, logger):
   '''Starts bam processing of input files'''
   
   # Filter for mapping quality
   # Sort
   # Merge to libraries
   # Remove pcr duplicates
   # Merge to final bam
   #
   # Input either:
   #    input as bam, mapq and libs
   # OR
   #    input as args.se, args.pe1, args.pe2 and mapq + libs
   # OR
   #    a lib file: bam_file\tqual_threshold\tlib_name\n      
   
   
   import subprocess
   import pipelinemod
   import os
   
   # set queueing
   paths = pipelinemod.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
   
   # create lib_file if not given
   if lib_file:
      (bam2lib,lib2bam) = read_libs(lib_file)
   else:
      lib_file = make_lib_file(bams, mapq, libs)
      (bam2lib,lib2bam) = read_libs(lib_file)
   
   ## CREATE CALLS ##
   
   # filter bam and sort
   (filter_sort_calls, filter_sort_files) = bam_filter_sort(lib2bam, bam2lib, 1500000000)
   
   # merge to libs
   (merge_lib_calls, lib_files) = merge_bam(lib2bam.keys(), lib2bam.values(), add_suffix=True, final_suffix='.flt.sort.bam')
   
   # rmdup on libs
   (rmdup_calls, rmdup_files) = rmdup(lib_files, tmpdir)
   
   # merge to final file
   (merge_final_call, final_file) = merge_bam([final_bam], [rmdup_files], add_suffix=False)
   
   
   ## SUBMIT JOBS ##
   
   print "Submitting jobs"
   filter_sort_ids = pipelinemod.submitjob(filter_sort_calls, home, paths, logger, 'run_genobox_filtersort', queue, cpuC, False)
   merge_lib_ids = pipelinemod.submitjob(merge_lib_calls, home, paths, logger, 'run_genobox_lib_merge', queue, cpuC, True, 'complex', map(len, lib2bam.values()), True, *filter_sort_ids)
   rmdup_ids = pipelinemod.submitjob(rmdup_calls, home, paths, logger, 'run_genobox_rmdup', queue, cpuC, True, 'one2one', 1, True, *merge_lib_ids)
   merge_final_ids = pipelinemod.submitjob(merge_final_call, home, paths, logger, 'run_genobox_final_merge', queue, cpuC, True, 'conc', len(rmdup_ids), True, *rmdup_ids)
   
   
   # release jobs #
   allids = []
   allids.extend(filter_sort_ids) ; allids.extend(merge_lib_ids) ; allids.extend(rmdup_ids) ; allids.extend(merge_final_ids)
   releasemsg = pipelinemod.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..." 
   pipelinemod.wait_semaphore(merge_final_ids, home, 'bam_processing', queue, 20, 2*86400)
   print "--------------------------------------"
   
   # return final bamfile
   return final_bam