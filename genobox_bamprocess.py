#!/usr/bin/python

# filter bam and sort
def bam_filter_sort(lib2bam, bam2lib, m=500000000):
   '''Filter bam on quality and sort on stream'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_samFilterSort.py'
   calls = []
   
   # set infiles and outfiles
   infiles = lib2bam.values()
   
   calls = []
   outfiles = []
   for libfiles in infiles:
      for f in libfiles:
         sort_prefix = f + '.flt.sort'
         out_bam = f + '.flt.sort.bam'
         outfiles.append(out_bam)
         arg = ' --i %s --q %s --m %i --o %s' % (f, bam2lib[f][0], m, sort_prefix)
         calls.append(cmd+arg)
   return (calls, outfiles)

# merge to libraries
def merge_bam(libs, lib_infiles, add_suffix=False, final_suffix=''):
   '''Merge bam files to libraries'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
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
      
      # add suffix to outfile if set and add alignment to path if it is not already there
      if lib.startswith('alignment/'):
         out_bam = lib + final_suffix
      else:
         out_bam = 'alignment/' + lib + final_suffix
      outfiles.append(out_bam)
      
      if len(list_bams) == 1:
         call = 'cp %s %s' % (' '.join(list_bams), out_bam)
      else:
         sam_cmd = paths['samtools_home'] + 'samtools merge'
         sam_arg = ' %s %s' % (out_bam, ' '.join(list_bams))
         call = sam_cmd+sam_arg
      calls.append(call)
   return (calls, outfiles)

# rmdup
def rmdup(infiles, tmpdir):
   '''Run rmdup on bam-files'''
   
   import genobox_modules   
   paths = genobox_modules.setSystem()
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


def start_bamprocess(libfile, bams, mapq, libs, tmpdir, queue, final_bam, sample, logger):
   '''Starts bam processing of input files'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore
   import os
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   # create libfile if not given
   if not libfile:
      libfile = genobox_modules.library_from_input(bams, sample, mapq, libs)
   (bam2lib,lib2bam) = genobox_modules.read_bam_libs(libfile)
   
   ## CREATE CALLS ##
   
   # filter bam and sort
   (filter_sort_calls, filter_sort_files) = bam_filter_sort(lib2bam, bam2lib, 1500000000)
   
   # merge to libs
   (merge_lib_calls, libfiles) = merge_bam(lib2bam.keys(), lib2bam.values(), add_suffix=True, final_suffix='.flt.sort.bam')
   
   # rmdup on libs
   (rmdup_calls, rmdup_files) = rmdup(libfiles, tmpdir)
   
   # merge to final file
   (merge_final_call, final_file) = merge_bam([final_bam], [rmdup_files], add_suffix=False)
   
   
   ## SUBMIT JOBS ##
   
   print "Submitting jobs"
   filtersort_moab = Moab(filter_sort_calls, logfile=logger, runname='run_genobox_filtersort', queue=queue, cpu=cpuC)
   mergelib_moab = Moab(merge_lib_calls, logfile=logger, runname='run_genobox_lib_merge', queue=queue, cpu=cpuC, depend=True, depend_type='complex', depend_val=map(len, lib2bam.values()), depend_ids=filtersort_moab.ids)
   rmdup_moab = Moab(rmdup_calls, logfile=logger, runname='run_genobox_rmdup', cpu=cpuC, depend=True, depend_type='one2one', depend_val=[1], depend_ids=mergelib_moab.ids)
   mergefinal_moab = Moab(merge_final_call, logfile=logger, runname='run_genobox_final_merge', queue=queue, cpu=cpuC, depend=True, depend_type='conc', depend_val=[len(rmdup_moab.ids)], depend_ids=rmdup_moab.ids)
   
   # release jobs #
   print "Releasing jobs"
   filtersort_moab.release()
   mergelib_moab.release()
   rmdup_moab.release()
   mergefinal_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..." 
   s = Semaphore(mergefinal_moab.ids, home, 'bam_processing', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   # return final bamfile
   return final_bam
