#!/panvol1/simon/bin/python2.7

# filter bam and sort
def bam_filter_sort(lib2bam, bam2lib, m=500000000):
   '''Filter bam on quality and sort on stream'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   cmd = paths['genobox_home'] + 'genobox_samFilterSort.py'
   calls = []
   
   # set infiles and outfiles
   infiles = lib2bam.values()
   
   calls = []
   outfiles = []
   for library in infiles:
      for f in library:
         sort_prefix = f + '.flt.sort.bam'
         out_bam = f + '.flt.sort.bam'
         outfiles.append(out_bam)
         arg = ' --i %s --q %s --m %i --o %s' % (f, bam2lib[f][0], m, sort_prefix)
         calls.append(cmd+arg)
   return (calls, outfiles)

# merge to libraries
def merge_bam(libs, lib_infiles, add_suffix=False, final_suffix='', tmpdir='/panvol1/simon/tmp/'):
   '''Merge bam files to libraries'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   calls = []
   outfiles = []
   java_call = paths['java_home']+'java -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms4500m -Xmx4500m -jar '
   picard_cmd = paths['picard_home'] + 'MergeSamFiles.jar'
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
         #sam_cmd = paths['samtools_home'] + 'samtools merge'
         #sam_arg = ' %s %s' % (out_bam, ' '.join(list_bams))
         #call = sam_cmd+sam_arg
         arg = ' INPUT=%s OUTPUT=%s TMP_DIR=%s ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT' % (' INPUT='.join(list_bams), out_bam, tmpdir)
         call = java_call + picard_cmd + arg
      calls.append(call)
   return (calls, outfiles)

# rmdup
def rmdup(infiles, tmpdir):
   '''Run rmdup on bam-files'''
   
   import genobox_modules   
   paths = genobox_modules.setSystem()
   calls = []
   outfiles = []
   java_call = paths['java_home']+'java -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms4500m -Xmx4500m -jar '
   picard_cmd = paths['picard_home'] + 'MarkDuplicates.jar'
   
   for i,f in enumerate(infiles):
      out_bam = f + '.rmdup.bam'
      outfiles.append(out_bam)
      metrics_file = f + '.metrics.txt'
      picard_arg = ' INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=%s VALIDATION_STRINGENCY=LENIENT' % (f, out_bam, metrics_file, tmpdir)
      call = java_call + picard_cmd + picard_arg
      calls.append(call)
   return (calls, outfiles)

def realign_bam(in_bam, out_bam, fa, known=None):
   '''Runs realignment of bam'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   calls = []
   
   gatk_cmd = paths['GATK_home'] + 'GenomeAnalysisTK.jar'
   java_call = paths['java_home']+'java -Djava.io.tmpdir=/panvol1/simon/tmp/ -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms4500m -Xmx4500m -jar %s ' % gatk_cmd
   
   realign_bam = out_bam.replace('.bam', '.realign.bam')
   
   # index bam
   cmd = paths['samtools_home'] + 'samtools '
   # adding pipe to make it being written as a shell-file so all commands are submitted at the same time (fix dependencies)
   arg = 'index %s | cat - ' % in_bam
   c = cmd+arg
   #calls.append(cmd+arg)
   
   # realigner target creator
   if known:
      arg = '-I %s -R %s -T RealignerTargetCreator -known %s -o %s' % (in_bam, fa, known, in_bam+'.intervals')
   else:
      arg = '-I %s -R %s -T RealignerTargetCreator -o %s' % (in_bam, fa, in_bam+'.intervals')
   #calls.append(java_call+arg)
   c = '%s\n\n%s%s' % (c, java_call, arg)
   
   # realignment step
   arg = '-I %s -T IndelRealigner -R %s -targetIntervals %s -o %s' % (in_bam, fa, in_bam+'.intervals', realign_bam)
   #calls.append(java_call+arg)
   c = '%s\n\n%s%s' % (c, java_call, arg)
   
   calls = [c]
   
   return (calls, realign_bam)


def start_bamprocess(library_file, bams, mapq, libs, tmpdir, queue, final_bam, realignment, known, fa, sample, partition, logger):
   '''Starts bam processing of input files'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab, Semaphore, Library
   import os
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=345600'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=345600'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=345600'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=345600'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=345600'
   cpuG = 'nodes=1:ppn=1,mem=6gb,walltime=345600'
   cpuH = 'nodes=1:ppn=2,mem=7gb,walltime=345600'
   
   # create library instance
   if library_file and library_file != 'None':
      if isinstance(library_file, Library):
         library = library_file
      else:
         library = Library(library_file)
         library.read()
   else:
      library = genobox_modules.initialize_library(libfile=library_file, sample=sample, mapq=mapq, libs=libs, bams=bams)
   
   (bam2lib, lib2bam) = library.getBamLibs()
      
   ## CREATE CALLS ##
   
   # filter bam and sort
   (filter_sort_calls, filter_sort_files) = bam_filter_sort(lib2bam, bam2lib, 1500000000)
   
   # merge to libs
   (merge_lib_calls, librarys) = merge_bam(lib2bam.keys(), lib2bam.values(), add_suffix=True, final_suffix='.flt.sort.bam', tmpdir=tmpdir)
   
   # rmdup on libs
   (rmdup_calls, rmdup_files) = rmdup(librarys, tmpdir)
   
   # optional: realignment
   if realignment:
      (merge_final_call, sample_file) = merge_bam([final_bam], [rmdup_files], add_suffix=False)
      (realign_calls, final_file) = realign_bam(final_bam, final_bam, fa, known)
   else:
      # merge to final file
      (merge_final_call, final_file) = merge_bam([final_bam], [rmdup_files], add_suffix=False)
   
   
   ## SUBMIT JOBS ##
   
   print "Submitting jobs"
   filtersort_moab = Moab(filter_sort_calls, logfile=logger, runname='run_genobox_filtersort', queue=queue, cpu=cpuH, partition=partition)
   mergelib_moab = Moab(merge_lib_calls, logfile=logger, runname='run_genobox_lib_merge', queue=queue, cpu=cpuE, depend=True, depend_type='complex', depend_val=map(len, lib2bam.values()), depend_ids=filtersort_moab.ids, partition=partition)
   rmdup_moab = Moab(rmdup_calls, logfile=logger, runname='run_genobox_rmdup', queue=queue, cpu=cpuG, depend=True, depend_type='one2one', depend_val=[1], depend_ids=mergelib_moab.ids, partition=partition)          # NB: If memory should be changed, also change java memory spec in rmdup function
   mergefinal_moab = Moab(merge_final_call, logfile=logger, runname='run_genobox_final_merge', queue=queue, cpu=cpuC, depend=True, depend_type='conc', depend_val=[len(rmdup_moab.ids)], depend_ids=rmdup_moab.ids, partition=partition)
   if realignment:
      realign_moab = Moab(realign_calls, logfile=logger, runname='run_genobox_realignment', queue=queue, cpu=cpuE, depend=True, depend_type='one2one', depend_val=[1], depend_ids=mergefinal_moab.ids, partition=partition)
   # realignment calls needs to be written together in a shell-file or dependent on each other #
   
   # release jobs #
   print "Releasing jobs"
   #filtersort_moab.release()
   #mergelib_moab.release()
   #rmdup_moab.release()
   #mergefinal_moab.release()
   #if realignment: realign_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..." 
   if realignment:
      s = Semaphore(realign_moab.ids, home, 'bam_processing', queue, 20, 345600)
   else:
      s = Semaphore(mergefinal_moab.ids, home, 'bam_processing', queue, 20, 345600)
   s.wait()
   print "--------------------------------------"
   
   # return final bamfile
   return final_bam
