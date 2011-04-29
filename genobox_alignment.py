#!/usr/bin/python

def check_formats_fq(i):
   '''Checks format of fastq file and returns it'''
   
   import genobox_moab
   
   # check if fastq and if so mode
   format = genobox_moab.set_filetype(i)
   if format != 'fastq':
      raise ValueError('Input must be fastq')
   else:
      fqtype = genobox_moab.set_fqtype(i)   
   return fqtype

def all_same(items):
   '''Check if all items in list/tuple are the same type'''
   return all(x == items[0] for x in items)


def bwa_se_align(fastqs, bwaindex, fqtypes, qtrim, alignpath, r, threads, queue, logger):
   '''Start alignment using bwa of fastq reads on index'''
   
   import subprocess
   import genobox_moab
   import os
   paths = genobox_moab.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=5gb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   if threads != 1:
      cpuB = 'nodes=1:ppn=%s,mem=5gb' % threads
   else:
      cpuB = cpuA
   
   # align
   cmd = paths['bwa_home'] + 'bwa'
   bwa_align = []
   saifiles = []
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      saifile = alignpath + f + '.sai'
      saifiles.append(saifile)
      if fqtypes[i] == 'Illumina':
         arg = ' aln -I -t %i -q %i %s %s > %s' % (threads, qtrim, bwaindex, fq, saifile)
      elif fqtypes[i] == 'Sanger':
         arg = ' aln -t %i -q %i %s %s > %s' % (threads, qtrim, bwaindex, fq, saifile)   
      bwa_align.append(cmd+arg)
   
   # samse
   bwa_samse = []
   bamfiles = []
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      bamfile = alignpath + f + '.bam'
      bamfiles.append(bamfile)
      call = '%sbwa samse -r %s %s %s %s | %ssamtools view -Sb - > %s' % (paths['bwa_home'], r, bwaindex, saifiles[i], fq, paths['samtools_home'], bamfile)
      bwa_samse.append(call)
      
   # submit jobs
   allids = []
   
   bwa_alignids = genobox_moab.submitjob(bwa_align, home, paths, logger, 'run_genobox_bwaalign', queue, cpuB, False)
   bwa_samseids = genobox_moab.submitjob(bwa_samse, home, paths, logger, 'run_genobox_bwasamse', queue, cpuA, True, 'one2one', 1, True, *bwa_alignids)
   
   allids.extend(bwa_alignids) ; allids.extend(bwa_samseids)
   
   # release jobs
   releasemsg = genobox_moab.releasejobs(allids)
   
   return (bwa_samseids, bamfiles)

def bwa_pe_align(pe1, pe2, bwaindex, fqtypes_pe1, fqtypes_pe2, qtrim, alignpath, a, r, threads, queue, logger):
   '''Start alignment using bwa of paired end fastq reads on index'''
   
   import subprocess
   import genobox_moab
   import os
   paths = genobox_moab.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=5gb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   if threads != 1:
      cpuB = 'nodes=1:ppn=%s,mem=5gb' % threads
   else:
      cpuB = cpuA
   
   # align and sampe
   cmd = paths['bwa_home'] + 'bwa'
   bwa_align = []
   sam2bam_calls = []
   bwa_align1_calls = []
   bwa_align2_calls = []
   bwa_sampe_calls = []
   
   saifiles1 = []
   saifiles2 = []
   bamfiles = []
   
   for i,fq in enumerate(pe1):
      # set input fastq format
      if fqtypes_pe1[i] != fqtypes_pe2[i]:
         raise ValueError('Fastq formats are not the same for %s and %s' % (pe1[i], pe2[i]))
      elif fqtypes_pe1[i] == 'Sanger':
         bwa_cmd = '%sbwa aln ' % paths['bwa_home']
      elif fqtypes_pe1[i] == 'Illumina':
         bwa_cmd = '%sbwa aln -I ' % paths['bwa_home']
      else:
         raise ValueError('fqtype must be Sanger or Illumina')
      
      # set filenames
      f1 = os.path.split(pe1[i])[1]
      f2 = os.path.split(pe2[i])[1]
      saifile1 = alignpath + f1 + '.sai'
      saifile2 = alignpath + f2 + '.sai'
      saifiles1.append(saifile1)
      saifiles2.append(saifile2)
      bamfiles.append(alignpath + f1 + '.bam')
      
      # generate calls
      bwa_align1 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, threads, qtrim, bwaindex, saifiles1[i], pe1[i])
      bwa_align2 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, threads, qtrim, bwaindex, saifiles2[i], pe2[i])
      sampecall = '%sbwa sampe -a %i -r %s %s %s %s %s %s | %ssamtools view -Sb - > %s' % (paths['bwa_home'], a, r, bwaindex, saifiles1[i], saifiles2[i], pe1[i], pe2[i], paths['samtools_home'], bamfiles[i])
      bwa_align1_calls.append(bwa_align1)
      bwa_align2_calls.append(bwa_align2)
      bwa_sampe_calls.append(sampecall)
      
   
   # submit jobs
   allids = []
   bwa_alignids1 = genobox_moab.submitjob(bwa_align1_calls, home, paths, logger, 'run_genobox_bwaalign1', queue, cpuB, False)
   bwa_alignids2 = genobox_moab.submitjob(bwa_align2_calls, home, paths, logger, 'run_genobox_bwaalign2', queue, cpuB, False)
   
   # set jobids in the correct way
   bwa_alignids = []
   for i in range(len(bwa_alignids1)):
      bwa_alignids.append(bwa_alignids1[i])
      bwa_alignids.append(bwa_alignids2[i])
   
   # submit sampe
   bwa_sampeids = genobox_moab.submitjob(bwa_sampe_calls, home, paths, logger, 'run_genobox_bwasampe', queue, cpuA, True, 'conc', len(bwa_alignids), True, *bwa_alignids)
   
   # release jobs
   allids.extend(bwa_alignids) ; allids.extend(bwa_sampeids)
   releasemsg = genobox_moab.releasejobs(allids)
   
   return (bwa_sampeids, bamfiles)
   

def start_alignment(args, logger):
   '''Start alignment of fastq files using BWA'''
   
   import genobox_moab
   import subprocess
   import os
   paths = genobox_moab.setSystem()
   home = os.getcwd()
   semaphore_ids = []
   bamfiles = []
   
   if not os.path.exists('alignment'):
      os.makedirs('alignment')
   
   # start single end alignments
   if args.se:
      
      # set fqtypes
      fqtypes_se = map(check_formats_fq, args.se)
      print "Submitting single end alignments"
      (se_align_ids, bamfiles_se) = bwa_se_align(args.se, args.bwaindex, fqtypes_se, args.qtrim, 'alignment/', args.r, args.n, args.queue, logger)
      semaphore_ids.extend(se_align_ids)
      bamfiles.extend(bamfiles_se)
      
   # start paired end alignments
   if args.pe1:
      if len(args.pe1) != len(args.pe2):
         raise ValueError('Same number of files must be given to --pe1 and --pe2')
            
      # set fqtypes
      fqtypes_pe = []
      fqtypes_pe1 = map(check_formats_fq, args.pe1)
      fqtypes_pe2 = map(check_formats_fq, args.pe2)
      fqtypes_pe.extend(fqtypes_pe1)
      fqtypes_pe.extend(fqtypes_pe2)
      
      print "Submitting paired end alignments"
      (pe_align_ids, bamfiles_pe) = bwa_pe_align(args.pe1, args.pe2, args.bwaindex, fqtypes_pe1, fqtypes_pe2, args.qtrim, 'alignment/', args.a, args.r, args.n, args.queue, logger)            
      semaphore_ids.extend(pe_align_ids)
      bamfiles.extend(bamfiles_pe)
   
   # wait for jobs to finish
   print "Waiting for jobs to finish ..." 
   genobox_moab.wait_semaphore(semaphore_ids, home, 'bwa_alignment', args.queue, 60, 86400)
   print "--------------------------------------"
      
   # return bamfiles   
   return bamfiles