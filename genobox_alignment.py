#!/panvol1/simon/bin/python2.7

def check_fa(fa, bwa6):
   '''Checks for a fa of the input fasta file. If not present creates it'''
   
   import genobox_modules
   import subprocess
   import os
   import sys
   
   paths = genobox_modules.setSystem()
   
   # check if fa exists
   if bwa6:
      index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
   else:
      index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.rbwt', '.rpac', '.rsa', '.sa']
   
   for suf in index_suffixes:
      f = fa + suf
      if os.path.exists(f):
         pass
      else:
         sys.stderr.write('%s not found, creating bwa index\n' % fa)
         if bwa6:
            call = paths['bwa_6_2_home'] + 'bwa index -a is %s' % fa
         else:
            call = paths['bwa_home'] + 'bwa index -a is %s' % fa
         try: 
            subprocess.check_call(call, shell=True)
         except:
            sys.stderr.write('bwa index -a is failed, trying bwa index -a bwtsw\n')
            if bwa6:
               call = paths['bwa_6_2_home'] + 'bwa index -a bwtsw %s' % fa
            else:
               call = paths['bwa_home'] + 'bwa index -a bwtsw %s' % fa
            try:
               subprocess.check_call(call, shell=True)
            except: 
               raise TypeError('bwa index could not be created from %s' % fa)
         break
 

def check_formats_fq(i, gz, bwa6):
   '''Checks format of fastq file and returns it'''
   
   import genobox_modules
   
   # check if fastq and if so mode
   format = genobox_modules.set_filetype(i, gz)
   if format != 'fastq':
      raise ValueError('Input must be fastq')
   else:
      fqtype = genobox_modules.set_fqtype(i, gz)   
   return fqtype


def all_same(items):
   '''Check if all items in list/tuple are the same type'''
   return all(x == items[0] for x in items)


def check_trim(args):
   '''Check if args.no_trim attribute is present and return files to use'''
   
   # check if args has no_trim attribute (will not if the module is alignment)
   present = False
   try:
      args.no_trim
   except: return ([args.se, args.pe1, args.pe2])
   
   # return if it was set
   if args.no_trim == True:
      return ([args.se, args.pe1, args.pe2])
   else:
      return ([args.se, args.pe1, args.pe2])
      #return args.trimmed_files     
      

def bwa_se_align(fastqs, fa, fqtypes, qtrim, N, alignpath, bwa6, library, threads, queue, add_aln, partition, logger):
   '''Start alignment using bwa of fastq reads on index'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   import os
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=7gb,walltime=345600'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=345600'
   if threads != 1:
      if partition == 'uv' or partition == 'uv2':
         cpuB = 'procs=%s,mem=5gb,walltime=345600,flags=sharedmem' % threads
      else:
         if threads > 8:
            cpuB = 'nodes=1:ppn=%s,mem=7gb,walltime=345600' % threads
         else:
            cpuB = 'nodes=1:ppn=%s,mem=5gb,walltime=345600' % threads
   else:
      cpuB = cpuA
   
   # get readgroups
   RG=library.getRG('Data')
   #RG = genobox_modules.read_groups_from_libfile('Data', library)
   
   # align
   if bwa6:
      cmd = paths['bwa_6_2_home'] + 'bwa aln '
   else:
      cmd = paths['bwa_home'] + 'bwa aln '
   
   if add_aln: cmd = cmd + add_aln
   bwa_align = []
   saifiles = []
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      saifile = alignpath + f + '.sai'
      saifiles.append(saifile)
      if fqtypes[i] == 'Illumina':
         arg = ' -I -t %i -q %i %s %s > %s' % (threads, qtrim, fa, fq, saifile)
      elif fqtypes[i] == 'Sanger':
         arg = ' -t %i -q %i %s %s > %s' % (threads, qtrim, fa, fq, saifile)   
      elif fqtypes[i] == 'Solexa':
         raise ValueError('File %s is in Solexa format, convert to Sanger first\n' % fq)
      bwa_align.append(cmd+arg)
   
   # samse
   bwa_samse = []
   bamfiles = []
   bamfiles_dict = dict()
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      bamfile = alignpath + f + '.bam'
      bamfiles.append(bamfile)
      bamfiles_dict[fq] = bamfile
      if bwa6:
         p = paths['bwa_6_2_home']
      else:
         p = paths['bwa_home']
      
      call = '%sbwa samse -n %i -r \"%s\" %s %s %s | %ssamtools view -Sb - > %s' % (p, N, '\\t'.join(RG[fq]), fa, saifiles[i], fq, paths['samtools_home'], bamfile)
      bwa_samse.append(call)
      
   
   # submit jobs
   # create moab instance for the align_calls and dispatch to queue
   bwa_align_moab = Moab(bwa_align, logfile=logger, runname='run_genobox_bwaalign', queue=queue, cpu=cpuB, partition=partition)
   bwa_samse_moab = Moab(bwa_samse, logfile=logger, runname='run_genobox_bwasamse', queue=queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=bwa_align_moab.ids, partition=partition)
      
   # release jobs
   print "Releasing jobs"
   #bwa_align_moab.release()
   #bwa_samse_moab.release()
      
   return (bwa_samse_moab.ids, bamfiles_dict)


def bwa_pe_align(pe1, pe2, fa, fqtypes_pe1, fqtypes_pe2, qtrim, N, alignpath, bwa6, a, library, threads, queue, add_aln, partition, logger):
   '''Start alignment using bwa of paired end fastq reads on index'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   import os
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=7gb,walltime=345600'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=345600'
   if threads != 1:
      if partition == 'uv' or partition == 'uv2':
         cpuB = 'procs=%s,mem=5gb,walltime=172800,flags=sharedmem' % threads
      else:
         if threads > 8:
            cpuB = 'nodes=1:ppn=%s,mem=7gb,walltime=345600' % threads
         else:
            cpuB = 'nodes=1:ppn=%s,mem=5gb,walltime=345600' % threads
   else:
      cpuB = cpuA
   
   # get readgroups
   RG=library.getRG('Data')
   #RG = genobox_modules.read_groups_from_libfile('Data', library)
   
   # align and sampe
   if bwa6:
      cmd = paths['bwa_6_2_home'] + 'bwa '
   else:
      cmd = paths['bwa_home'] + 'bwa '
   
   bwa_align = []
   sam2bam_calls = []
   bwa_align1_calls = []
   bwa_align2_calls = []
   bwa_sampe_calls = []
   
   saifiles1 = []
   saifiles2 = []
   bamfiles = []
   bamfiles_dict = dict()
   
   for i,fq in enumerate(pe1):
      # set input fastq format
      if fqtypes_pe1[i] != fqtypes_pe2[i]:
         raise ValueError('Fastq formats are not the same for %s and %s' % (pe1[i], pe2[i]))
      elif fqtypes_pe1[i] == 'Sanger':
         bwa_cmd = '%s aln' % cmd
      elif fqtypes_pe1[i] == 'Illumina':
         bwa_cmd = '%s aln -I ' % cmd
      else:
         raise ValueError('fqtype must be Sanger or Illumina')
      
      if add_aln: bwa_cmd = bwa_cmd + add_aln
      
      # set filenames
      f1 = os.path.split(pe1[i])[1]
      f2 = os.path.split(pe2[i])[1]
      saifile1 = alignpath + f1 + '.sai'
      saifile2 = alignpath + f2 + '.sai'
      saifiles1.append(saifile1)
      saifiles2.append(saifile2)
      bamfiles.append(alignpath + f1 + '.bam')
      
      bamfiles_dict[pe1[i]] = alignpath + f1 + '.bam'
      bamfiles_dict[pe2[i]] = alignpath + f1 + '.bam'
      
      if bwa6:
         p = paths['bwa_6_2_home']
      else:
         p = paths['bwa_home']
      
      # generate calls
      bwa_align1 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, threads, qtrim, fa, saifiles1[i], pe1[i])
      bwa_align2 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, threads, qtrim, fa, saifiles2[i], pe2[i])
      sampecall = '%sbwa sampe -n %i -a %i -r \"%s\" %s %s %s %s %s | %ssamtools view -Sb - > %s' % (p, N, a, '\\t'.join(RG[fq]), fa, saifiles1[i], saifiles2[i], pe1[i], pe2[i], paths['samtools_home'], bamfiles[i])
      bwa_align1_calls.append(bwa_align1)
      bwa_align2_calls.append(bwa_align2)
      bwa_sampe_calls.append(sampecall)
   
   # submit jobs
   # create moab instance for the align_calls and dispatch to queue
   bwa_align1_moab = Moab(bwa_align1_calls, logfile=logger, runname='run_genobox_bwaalign1', queue=queue, cpu=cpuB, partition=partition)
   bwa_align2_moab = Moab(bwa_align2_calls, logfile=logger, runname='run_genobox_bwaalign2', queue=queue, cpu=cpuB, partition=partition)
   
   # set jobids in the correct way
   bwa_alignids = []
   for i in range(len(bwa_align1_moab.ids)):
      bwa_alignids.append(bwa_align1_moab.ids[i])
      bwa_alignids.append(bwa_align2_moab.ids[i])
   
   # submit sampe
   bwa_sampe_moab = Moab(bwa_sampe_calls, logfile=logger, runname='run_genobox_bwasampe', queue=queue, cpu=cpuA, depend=True, depend_type='conc', depend_val=[2], depend_ids=bwa_alignids, partition=partition)
   
   # release jobs
   print "Releasing jobs"
   #bwa_align1_moab.release()
   #bwa_align2_moab.release()
   #bwa_sampe_moab.release()
      
   return (bwa_sampe_moab.ids, bamfiles_dict)


def bwasw_pacbio(fastqs, fa, fqtypes, alignpath, bwa6, library, threads, queue, partition, logger):
   '''Start alignment of fastq files using BWA-SW Pacific Biosciences data'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   import os
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=7gb,walltime=345600'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=345600'
   if threads != 1:
      if partition == 'uv' or partition == 'uv2':
         cpuB = 'procs=%s,mem=5gb,walltime=345600,flags=sharedmem' % threads
      else:
         if threads > 8:
            cpuB = 'nodes=1:ppn=%s,mem=7gb,walltime=345600' % threads
         else:
            cpuB = 'nodes=1:ppn=%s,mem=5gb,walltime=345600' % threads
   else:
      cpuB = cpuA
   
   # align
   if bwa6:
      cmd = paths['bwa_6_2_home'] + 'bwa'
   else:
      cmd = paths['bwa_home'] + 'bwa'
   
   bwa_align = []
   bamfiles = []
   bamfiles_dict = dict()
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      bamfile = alignpath + f + '.bam'
      bamfiles.append(bamfile)
      bamfiles_dict[fq] = bamfile      
      if fqtypes[i] == 'Illumina':
         raise ValueError('BWA-SW should not align reads with Illumina Qualities')
      elif fqtypes[i] == 'Sanger':
         arg = ' bwasw -t %i -b5 -q2 -r1 -z10 %s %s |  %ssamtools view -Sb - > %s' % (threads, fa, fq, paths['samtools_home'], bamfile)
      bwa_align.append(cmd+arg)
   
   # submit jobs
   # create moab instance for the align_calls and dispatch to queue
   bwa_align_moab = Moab(bwa_align, logfile=logger, runname='run_genobox_bwaalign', queue=queue, cpu=cpuB, partition=partition)
   
   # release jobs
   print "Releasing jobs"
   #bwa_align_moab.release()
   
   return (bwa_align_moab.ids, bamfiles_dict)

def bwasw_iontorrent(fastqs, fa, fqtypes, alignpath, bwa6, library, threads, queue, partition, logger):
   '''Start alignment of fastq files using BWA-SW Iontorrent data'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   import os
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   
   # setting cpus
   cpuA = 'nodes=1:ppn=1,mem=7gb,walltime=345600'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=345600'
   if threads != 1:
      if partition == 'uv' or partition == 'uv2':
         cpuB = 'procs=%s,mem=5gb,walltime=345600,flags=sharedmem' % threads
      else:
         if threads > 8:
            cpuB = 'nodes=1:ppn=%s,mem=7gb,walltime=345600' % threads
         else:
            cpuB = 'nodes=1:ppn=%s,mem=5gb,walltime=345600' % threads
   else:
      cpuB = cpuA
   
   # align
   if bwa6:
      cmd = paths['bwa_6_2_home'] + 'bwa '
   else:
      cmd = paths['bwa_home'] + 'bwa '
   
   bwa_align = []
   bamfiles = []
   bamfiles_dict = dict()
   for i,fq in enumerate(fastqs):
      f = os.path.split(fq)[1]
      bamfile = alignpath + f + '.bam'
      bamfiles.append(bamfile)
      bamfiles_dict[fq] = bamfile      
      if fqtypes[i] == 'Illumina':
         raise ValueError('BWA-SW should not align reads with Illumina Qualities')
      elif fqtypes[i] == 'Sanger':
         arg = ' bwasw -t %i %s %s |  %ssamtools view -Sb - > %s' % (threads, fa, fq, paths['samtools_home'], bamfile)
      bwa_align.append(cmd+arg)
   
   # submit jobs
   # create moab instance for the align_calls and dispatch to queue
   bwa_align_moab = Moab(bwa_align, logfile=logger, runname='run_genobox_bwaalign', queue=queue, cpu=cpuB, partition=partition)
   
   # release jobs
   print "Releasing jobs"
   #bwa_align_moab.release()
   
   return (bwa_align_moab.ids, bamfiles_dict)


def start_alignment(args, logger):
   '''Start alignment of fastq files using BWA'''
   
   import genobox_modules
   from genobox_classes import Semaphore, Library
   import subprocess
   import os
   import random
   import string
   
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   semaphore_ids = []
   bamfiles = dict()
   
   if not os.path.exists('alignment'):
      os.makedirs('alignment')
   
   # initialize library file from given arguments (if args.mapq is defined then its called from abgv, else it is called from alignment)
   if hasattr(args, 'mapq'):
      library = genobox_modules.initialize_library(args.libfile, args.se, args.pe1, args.pe2, args.sample, args.mapq, args.libs, args.pl)
   else:
      library = genobox_modules.initialize_library(args.libfile, args.se, args.pe1, args.pe2, args.sample, [30], args.libs, args.pl)
   
   # check for fa
   check_fa(args.fa, args.bwa6)
   
   # check for if trimming was performed (abgv only) and set correct files
   #(se_files, pe1_files, pe2_files) = check_trim(args)
   
   # start single end alignments
   if args.se:
            
      # get platform info
      (PL, PL2data) = library.getPL('Data')      
      
      print "Submitting single end alignments"
      for key,value in PL2data.items():
         if key == 'ILLUMINA' or key == 'HELICOS':
            fqtypes_se = []
            # filter to only contain single end files
            toalign = []
            for v in value:
               if v in args.se: toalign.append(v)
            for fq in toalign: 
               if args.quals:
                  fqtypes_se.append(args.quals)
               else:
                  fqtypes_se.append(check_formats_fq(fq, args.gz, args.bwa6))
            
            # submit
            (se_align_ids, bamfiles_se) = bwa_se_align(toalign, args.fa, fqtypes_se, args.qtrim, args.N, 'alignment/', args.bwa6, library, args.n, args.queue, args.add_aln, args.partition, logger)
            semaphore_ids.extend(se_align_ids)
            bamfiles.update(bamfiles_se)
         elif key == 'PACBIO':
            toalign = []
            for v in value:
               if v in args.se: toalign.append(v)
            fqtypes_se = []
            for fq in toalign: 
               if args.quals:
                  fqtypes_se.append(args.quals)
               else:
                  fqtypes_se.append(check_formats_fq(fq, args.gz, args.bwa6))
            
            # submit
            (se_align_ids, bamfiles_se) = bwasw_pacbio(toalign, args.fa, fqtypes_se, 'alignment/', args.bwa6, library, args.n, args.queue, args.partition, logger)
            semaphore_ids.extend(se_align_ids)
            bamfiles.update(bamfiles_se)
         elif key == 'IONTORRENT' or key == '454':
            toalign = []
            for v in value:
               if v in args.se: toalign.append(v)
            fqtypes_se = []
            for fq in toalign: 
               if args.quals:
                  fqtypes_se.append(args.quals)
               else:
                  fqtypes_se.append(check_formats_fq(fq, args.gz, args.bwa6))
            
            # submit
            (se_align_ids, bamfiles_se) = bwasw_iontorrent(toalign, args.fa, fqtypes_se, 'alignment/', args.bwa6, library, args.n, args.queue, args.partition, logger)
            semaphore_ids.extend(se_align_ids)
            bamfiles.update(bamfiles_se)
   
   # start paired end alignments
   if args.pe1:
      if len(args.pe1) != len(args.pe2):
         raise ValueError('Same number of files must be given to --pe1 and --pe2')
      
      # set fqtypes
      fqtypes_pe1 = []
      fqtypes_pe2 = []
      for fq in args.pe1: 
         if args.quals:
            fqtypes_pe1.append(args.quals)
         else:
            fqtypes_pe1.append(check_formats_fq(fq, args.gz, args.bwa6))
      
      for fq in args.pe2:
         if args.quals:
            fqtypes_pe2.append(args.quals)
         else:
            fqtypes_pe2.append(check_formats_fq(fq, args.gz, args.bwa6))
      
      # submit
      print "Submitting paired end alignments"
      (pe_align_ids, bamfiles_pe) = bwa_pe_align(args.pe1, args.pe2, args.fa, fqtypes_pe1, fqtypes_pe2, args.qtrim, args.N, 'alignment/', args.bwa6, args.a, library, args.n, args.queue, args.add_aln, args.partition, logger)            
      semaphore_ids.extend(pe_align_ids)
      bamfiles.update(bamfiles_pe)
   
   # update library
   library.update_with_tag('Data', 'BAM', bamfiles, True)
   
   # wait for jobs to finish
   print "Waiting for jobs to finish ..." 
   
   s = Semaphore(semaphore_ids, home, 'bwa_alignment', args.queue, 60, 345600)
   s.wait()
   
   print "--------------------------------------"
   
   # return bamfiles   
   return (bamfiles, library)
