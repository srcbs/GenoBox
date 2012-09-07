#!/panvol1/simon/bin/python2.7

# functions for use in the genobox pipeline #


import platform
import os
import sys
import re
import subprocess
import time

###   SET SYSTEM PATH   ###

def setSystem():
   '''Identifies architecture and sets system specific
      variables such as executable paths'''
      
   plat = platform.uname()
   
   sys = {}
   
   if plat[1] == 'amyloid.local':
      sys['host'] = 'amyloid.local'
      sys['bowtie_home'] = '/Users/simon/Work/bin/bowtie-0.12.5/'
      sys['bwa_home'] = '/Users/simon/Work/bin/bwa-0.5.8a/'
      sys['euler_home'] = '/Users/simon/Work/bin/euler-sr.1.1.2/'
      sys['hmm-3.0_home'] = '/Users/simon/Work/bin/hmmer-3.0/'
      sys['mothur_home'] = '/Users/simon/Work/bin/mothur/'
      sys['pyroNoise2_home'] = '/Users/simon/Work/bin/pyroNoise2/bin/'
      sys['rna_hmm3_home'] = '/Users/simon/Work/bin/rna_hmm3/'
      sys['samtools_home'] = '/Users/simon/Work/bin/samtools/'
      sys['velvet_home'] = '/Users/simon/Work/bin/velvet_0.7.63/'
      sys['pyscripts_home'] = '/Users/simon/Work/bin/python/pipeline/'
      sys['greengenes_home'] = '/Users/simon/Work/projects/metagenomics/alignment_databases/greengenes/'
      
   elif plat[1] == 'sbiology' or plat[1] == 'ibiology' or plat[1] == 'life' or plat[1] == 'interaction':
      sys['host'] = plat[1]
      sys['bowtie_home'] = '/home/people/simon/bin/bowtie-0.12.7/'
      sys['bwa_home'] = '/home/people/simon/bin/bwa-0.5.8a/'
      sys['euler_home'] = '/home/people/simon/bin/euler-sr.1.1.2/'
      sys['hmm-3.0_home'] = '/home/people/simon/bin/hmmer-3.0/'
      sys['mothur_home'] = '/home/people/simon/bin/mothur/'
      sys['pyroNoise2_home'] = '/home/people/simon/bin/pyroNoise2/bin/'
      #sys['rna_hmm3_home'] = '/home/people/simon/bin/rna_hmm3/'
      sys['samtools_home'] = '/home/people/simon/bin/samtools-0.1.12a/'
      sys['velvet_home'] = '/home/people/simon/bin/velvet_0.7.63/'
      sys['pyscripts_home'] = '/home/people/simon/bin/python/pipeline/'
      sys['maq_home'] = '/home/people/simon/bin/maq-0.7.1/'
      sys['greengenes_home'] = '/home/people/simon/projects/metagenomics/alignment_databases/greengenes/'
      sys['prodigal_home'] = ''
      sys['chain_index'] = '/panvol1/simon/databases/chainmap/'
      sys['hs_ref_GRCh37'] = '/panvol1/simon/databases/hs_ref37/'
      sys['hs_ref_ncbi36'] = '/panvol1/simon/databases/hs_ref36/'
      sys['java_home'] = '/usr/java/jre1.6.0_22/bin/'
      sys['picard_home'] = '/panvol1/simon/bin/picard-tools-1.26/'
      sys['blastdb'] = '/home/people/simon/projects/databases/blastdb/'
      sys['blastall_home'] = '/usr/cbs/bio/bin/linux64/'
      sys['taxonomy_ncbi'] = '/panvol1/simon/databases/taxonomy/'
      
   elif plat[4] == 'x86_64':
      sys['host'] = 'protein-s0'
      sys['bin_home'] = '/panvol1/simon/bin/'
      sys['python2.7_home'] = '/panvol1/simon/bin/'
      sys['java_home'] = '/panvol1/simon/bin/jre1.6.0_21/bin/'
      sys['phymm_home'] = '/panvol1/simon/bin/phymm/'
      sys['pyscripts_home'] = '/panvol1/simon/bin/pipeline/'
      sys['genobox_home'] = '/panvol1/simon/bin/genobox/'
      sys['blastall_home'] = '/panfs/saqqaq/bin/'
      sys['blastdb'] = '/panvol1/simon/databases/blastdb/'
      sys['blat_home'] = '/panfs/saqqaq/bin/'
      sys['fastxtoolkit_home'] = '/panvol1/simon/bin/'
      sys['bowtie_home'] = '/panvol1/simon/bin/bowtie-0.12.5/'
      sys['bwtindexes_home'] = '/panvol1/simon/bin/bowtie-0.12.5/indexes/'
      sys['hs_ref_GRCh37'] = '/panvol1/simon/databases/hs_ref37/'
      sys['hs_ref_ncbi36'] = '/panvol1/simon/databases/hs_ref36/'
      sys['bwa_home'] = '/panvol1/simon/bin/bwa-0.5.9/'
      sys['bwaindexes_home'] = '/panvol1/simon/bin/bwa-0.5.8a/indexes/'
      #sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.12.a/'
      sys['samtools_svn_home'] = '/panvol1/simon/bin/samtools_svn/'
      sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.18/'
      sys['picard_home'] = '/panvol1/simon/bin/picard-tools-1.26/'
      sys['bedtools_home'] = '/panvol1/simon/bin/BEDTools-Version-2.12.0/'
      sys['bedtools_bin'] = '/panvol1/simon/bin/'                          # 2-16.2
      sys['greengenes_home'] = '/panvol1/simon/databases/greengenes/'
      sys['mothur_home'] = '/panvol1/simon/bin/Mothur-1.11.0/'
      sys['rna_hmm3_home'] = '/panvol1/simon/bin/rna_hmm3/'
      sys['hmm-3.0_home'] = '/panvol1/simon/bin/hmmer-3.0/'
      sys['phymm_home'] = '/panvol1/simon/bin/phymm/'
      sys['prodigal_home'] = '/panvol1/simon/bin/prodigal/'
      sys['rnammer_home'] = '/panvol1/simon/bin/rnammer-1.2/'
      sys['cge_genomes'] = '/panvol1/simon/databases/cge_genomes/'
      sys['velvet_home'] = '/panvol1/simon/bin/velvet_1.1.02/'
      #sys['newbler'] = '/panfs/saqqaq/bin/newbler/bin/'
      sys['newbler'] = '/panvol1/simon/bin/454/bin/'
      sys['stampy'] = '/panvol1/simon/bin/stampy-1.0.6/'
      sys['chain_index'] = '/panvol1/simon/databases/chainmap/'
      sys['taxonomy_ncbi'] = '/panvol1/simon/databases/taxonomy/'
      sys['R_home'] = '/tools/bin/'
      sys['bwa_6_2_home'] = '/panvol1/simon/bin/bwa-0.6.2/'
      sys['GATK_home'] = '/panvol1/simon/bin/GenomeAnalysisTK-2.0-34/'
   else:
      raise ValueError('Platform not identified')
   
   return(sys)




#################################

###      JOB SUBMISSION       ###

#################################

def submit_create_dependencies(calls, depend, dependtype, n, a):
   '''Create job dependencies
      
      'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
      'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
      'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
      'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
   '''
   
   if not depend:
      depends = None
   else:
         if dependtype == 'one2one':
            depends = list(a)
         
         elif dependtype == 'expand':
            a = list(a)
            depends = []
            c = 0
            for j in range(len(calls)):
               c = c + 1
               if c < int(n):
                  depends.append(a[0])
               if c == int(n):
                  depends.append(a.pop(0))
                  c = 0
         
         elif dependtype == 'conc':
            a = list(a)
            depends = []
            for j in range(len(calls)):
               s = ':'.join(a[:int(n)])
               depends.append(s)
               a = a[int(n):]
         
         elif dependtype == 'complex':
            old_index = 0
            depends = []
            a = list(a)
            for index in n:
               s = ':'.join(a[old_index:(index+old_index)])
               depends.append(s)
               old_index=index
         else:
            raise AttributeError('dependtype not recognized: %s' %dependtype)
   return depends

def submit_xmsub(calls, home, paths, logger, runname, queue, cpu, depend, hold, depends):
   '''Submits jobs using xmsub'''
   
   import re
   import subprocess
   import time
   
   ids = []
   for i in range(len(calls)):
      call = calls[i]
      stdout = '%s/log/%s%i.o' % (home, runname, i)
      stderr = '%s/log/%s%i.e' % (home, runname, i)
      
      # catch stdouts if call includes 'program infile > outfile', needs to be directed as -O instead of >
      pattern = re.compile(r'(^.+)>\s(.+)$')
      match = pattern.search(call)
      if match:
         call = match.group(1)
         stdout = '%s/%s' % (home, match.group(2))
      
      # create xmsub commands
      cmd = '%sxmsub' % (paths['pyscripts_home'])
      
      # toggle if job should be on hold
      if hold:
         cmd = '%s -h ' % cmd
      if not depend:
         xmsub = cmd+' -d %s -l %s,walltime=172800 -O %s -E %s -r y -q %s -N %s -t %s' % (home, cpu, stdout, stderr, queue, runname, call)
      else:
         xmsub = cmd+' -d %s -l %s,walltime=172800,depend=%s -O %s -E %s -r y -q %s -N %s -t %s' % (home, cpu, depends[i], stdout, stderr, queue, runname, call)
      
      time.sleep(0.1)
      logger.info(xmsub)
                
      # submit
      try:
         id = subprocess.check_output(xmsub, shell=True)
      except:
         print 'Job error, waiting 1m'
         time.sleep(60)
         id = subprocess.check_output(xmsub, shell=True)
      ids.append(id.split('\n')[1])
   return(ids)


def submit_wrapcmd(calls, home, paths, logger, runname, queue, cpu, depend, hold, depends):
   '''Take input as command and submit to msub (this way pipes and redirects can be done)'''
   
   import subprocess
   import random
   import string
   import time
   
   ids = []
   for i in range(len(calls)):
      call = calls[i]
      N = 10
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))
      filename = 'pbsjob.tmp%s' % rand
      
      stdout = '%s/log/%s%i.o' % (home, runname, i)
      stderr = '%s/log/%s%i.e' % (home, runname, i)
      
      # write pbsjob file
      fh = open(filename, 'w')
      fh.write('#!/bin/sh\n\n')
      fh.write('%s\n' % call)
      fh.close()
      
      # create msub command
      cmd = 'msub'
      if hold:
         cmd = '%s -h' % cmd
      if not depend:
         msub = '%s -d %s -l %s,walltime=172800 -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, cpu, stdout, stderr, queue, runname, filename)
      else:
         msub = '%s -d %s -l %s,walltime=172800,depend=%s -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, cpu, depends[i], stdout, stderr, queue, runname, filename)
      
      time.sleep(0.1)
      logger.info(msub)
      
      # submit
      try:
         id = subprocess.check_output(msub, shell=True)
      except:
         print 'Job error, waiting 1m'
         time.sleep(60)
         id = subprocess.check_output(xmsub, shell=True)
      ids.append(id.split('\n')[1])
      
      # remove pbsjob file
      #rm_files([filename])
   return ids


def submitjob(calls, home, paths, logger, runname, queue='cbs', cpu='nodes=1:ppn=16,mem=8gb', depend=False, dependtype='one2one', n=[1], hold=True, *a):
   '''Submit job to cluster, if depend=True dependent jobids should be given as *a. 
      Dependtype can be 'one2one', 'expand' or 'conc', this allows flexibility for 
      when running jobs where files have been split (eg. blast runs)
      
      'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
      'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
      'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
      'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...'''
   
   # create job dependencies
   depends = submit_create_dependencies(calls, depend, dependtype, n, a)
   
   # temporary fix to group problems on cge-s2 cluster
   if calls[0].find('|') > -1 or calls[0].find('<') > -1:
      # perform wrapcmd if calls includes pipes / left-redirects
      ids = submit_wrapcmd(calls, home, paths, logger, runname, queue, cpu, depend, hold, depends)
   else:
      # perform xmsub if calls does not include pipes (can have right-redirects)
      ids = submit_xmsub(calls, home, paths, logger, runname, queue, cpu, depend, hold, depends)
   return(ids)


def releasejobs(jobids):
   '''Release jobids from hold and start running'''
   
   print  "Releasing jobs"
   while len(jobids) > 0:
      if len(jobids) > 199:
         cmd = 'mjobctl -u user \"%s\"' % (' ').join(jobids[:199])
         del jobids[:199]
         out = subprocess.check_output(cmd, shell=True)
      else:
         cmd = 'mjobctl -u user \"%s\"' % (' ').join(jobids)
         out = subprocess.check_output(cmd, shell=True)
         break
   return(out)

def submitdummy(home, paths, logger, jobids, queue='cbs'):
   '''Submit dummy job to wait for all other jobs to be finished'''
   
   print 'waiting for jobs to finish...'
   xmsub = '%sxmsub -l ncpus=1,mem=10mb,walltime=180,depend=%s -O %s/log/dummy -K -q %s -r y -t saqqaq_dummy.py' % (paths['pyscripts_home'], jobids[0], home, queue)
   logger.info(xmsub)
   try:
      lastjobid = subprocess.check_output(xmsub, shell=True)
   except:
      print 'Dummy job error, waiting 10m'
      time.sleep(600)
      lastjobid = subprocess.check_output(xmsub, shell=True)
   return lastjobid


##############################

###      SEMAPHORE       ###

##############################


def wait_semaphore(semaphore_ids, home, file_prefix, queue, check_interval, max_time):
   '''Wait for files to be created, times are in seconds'''
   
   from time import sleep
   import string
   import random
   import os
   import genobox_modules
   
   paths = genobox_modules.setSystem()
   
   # add directory and set semaphore filename
   if not os.path.exists('semaphores/'):
      os.makedirs('semaphores/')
   
   rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
   semaphore_file = 'semaphores/' + file_prefix + '.' + rand
   semaphore_file_err = 'log/' + file_prefix + '.' + rand + '.err'
   
   # submit job 
   depends = ':'.join(semaphore_ids)
   xmsub = '%sxmsub -d %s -l ncpus=1,mem=10mb,walltime=180,depend=%s -O %s -q %s -N semaphores -E %s -r y -t echo done' % (paths['pyscripts_home'], home, depends, semaphore_file, queue, semaphore_file_err)
   dummy_id = subprocess.check_output(xmsub, shell=True)
   
   # check for file to appear
   cnt = max_time
   while cnt > 0:
      if os.path.isfile(semaphore_file):
         break
      cnt -= check_interval
      sleep(check_interval)
   if cnt <= 0:
      raise SystemExit('%s did not finish in %is' % ())


##################################

###      FILE OPERATIONS       ###

##################################


def rm_files(patterns):
   '''Remove files using glob given as list of patterns'''
   
   import glob
   import os
   
   for p in patterns:
      files = glob.glob(p)
      if len(files) == 0:
         pass
      else:
         map(os.remove, files)






################################

###      LIBRARY FILES       ###

################################

# format is tab separated:
# ID, Data, Mapq, Read group information
# eg: ID Data	MAPQ	LB	PL	SM	CN

# test
#se = ['SRR002081se.recal.fastq', 'SRR002082se.recal.fastq']
#pe1 = ['SRR002137pe_1.recal.fastq', 'SRR002137pe_1.recal.fastq']
#pe2 = ['SRR002137pe_2.recal.fastq', 'SRR002137pe_2.recal.fastq']
#input = se + pe1 + pe2
#sample = 'NA12891'
#mapq = [30]
#libs = ['libSe']


## DEPRECATED ##

def read_library(library_file):
   '''Reads in library file and returns dict with ID as keys'''
   
   fh = open(library_file, 'r')
   header = fh.readline().rstrip().split('\t')
   
   # fill dict with data
   id_dict = dict()
   for line in fh:
      fields = line.rstrip().split('\t')
      if len(header) == len(fields):
         id_dict[fields[0]] = fields[0:]
      else:
         raise ValueError('Line \"%s\" do not have same number of fields as header')
   
   id_dict['header'] = header
   return id_dict

def library_from_input(input, sample='sample', mapq=[30], libs=['lib']):
   '''Create library dict from inputs'''
   
   import sys
   
   # check lengths of mapq and libs
   if len(mapq) < len(input):
      sys.stderr.write("warning: length of mapq is shorter than length of inputfiles, reusing %i\n" % mapq[0])   
      while len(mapq) < len(input):
         mapq.append(mapq[0])
   if len(libs) < len(input):
      sys.stderr.write("warning: length of libs is shorter than length of inputfiles, reusing %s\n" % libs[0])   
      while len(libs) < len(input):
         libs.append(libs[0])
   
   # create dict
   id_dict = dict()
   c = 0
   for i,f in enumerate(input):
      id = '%s_%i' % (sample, c + 1)
      id_dict[id] = [id, f, str(mapq[c]), libs[c], 'ILLUMINA', sample]
      c = c + 1
   
   id_dict['header'] = ['ID', 'Data', 'MAPQ', 'LB', 'PL', 'SM']
   
   # write to disk
   fh = open('libs.%s.txt' % sample, 'w')
   fh.write('%s\n' % '\t'.join(id_dict['header']))
   for key in sorted(id_dict.keys()):
      if key == 'header':
         continue
      fh.write('%s\n' % '\t'.join(id_dict[key]))
   return (id_dict, 'libs.%s.txt' % sample)

def update_libfile(libfile, key_col, new_col, val_dict, force=False):
   '''Update libfile with a new column of data'''
   
   import sys
   
   # read in current file
   fh = open(libfile, 'r')
   header = fh.readline().rstrip().split('\t')
   data = fh.read().split('\n')
   if data[-1] == '':
      data = data[:-1]
   fh.close()
   
   # check length of data and new val_dict
   if len(data) == len(val_dict):
      pass
   else:
      raise ValueError('Rows of libfile (%i) does not match length of new val_dict (%i)' % (len(data), len(val_dict)))
   
   # check if column already exist, if forced then remove old column
   if new_col in header:
      if force == False:
         sys.stderr.write('column %s already exists, %s not updated\n' % (new_col, libfile))
         return 
      else:
         n = header.index(new_col)
         header.remove(new_col)
         new_data = []
         for line in data:
            fields = line.split('\t')
            new_fields = fields[:n] + fields[(n+1):]
            new_data.append('\t'.join(new_fields))
         data = new_data
   
   # append data
   fh = open(libfile, 'w')
   header.append(new_col)
   fh.write('%s\n' % '\t'.join(header))
   for line in data:
      fields = line.split('\t')
      for key in val_dict.keys():
         if fields[header.index(key_col)] == key:
            new_line = '%s\t%s\n' % (line, val_dict[key])
            fh.write(new_line)
            break
   
   fh.close()
   return

def read_groups_from_libfile(index, libfile):
   '''Return read group from libfile with index as key'''
   
   allowedRGs = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
   
   RG = dict()
   #RG = '\@RG\tID:foo\tSM:bar'
   header = libfile['header']
   for key,value in libfile.items():
      if key == 'header':
         continue
      
      currRG = ['@RG']
      for h in header:
         if h in allowedRGs:
            currRG.append('%s:%s' % (h, value[header.index(h)]))
      RG[value[header.index(index)]] = currRG
   return RG

# read library file 
def read_bam_libs(libfile):
   '''Read library file and return dict'''
   
   from collections import defaultdict
   
   fh = open(libfile, 'r')
   header = fh.readline().rstrip().split('\t')
   
   # check if bamfiles are written in libfile
   if 'BAM' in header:
      pass
   else:
      raise IndexError('There is no column \"BAM\" in the libfile (%s)' % libfile)
   
   bam2lib = {}
   lib2bam = defaultdict(list)
   for line in fh:
      fields = line.rstrip().split('\t')
      bam2lib[fields[header.index('BAM')]] = (fields[header.index('MAPQ')], fields[header.index('LB')])
      
      lib2bam[fields[header.index('LB')]].append(fields[header.index('BAM')])
   
   # unique on lib2bam
   for key,values in lib2bam.items():
      lib2bam[key] = unique(values)
   
   return (bam2lib, lib2bam)

## DEPRECATED END ##



# Library inititation
def initialize_library(libfile, se=[], pe1=[], pe2=[], sample='sample', mapq=[30], libs=['A'], pl=['ILLUMINA'], bams=None):
   '''Initiates library file from arguments'''
   
   from genobox_classes import Library
   import random
   import string
   
   def try_append(index, from_list, target_list):
      '''Try to append value (indexed) from list to another list
         if the value does not exist reuse first value of list
         Converts all input values to strings
      '''
      try:
         target_list.append(str(from_list[index]))
      except:
         target_list.append(str(from_list[0]))
   
   if libfile:
      # copy library file so that it can be edited
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
      newlibfile = os.getcwd() + '/' + os.path.split(libfile)[1] + '.' + rand
      returnmsg = subprocess.check_call('cp %s %s' % (libfile, newlibfile), shell=True)
      libfile = newlibfile
      
      # create instance and read in library file (Library(libfile) ; .read())
      library = Library(libfile)
      library.read()
      
      # remove all non-input lines from library file
      library.keep('Data', se+pe1+pe2)      
   else:
      # else create new from input
      library = Library('libs.%s.txt' % sample)
      
      # check if sample is None
      if not sample:
         sample = 'sample'
      
      # create the library file
      f_count = 0
      (ID, Data, SM, MAPQ, LB, PL, BAM) = ([], [], [], [], [], [], [])
      if se and se != 'None':
         for i,f in enumerate(se):
            ID.append(sample + '_%i' % f_count)
            Data.append(f)
            SM.append(sample)
            try_append(f_count, mapq, MAPQ)
            try_append(f_count, libs, LB)
            try_append(f_count, pl, PL)
            f_count += 1
      
      if pe1 and pe1 != 'None':
         for i,f in enumerate(pe1):
            ID.append(sample + '_%i' % f_count)
            Data.append(f)
            SM.append(sample)
            try_append(f_count, mapq, MAPQ)
            try_append(f_count, libs, LB)
            try_append(f_count, pl, PL)
            f_count += 1
      
      if pe2 and pe2 != 'None':
         for i,f in enumerate(pe2):
            ID.append(sample + '_%i' % f_count)
            Data.append(f)
            SM.append(sample)
            try_append(f_count, mapq, MAPQ)
            try_append(f_count, libs, LB)
            try_append(f_count, pl, PL)
            f_count += 1
      
      if bams and bams != 'None':
         for i,f in enumerate(bams):
            ID.append(sample + '_%i' % f_count)
            Data.append(f)
            BAM.append(f)
            SM.append(sample)
            try_append(f_count, mapq, MAPQ)
            try_append(f_count, libs, LB)
            try_append(f_count, pl, PL)
            f_count += 1
      
      if bams and bams != 'None':
         library.create(ID=ID, Data=Data, SM=SM, MAPQ=MAPQ, LB=LB, PL=PL, BAM=BAM)
      else:
         library.create(ID=ID, Data=Data, SM=SM, MAPQ=MAPQ, LB=LB, PL=PL)
   
   return library



####################################

###      CHECK GENOME FILE       ###

####################################


def check_genome(genome):
   '''Check if genome file is valid'''
   
   # check if it can be opened
   try:
      fh = open(genome, 'r')
   except:
      sys.stdout.write('ERROR: Could not open genome-file\n')
      sys.exit(1)
   
   # check fields
   for line in fh:
      line = line.rstrip()
      fields = line.split('\t')
      for i,f in enumerate(fields):
         if i==1:
            try: 
               int(f)
            except: 
               sys.stdout.write('ERROR: Genome length in genome-file is not an integer\n')
               sys.exit(1)
         if i==3:
            if f == "haploid" or f=="diploid" or f=="na": pass
            else:
               sys.stdout.write('ERROR: Genome file ploidy is not set as diploid, haploid or na: %s\n' % f)
               sys.exit(1)
         if i==4:
            try:
               int(f)
            except:
               sys.stdout.write('ERROR: Genome file min coverage not an integer: %s\n' % f)
               sys.exit(1)
         if i ==5:
            try:
               int(f)
            except:
               sys.stdout.write('ERROR: Genome file max coverage not an integer: %s\n' % f)
               sys.exit(1)
   



##################################

###      FASTQ DETECTION       ###

##################################

def set_filetype(f, gz):
   '''Detects filetype from fa, fq and sff input file'''
   
   import gzip
   
   if gz:
      inhandle = gzip.open(f, 'rb')
   else:
      inhandle = open(f, "r")
   
   # parse
   line = inhandle.readline()
   if line.startswith(">"):
      out = 'fasta'
   elif line.startswith("@"):
      out = 'fastq'
   else:
      inhandle = open(f, "rb")
      line = inhandle.readline()
      if line.startswith(".sff"):
         out = 'sff'
      else:
         raise ValueError('Input must be fasta, fastq or sff')
   
   return out

def set_fqtype(f, gz):
   '''Detects sanger or illumina format from fastq'''
   
   # Open fastq, convert ASCII to number, check if number is above or below certain thresholds to determine format
   
   import gzip
   from Bio.SeqIO.QualityIO import FastqGeneralIterator
   
   if gz:
      inhandle = gzip.open(f, 'rb')
   else:
      inhandle = open(f, "r")
   
   count = 0
   type = 'not determined'
   for (title, sequence, quality) in FastqGeneralIterator(inhandle):
      qs = map(ord, quality)
      pos = None
      count +=1
      if title.startswith('ERR') or title.startswith('SRR'):
         for q in qs:
            if pos == 'Illumina_Solexa':
               if q < 64:
                  type = 'Solexa'
                  break
               if count >= 5000:
                  type = 'Illumina'
                  break
            if q > 83:
               pos = 'Illumina_Solexa'
            elif q < 59:
               type = 'Sanger'
               break
      else:
         for q in qs:
            if pos == 'Illumina_Solexa':
               if q < 64:
                  type = 'Solexa'
                  break
               if count >= 5000:
                  type = 'Illumina'
                  break
            if q > 83:
               pos = 'Illumina_Solexa'
            elif q < 59:
               type = 'Sanger'
               break
      if type != 'not determined':
         break
   
   if type == 'not determined' and count < 5000 and pos == 'Illumina_Solexa':
      type = 'Illumina'
   elif type == 'not determined':
      sys.stderr.write('Fastq format not identified, are you sure it is sanger/illumina/solexa?')
      return None
   return type


###################################

###      GENERAL FUNCTIONS      ###

###################################


def unique(seq):
   '''Return unique and ordered list'''
   
   # order preserving 
   checked = [] 
   for e in seq: 
      if e not in checked: 
         checked.append(e) 
   return checked


def rm_files(patterns):
   '''Remove files using glob given as list of patterns'''
   
   import glob
   import os
   
   for p in patterns:
      files = glob.glob(p)
      if len(files) == 0:
         pass
      else:
         map(os.remove, files)

