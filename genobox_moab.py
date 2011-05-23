#!/usr/bin/python

import platform
import os
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
      sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.16/'
      sys['picard_home'] = '/panvol1/simon/bin/picard-tools-1.26/'
      sys['bedtools_home'] = '/panvol1/simon/bin/BEDTools-2.8.3/'
      sys['greengenes_home'] = '/panvol1/simon/databases/greengenes/'
      sys['mothur_home'] = '/panvol1/simon/bin/Mothur-1.11.0/'
      sys['rna_hmm3_home'] = '/panvol1/simon/bin/rna_hmm3/'
      sys['hmm-3.0_home'] = '/panvol1/simon/bin/hmmer-3.0/'
      sys['phymm_home'] = '/panvol1/simon/bin/phymm/'
      sys['prodigal_home'] = '/panvol1/simon/bin/prodigal/'
      sys['rnammer_home'] = '/panvol1/simon/bin/rnammer-1.2/'
      sys['cge_genomes'] = '/panvol1/simon/databases/cge_genomes/'
      sys['velvet_home'] = '/panvol1/simon/bin/velvet_1.1.02/'
      sys['newbler'] = '/panfs/saqqaq/bin/newbler/bin/'
      sys['stampy'] = '/panvol1/simon/bin/stampy-1.0.6/'
      sys['chain_index'] = '/panvol1/simon/databases/chainmap/'
      sys['taxonomy_ncbi'] = '/panvol1/simon/databases/taxonomy/'
      sys['R_home'] = '/tools/bin/'
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
   
   group = 'cdrom'
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
         xmsub = cmd+' -d %s -l %s,walltime=172800 -O %s -E %s -r y -q %s -N %s -W group_list=%s -t %s' % (home, cpu, stdout, stderr, queue, runname, group, call)
      else:
         xmsub = cmd+' -d %s -l %s,walltime=172800,depend=%s -O %s -E %s -r y -q %s -N %s -W group_list=%s -t %s' % (home, cpu, depends[i], stdout, stderr, queue, runname, group, call)
      
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
   
   group = 'cdrom'
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
         msub = '%s -d %s -l %s,walltime=172800 -o %s -e %s -q %s -r y -N %s -W group_list=%s %s' % (cmd, home, cpu, stdout, stderr, queue, runname, group, filename)
      else:
         msub = '%s -d %s -l %s,walltime=172800,depend=%s -o %s -e %s -q %s -r y -N %s -W group_list=%s %s' % (cmd, home, cpu, depends[i], stdout, stderr, queue, runname, group, filename)
      
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




