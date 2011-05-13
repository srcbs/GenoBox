#!/usr/bin/python

## For testing ##

# create calls
#align_calls = ['bwa aln -t 4 input1.file reference.fasta > output1.sai', 'bwa aln -t 4 input2.file reference.fasta > output2.sai']
#samse_calls = ['bwa samse output1.sai reference.fasta input1.file > input1.aln.bam', 'bwa samse output2.sai reference.fasta input2.file > input2.aln.bam']

# create moab instance for the align_calls and dispatch to queue
#align_moab = Moab(align_calls, runname='run_align', cpu='nodes=1:ppn=4,mem=8gb,walltime=43200')
#align_ids = align_moab.dispatch()
#samse_moab = Moab(samse_calls, runname='run_samse', cpu='nodes=1:ppn=1,mem=2gb,walltime=43200', depend=True, depend_type='one2one', depend_val = [1], align_ids)
#samse_ids = samse_moab.dispatch()

# release jobs
#align_moab.releasejobs()
#samse_moab.releasejobs()


# class
class Moab:
   '''Submits a list of calls to the scheduler using msub/xmsub. Job dependencies are controlled using depend (logical), depend_type ('one2one', 'expand', 'conc', 'complex'), depend_val (list of integers) and ids (list of jobids):
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes several integers in a list as input and makes the jobs dependent on the number of ids given. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
      
      Jobs are submitted as hold by default and should be released using Moab.releasejobs().
   '''
   
   def __init__(self, calls, logfile=None, runname='run_test', queue='cbs', group='cdrom', cpu='nodes=1:ppn=1,mem=2gb,walltime=43200', depend=False, depend_type='one2one', depend_val=[], hold=True, ids=[]):
      '''Constructor for Moab class'''
      self.calls = calls
      self.runname = runname
      self.logfile = logfile
      self.queue = queue
      self.group = group
      self.cpu = cpu
      self.depend = depend
      self.depend_type = depend_type
      self.depend_val = depend_val
      self.hold = hold
      self.ids = ids
   
   def __repr__(self):
      '''Return string of attributes'''
      msg = 'Moab(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' % ("["+", ".join(self.calls)+"]", self.logfile, self.runname, self.queue, self.group, self.cpu, str(self.depend), self.depend_type, "["+", ".join(map(str,self.depend_val))+"]", str(self.hold), "["+", ".join(self.ids)+"]")
      return msg
   
   def get_logger(self):
      '''Return logger object'''
      
      if self.logfile:
         import logging
         
         logger = logging.getLogger('moab_submissions')
         hdlr = logging.FileHandler('%s' % self.logfile)
         formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
         hdlr.setFormatter(formatter)
         logger.addHandler(hdlr) 
         logger.setLevel(logging.INFO)
         
         return logger
      else:
         return None
   
   def create_dependencies(self):
      '''Create job dependencies
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
      '''
      
      if not self.depend:
         depends = None
      else:
            if self.depend_type == 'one2one':
               depends = self.ids
            
            elif self.depend_type == 'expand':
               n = int(self.depend_val[0])
               depends = []
               c = 0
               for j in range(len(self.calls)):
                  c = c + 1
                  if c < int(n):
                     depends.append(ids[0])
                  if c == int(n):
                     depends.append(ids.pop(0))
                     c = 0
            
            elif self.depend_type == 'conc':
               n = int(self.depend_val[0])
               depends = []
               for j in range(len(self.calls)):
                  s = ':'.join(a[:int(n)])
                  depends.append(s)
                  a = a[int(n):]
            
            elif self.depend_type == 'complex':
               old_index = 0
               depends = []
               for index in self.depend_val:
                  s = ':'.join(a[old_index:(index+old_index)])
                  depends.append(s)
                  old_index=index
            else:
               raise AttributeError('depend_type not recognized: %s' % self.depend_type)
      return depends
   
   def submit_xmsub(self, depends, logger):
      '''Submits jobs using xmsub'''
      
      import re
      import subprocess
      import time
      import os
      
      home = os.getcwd()
      
      ids = []
      for i in range(len(self.calls)):
         call = self.calls[i]
         stdout = '%s/log/%s%i.o' % (home, self.runname, i)
         stderr = '%s/log/%s%i.e' % (home, self.runname, i)
         
         # catch stdouts if call includes 'program infile > outfile', needs to be directed as -O instead of >
         pattern = re.compile(r'(^.+)>\s(.+)$')
         match = pattern.search(call)
         if match:
            call = match.group(1)
            stdout = '%s/%s' % (home, match.group(2))
         
         # create xmsub commands
         cmd = '/panvol1/simon/bin/pipeline/xmsub'
         
         # toggle if job should be on hold
         if self.hold:
            cmd = '%s -h ' % cmd
         if not self.depend:
            xmsub = cmd+' -d %s -l %s -O %s -E %s -r y -q %s -N %s -W group_list=%s -t %s' % (home, self.cpu, stdout, stderr, self.queue, self.runname, self.group, call)
         else:
            xmsub = cmd+' -d %s -l %s,depend=%s -O %s -E %s -r y -q %s -N %s -W group_list=%s -t %s' % (home, self.cpu, depends[i], stdout, stderr, self.queue, self.runname, self.group, call)
         
         time.sleep(0.1)
         if logger:
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
   
   def submit_wrapcmd(self, depends, logger):
      '''Take input as command and submit to msub (this way pipes and redirects can be done)'''
      
      import subprocess
      import random
      import string
      import time
      import os
      
      home = os.getcwd()
            
      ids = []
      for i in range(len(self.calls)):
         call = self.calls[i]
         N = 10
         rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))
         filename = 'pbsjob.tmp%s' % rand
         
         stdout = '%s/log/%s%i.o' % (home, self.runname, i)
         stderr = '%s/log/%s%i.e' % (home, self.runname, i)
         
         # write pbsjob file
         fh = open(filename, 'w')
         fh.write('#!/bin/sh\n\n')
         fh.write('newgrp cge\n')      # temporary fix to group problems on cge-s2 cluster
         fh.write('%s\n' % call)
         fh.close()
         
         # create msub command
         cmd = 'msub'
         if self.hold:
            cmd = '%s -h' % cmd
         if not self.depend:
            msub = '%s -d %s -l %s -o %s -e %s -q %s -r y -N %s -W group_list=%s %s' % (cmd, home, self.cpu, stdout, stderr, self.queue, self.runname, self.group, filename)
         else:
            msub = '%s -d %s -l %s,depend=%s -o %s -e %s -q %s -r y -N %s -W group_list=%s %s' % (cmd, home, self.cpu, depends[i], stdout, stderr, self.queue, self.runname, self.group, filename)
         
         time.sleep(0.1)
         
         if logger:
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
     
   def releasejobs(self):
      '''Release submitted jobs from hold'''
      
      import subprocess
      
      print  "Releasing jobs"
      while len(self.ids) > 0:
         if len(self.ids) > 199:
            cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids[:199])
            del self.ids[:199]
            out = subprocess.check_output(cmd, shell=True)
         else:
            cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids)
            out = subprocess.check_output(cmd, shell=True)
            break
   
   def dispatch(self):
      '''Submit job to queue, if depend=True dependent jobids should be given as ids. 
         Dependtype can be 'one2one', 'expand' 'conc', 'complex'.
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...'''
      
      # start logger
      logger = self.get_logger()
      
      # create job dependencies
      depends = self.create_dependencies()
      
      if self.calls[0].find('|') > -1 or self.calls[0].find('<') > -1 or self.calls[0].find('>>') > -1:
         # perform wrapcmd if calls includes pipes / left-redirects
         self.ids = self.submit_wrapcmd(depends, logger)
      else:
         # perform xmsub if calls does not include pipes (can have right-redirects)
         self.ids = self.submit_xmsub(depends, logger)
      
      return self.ids
