#!/usr/bin/python


# class
class Moab:
   '''Submits a list of calls to the scheduler using msub/xmsub. Job dependencies are controlled using depend (logical), depend_type ('one2one', 'expand', 'conc', 'complex'), depend_val (list of integers) and ids (list of jobids):
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes several integers in a list as input and makes the jobs dependent on the number of ids given. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
      
      Jobs are submitted as hold by default and should be released using Moab.release().
   '''
   
   def __init__(self, calls, logfile=None, runname='run_test', queue='cbs', group='cdrom', cpu='nodes=1:ppn=1,mem=2gb,walltime=43200', depend=False, depend_type='one2one', depend_val=[], hold=True, depend_ids=[]):
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
      self.depend_ids = depend_ids
      
      # put jobs in queue upon construction
      self.dispatch()
   
   def __repr__(self):
      '''Return string of attributes'''
      msg = 'Moab(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' % ("["+", ".join(self.calls)+"]", self.logfile, self.runname, self.queue, self.group, self.cpu, str(self.depend), self.depend_type, "["+", ".join(map(str,self.depend_val))+"]", str(self.hold), "["+", ".join(self.depend_ids)+"]")
      return msg
   
   def get_logger(self):
      '''Return logger object'''
      
      if self.logfile:
         import logging
         
         # if already a logger do nothing and return
         # else create logger with given logfile
         if isinstance(self.logfile, logging.Logger):
            return self.logfile
         else:
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
               depends = self.depend_ids
            
            elif self.depend_type == 'expand':
               n = int(self.depend_val[0])
               depends = []
               depend_ids = self.depend_ids     # do not want to remove from self.depend_ids because it will remove the ids from input class instance (eg. other calls)
               c = 0
               for j in range(len(self.calls)):
                  c = c + 1
                  if c < int(n):
                     depends.append(depend_ids[0])
                  if c == int(n):
                     depends.append(depend_ids[0])
                     depend_ids = depend_ids[1:]
                     c = 0
            
            elif self.depend_type == 'conc':
               n = int(self.depend_val[0])
               depends = []
               for j in range(len(self.calls)):
                  s = ':'.join(self.depend_ids[:int(n)])
                  depends.append(s)
                  self.depend_ids = self.depend_ids[int(n):]
            
            elif self.depend_type == 'complex':
               old_index = 0
               depends = []
               for index in self.depend_val:
                  s = ':'.join(self.depend_ids[old_index:(index+old_index)])
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
            #print id
         except:
            print 'Job error, waiting 1m'
            time.sleep(60)
            id = subprocess.check_output(xmsub, shell=True)
         ids.append(id.split('\n')[1])
         #print ids
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
            #print id
         except:
            print 'Job error, waiting 1m'
            time.sleep(60)
            id = subprocess.check_output(msub, shell=True)
         ids.append(id.split('\n')[1])
         
         # remove pbsjob file
         #rm_files([filename])
         #print ids
      return ids
   
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
      
   
   def release(self):
      '''Release submitted jobs from hold'''
      
      import subprocess
      
      while len(self.ids) > 0:
         if len(self.ids) > 199:
            cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids[:199])
            del self.ids[:199]
            out = subprocess.check_output(cmd, shell=True)
         else:
            cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids)
            out = subprocess.check_output(cmd, shell=True)
            break
   


class Semaphore:
   '''Wait for files to be created, times are in seconds'''
   
   def __init__(self, semaphore_ids, home, file_prefix, queue, check_interval, max_time):
      '''Constructor for Semaphore class'''
      self.semaphore_ids = semaphore_ids
      self.home = home
      self.file_prefix = file_prefix
      self.queue = queue
      self.check_interval = check_interval
      self.max_time = max_time
   
   def wait(self):
      '''Wait for files to be created'''
      
      from time import sleep
      import string
      import random
      import os
      import genobox_modules
      import subprocess
      
      paths = genobox_modules.setSystem()
      
      # add directory and set semaphore filename
      if not os.path.exists('semaphores/'):
         os.makedirs('semaphores/')
      
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
      semaphore_file = 'semaphores/' + self.file_prefix + '.' + rand
      semaphore_file_err = 'log/' + self.file_prefix + '.' + rand + '.err'
      
      # submit job 
      depends = ':'.join(self.semaphore_ids)
      xmsub = '%sxmsub -d %s -l ncpus=1,mem=10mb,walltime=180,depend=%s -O %s -q %s -N semaphores -E %s -r y -t echo done' % (paths['pyscripts_home'], self.home, depends, semaphore_file, self.queue, semaphore_file_err)
      dummy_id = subprocess.check_output(xmsub, shell=True)
      
      # check for file to appear
      cnt = self.max_time
      while cnt > 0:
         if os.path.isfile(semaphore_file):
            break
         cnt -= self.check_interval
         sleep(self.check_interval)
      if cnt <= 0:
         raise SystemExit('%s did not finish in %is' % ())

class Lib:
   '''class for use, construct and read from library file'''
   
   def __init__(self, f, mapq, libs, pl, sample):
      '''Constructor for lib class'''
      self.f = f
      self.mapq = mapq
      self.libs = libs
      self.pl = pl
      self.sample = sample
   
   def __repr__(self):
      '''Return string of attributes'''
      msg = 'Lib(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' % (f, "["+", ".join(self.mapq)+"]", "["+", ".join(self.libs)+"]", "["+", ".join(self.pl)+"]", self.sample)
      return msg
   
   def read(self):
      '''Reads in library file and returns dict with ID as keys'''
      
      fh = open(self.f, 'r')
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
   
   def create(self):
      '''Create library dict from inputs'''
      
      import sys
      
      # check lengths of mapq and libs
      if len(self.mapq) < len(self.input):
         sys.stderr.write("warning: length of mapq is shorter than length of inputfiles, reusing %i\n" % self.mapq[0])   
         while len(self.mapq) < len(self.input):
            self.mapq.append(self.mapq[0])
      if len(self.libs) < len(self.input):
         sys.stderr.write("warning: length of libs is shorter than length of inputfiles, reusing %s\n" % self.libs[0])   
         while len(self.libs) < len(self.input):
            self.libs.append(self.libs[0])
      
      # create dict
      id_dict = dict()
      c = 0
      for i,input in enumerate(input):
         id = '%s_%i' % (self.sample, c + 1)
         id_dict[id] = [id, input, str(self.mapq[c]), self.libs[c], self.pl, self.sample]
         c = c + 1
      
      id_dict['header'] = ['ID', 'Data', 'MAPQ', 'LB', 'PL', 'SM']
      
      # write to disk
      fh = open('libs.%s.txt' % self.sample, 'w')
      fh.write('%s\n' % '\t'.join(id_dict['header']))
      for key in sorted(id_dict.keys()):
         if key == 'header':
            continue
         fh.write('%s\n' % '\t'.join(id_dict[key]))
      return (id_dict, 'libs.%s.txt' % self.sample)
   
   
   def update(self, key_col, new_col, val_dict, force=False):
      '''Update libfile with a new column of data'''
      
      import sys
      
      # read in current file
      fh = open(self.f, 'r')
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
            sys.stderr.write('column %s already exists, %s not updated\n' % (new_col, self.f))
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
      fh = open(self.f, 'w')
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
   
   def getRG(self, index):
      '''Return read group from lib_dict with index as key'''
      
      allowedRGs = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
      
      id_dict = self.read(self.f)
      
      RG = dict()
      #RG = '\@RG\tID:foo\tSM:bar'
      header = self.f['header']
      for key,value in id_dict.items():
         if key == 'header':
            continue
         
         currRG = ['@RG']
         for h in header:
            if h in allowedRGs:
               currRG.append('%s:%s' % (h, value[header.index(h)]))
         RG[value[header.index(index)]] = currRG
      return RG
   
   def get_bamlibs(self):
      '''Read library file and return dict'''
      
      from collections import defaultdict
      
      fh = open(self.f, 'r')
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
