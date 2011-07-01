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
   
   def __init__(self, calls, logfile=None, runname='run_test', queue='cbs', cpu='nodes=1:ppn=1,mem=2gb,walltime=43200', depend=False, depend_type='one2one', depend_val=[], hold=True, depend_ids=[]):
      '''Constructor for Moab class'''
      self.calls = calls
      self.runname = runname
      self.logfile = logfile
      self.queue = queue
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
      msg = 'Moab(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' % ("["+", ".join(self.calls)+"]", self.logfile, self.runname, self.queue,  self.cpu, str(self.depend), self.depend_type, "["+", ".join(map(str,self.depend_val))+"]", str(self.hold), "["+", ".join(self.depend_ids)+"]")
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
            xmsub = cmd+' -d %s -l %s -O %s -E %s -r y -q %s -N %s -t %s' % (home, self.cpu, stdout, stderr, self.queue, self.runname, call)
         else:
            xmsub = cmd+' -d %s -l %s,depend=%s -O %s -E %s -r y -q %s -N %s -t %s' % (home, self.cpu, depends[i], stdout, stderr, self.queue, self.runname, call)
         
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
         #fh.write('newgrp cge\n')      # temporary fix to group problems on cge-s2 cluster
         fh.write('%s\n' % call)
         fh.close()
         
         # create msub command
         cmd = 'msub'
         if self.hold:
            cmd = '%s -h' % cmd
         if not self.depend:
            msub = '%s -d %s -l %s -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, self.cpu, stdout, stderr, self.queue, self.runname, filename)
         else:
            msub = '%s -d %s -l %s,depend=%s -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, self.cpu, depends[i], stdout, stderr, self.queue, self.runname, filename)
         
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

class Library:
   '''Class for use, construct and read from library file'''
   
   def __init__(self, f):
      '''Constructor for lib class'''
      self.f = f
      self.data = {}
      self.allowedRGs = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
   
   def __repr__(self):
      '''Return string of attributes'''
      msg = 'Lib(%s)' % (self.f)
      return msg
   
   def read(self):
      '''Reads in library file and returns dict with tags as keys and data as values'''
      
      fh = open(self.f, 'r')
      header = fh.readline().rstrip().split('\t')
      
      # create dict with tags as key
      data = {}
      for key in header:
         data[key] = []
      
      # fill dict with data
      for line in fh:
         fields = line.rstrip().split('\t')
         for i in range(len(fields)):
            data[header[i]] = data[header[i]]+[fields[i]]
      
      self.data = data
   
   def create(self, **kwargs):
      '''Create library dict from inputs and write to file (self.f)
         input: list of tagA=[value1, value2], tagB=[value1, value2] ..., tag names will become the header
      '''
      
      # open handle
      fh = open(self.f, 'w')
      
      header = kwargs.keys()
      header.sort()
      fh.write('\t'.join(header)+'\n')
      for i in range(len(kwargs[kwargs.keys()[0]])):
         for key in header:
            fh.write(kwargs[key][i])
            if key == header[-1]:
               fh.write('\n')
            else:
               fh.write('\t')
      fh.close()
      
      # when created read into self
      self.read()
   
   def write(self, data):
      '''Writes library dictionary to file'''
      
      # open handle
      fh = open(self.f, 'w')
      
      header = data.keys()
      header.sort()
      fh.write('\t'.join(header)+'\n')
      for i in range(len(data[data.keys()[0]])):
         for key in header:
            fh.write(data[key][i])
            if key == header[-1]:
               fh.write('\n')
            else:
               fh.write('\t')
      fh.close()
   
   def update(self, **kwargs):
      '''Add column(s) to library file'''
      
      # check if new column is same length
      for key in kwargs.keys():
         try:
            if len(kwargs[key]) != len(self.data[self.data.keys()[0]]):
               raise ValueError('New column is not of same length as existing')
         except:
            pass
      
      # read in dict, update with new column and write to file
      print self.data
      print "-------------"
      current = self.data
      current.update(kwargs)
      print current
      self.write(current)
      self.data = current
   
   def update_with_tag(self, key_col, new_col, val_dict, force=True):
      '''Add column to library file using keys in dictionary as tag in key_col'''
      
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
      self.read()
   
   def keep(self, tag, values):
      '''Remove rows that does not contain any of the tag=values'''
      
      rowsToRemove = []
      if self.data.has_key(tag):
         for i in range(len(self.data[self.data.keys()[0]])):
            if self.data[tag][i] in values:
               pass
            else:
               rowsToRemove.append(i)
         # remove rows
         for key in self.data.keys():
            for element in rowsToRemove[::-1]:
               del self.data[key][element]
      
      self.write(self.data)
      
   def getRG(self, tag):
      '''Return read group from data with tag as key'''
      
      RG = {}
      header = self.data.keys()      
      for i in range(len(self.data[self.data.keys()[0]])):
         currRG = ['@RG']
         for key in self.data.keys():
            if key in self.allowedRGs:
               currRG.append('%s:%s' % (key, self.data[key][i]))              
         RG[self.data[tag][i]] = currRG
      return RG
   
   def getPL(self, tag):
      '''Return platform from data with tag as key'''
      
      from collections import defaultdict
      
      if 'PL' in self.data.keys():
         pass
      else:
         raise IndexError('There is no column \"PL\" in the libfile (%s)' % self.f)
      
      PL = {}
      for i in range(len(self.data[self.data.keys()[0]])):
         PL[self.data['Data'][i]] = self.data['PL'][i]
      
      PL2data = defaultdict(list)
      for i in range(len(self.data[self.data.keys()[0]])):
         PL2data[self.data['PL'][i]].append(self.data['Data'][i])
      
      return (PL, PL2data)
   
   def getBamLibs(self):
      '''Return two dictionaries:
         key=BAM, value=(MAPQ, LB)
         key=LB, value=[BAM1, BAM2,...]
      '''
      
      from collections import defaultdict
      
      def unique(seq):
         '''Return unique and ordered list'''
         
         # order preserving 
         checked = [] 
         for e in seq: 
            if e not in checked: 
               checked.append(e) 
         return checked
      
      # check if BAM-files are written
      if 'BAM' in self.data.keys():
         pass
      else:
         raise IndexError('There is no column \"BAM\" in the libfile (%s)' % self.f)
      
      bam2lib = {}
      for i in range(len(self.data[self.data.keys()[0]])):
         bam2lib[self.data['BAM'][i]] = (self.data['MAPQ'][i], self.data['LB'][i])
      
      lib2bam = defaultdict(list)
      for i in range(len(self.data[self.data.keys()[0]])):
         lib2bam[self.data['LB'][i]].append(self.data['BAM'][i])
      
      # unique on lib2bam
      for key,values in lib2bam.items():
         lib2bam[key] = unique(values)
      
      return (bam2lib, lib2bam)

