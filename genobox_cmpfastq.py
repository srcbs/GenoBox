#!/panvol1/simon/bin/python2.7

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import sys
import random
import string
import cdb
import multiprocessing
import gzip

class Cmpfastq():
   '''Compare reads from two paired end read files and output common and unique reads'''
   
   def __init__(self, i, gz):
      '''Constructor'''
      
      self.i = i
      self.gz = gz
      
      # set output files
      o_common = []
      o_unique = []
      for f in self.i:
         if gz: 
            o_common.append(f+'.common.fq.gz')
            o_unique.append(f+'.unique.fq.gz')
         else: 
            o_common.append(f+'.common.fq')
            o_unique.append(f+'.unique.fq')
      self.oc = o_common
      self.ou = o_unique
   
   def __repr__(self):
      msg = 'Cmpfastq(%s, %s)' % ("["+", ".join(self.i)+"]", self.gz)
      return msg
   
   def cmp(self):
      '''Compare two paired fastq files'''
      
      def intersect(a, b):
         '''Intesection between lists'''
         return list(set(a) & set(b))
         
      def write_out(db_common, f, o, gz):
         '''Write out reads'''
         
         if gz:
            fh = gzip.open(f, 'rb')
            out = gzip.open(o, 'wb')
         else:
            fh = open(f, 'r')
            out = open(o, 'w')
         
         written_count = 0
         total_count = 0
         for (title, sequence, quality) in FastqGeneralIterator(fh):
            total_count += 1
            if db_common.has_key(title[:-2]):
               out.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
               written_count += 1
         sys.stderr.write('%s: Total %i, Written %i (%.1f%%)\n' % (f, total_count, written_count, written_count/total_count*100))
         fh.close()
         out.close()
      
      def create_db(f, db_fname, gz):
         '''Write out db of headers'''
         
         if gz:
            fh = gzip.open(f, 'rb')
         else:
            fh = open(f, 'r')
         
         fh_headers = (x.strip()[1:-2] for i, x in enumerate(fh) if not (i % 4))
         
         db = cdb.cdbmake(db_fname, db_fname + '.tmp')
         for h in fh_headers:
            db.add(h, 'T')
         db.finish()
         del(db)
      
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
      
      ## get headers from both trimmed files ##
      # strip the /2 or /1 and grab only the headers
      # write in dbm to minimze memory usage
      
      # create db's (parallel)
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(36))
      db1_fname = 'db1_%s' % rand
      db2_fname = 'db2_%s' % rand
      
      jobs = []
      p = multiprocessing.Process(target=create_db, args=(self.i[0], db1_fname, self.gz, ))
      p.start()
      jobs.append(p)
      
      p = multiprocessing.Process(target=create_db, args=(self.i[1], db2_fname, self.gz, ))
      p.start()
      jobs.append(p)
      
      # wait for jobs to finish
      for job in jobs:
         job.join()
      
      ## get headers that are in both trimmed files ##
      db1 = cdb.init(db1_fname)
      db2 = cdb.init(db2_fname)
      common = intersect(db1.keys(), db2.keys())
      
      dbcommon_fname = 'dbcommon_%s' % rand
      db_common = cdb.cdbmake(dbcommon_fname, dbcommon_fname + '.tmp')
      for h in common:
         db_common.add(h, 'T')
      db_common.finish()
      del(db_common)
      
      ## get headers that are in only one trimmed file ##
      symdiff = set(db1.keys()).symmetric_difference(set(db2.keys()))
      
      dbdiff_fname = 'dbdiff_%s' % rand
      db_diff = cdb.cdbmake(dbdiff_fname, dbdiff_fname + '.tmp')
      for h in symdiff:
         db_diff.add(h, 'T')
      db_diff.finish()
      del(db_diff)
      
      
      ## open common db ##
      db_common = cdb.init(dbcommon_fname)
      jobs = []
      p = multiprocessing.Process(target=write_out, args=(db_common, self.i[0], self.oc[0], self.gz,))
      p.start()
      jobs.append(p)
      
      p = multiprocessing.Process(target=write_out, args=(db_common, self.i[1], self.oc[1], self.gz,))
      p.start()
      jobs.append(p)
      
      ## open single db ##
      
      db_diff = cdb.init(dbdiff_fname)
      p = multiprocessing.Process(target=write_out, args=(db_diff, self.i[0], self.ou[0], self.gz,))
      p.start()
      jobs.append(p)
      
      p = multiprocessing.Process(target=write_out, args=(db_diff, self.i[1], self.ou[1], self.gz,))
      p.start()
      jobs.append(p)
      
      # wait for jobs to finish
      for job in jobs:
         job.join()
      
      rm_files([db1_fname, db2_fname, dbcommon_fname, dbdiff_fname])
         
      
      
      

if __name__ == '__main__':
      
   parser = argparse.ArgumentParser(prog='genobox_cmpfastq.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=30,width=100),
                                 description='''Compare two fastq-files, output the common and unique reads in different files''', 
                                 usage='%(prog)s [options]')
   
   # add the arguments
   parser.add_argument('--i', help='input paired end files', nargs=2, required=True)
   parser.add_argument('--gz', help='input files are gzipped [False]', default=False, action='store_true')
   parser.add_argument('--log', help='log level [INFO]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--i tmp_1.fq.gz tmp_2.fq.gz --gz'.split())
   
   c = Cmpfastq(args.i, args.gz)
   c.cmp()
