#!/panvol1/simon/bin/python2.7

from __future__ import division
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import genobox_modules
import sys
import random
import string
import cdb

class FastqTrim:
   '''Trim/Filter paired end fastq:
      trim 3' trailing bases under (q)
      trim sequencing adaptors (Illumina) (a)
      filter reads with less than avg. Q quals (m)
      filter reads with less than L length (l)
      filter reads with Ns (keep_n)
   '''
   
   def __init__(self, f, o, l=25, q=20, m=20, keep_n=False, M=20, 
                a=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 
                'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 
                'ACACTCTTTCCCTACACGACGCTCTTCCGATCT']):
      self.f = f
      self.o = o
      self.l = l
      self.q = q
      self.m = m
      self.keep_n = keep_n
      self.M = M
      self.a = a
      
      # set adaptors and fastq format
      self.adaptors = self.set_adaptors()
      self.format = genobox_modules.set_filetype(self.f[0])
      if self.format == 'fastq':
         self.fqtype = genobox_modules.set_fqtype(self.f[0])
   
   
   def __repr__(self):
      msg = 'FastqTrim(%s, %s, %i, %i, %i, %s, %i)' % ("["+", ".join(self.f)+"]", "["+", ".join(self.o)+"]", self.l, self.q, self.m, self.keep_n, self.M)
      return msg
   
   def set_adaptors(self):
      '''Set adaptors according to adaptor length requirements'''
      
      adaptors = []
      for a in self.a:
         if self.M == 0:
            adaptor = a
         elif (self.M) > len(a):
            adaptor = a
         else:
            adaptor = a[:self.M]
         adaptors.append(adaptor)
      return adaptors
   
   def filter_adaptor(self, title, sequence, quality):
      '''Filters adaptor and Ns in read'''
      
      has_adaptor = False
      has_N = False
      is_short = False
      
      # search for adaptor in reads
      for adaptor in self.adaptors:
         hit = sequence.find(adaptor)
         if hit > -1:
            sequence = sequence[:hit]
            quality = quality[:hit]
            has_adaptor = True
      
      # check if read contains Ns
      if sequence.find('N') > -1:
         has_N = True
      
      # check if too small then do not print out
      if len(sequence) < self.l:
         is_short = True
      
      # return read
      if has_N or is_short:
         return (None, None, None)
      else:
         return (title, sequence, quality)
      
   # the actual read trimming function
   def trim_qual(self, title, sequence, quality):
      '''Trims read based on 3' qualities and average quality'''
      
      # functions for read trimming
      def average(values):
         '''Computes the arithmetic mean of a list of numbers.'''
         return sum(values, 0.0) / len(values)
      
      def illumina2qual(qual_string):
         '''Convert illumina qualities (offset 64) to values.'''
         qual = []
         for q in qual_string:
            qual.append(ord(q)-64)
         return qual
      
      def sanger2qual(qual_string):
         '''Convert sanger qualities (offset 33) to values.'''
         qual = []
         for q in qual_string:
            qual.append(ord(q)-33)
         return qual
      
      def trim_index(qual, m):
         '''Determine position to trim from'''
         qs = qual[::-1]
         index = None
         for i,q in enumerate(qs):
            if q > m:
               index = len(qual)-i
               break
         return index
      
      # check if read is not valid
      if title == None:
         return (None, None, None)
      
      if self.fqtype == 'Illumina':
         qual = illumina2qual(quality)
      elif self.fqtype == 'Sanger':
         qual = sanger2qual(quality)
      else:
         raise ValueError('Fastq must be Illumina or Sanger')
      
      # check trailing qualities, should it be trimmed?
      t = trim_index(qual, self.q)
      if t:
         sequence = sequence[:t]
         quality = quality[:t]
         # do not print out if average quality is less than m or length is less than l
         a = average(qual[:t])
         if a < self.m or len(sequence) < self.l:
            return (None, None, None)
         else:
            return (title, sequence, quality)
      else:
         # do not print out if average quality is less than m or length is less than l
         a = average(qual)
         if a < self.m or len(sequence) < self.l:
            return (None, None, None)
         else:
            return (title, sequence, quality)
   
   def write_pairs(self, f1, f2):
      '''Parse through two paired files and only write if both pairs are present'''
      
      def intersect(a, b):
         '''Intesection between lists'''
         return list(set(a) & set(b))
      
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
      
      # file 1
      fh1 = open(f1, 'r')
      fh1_headers = (x.strip()[1:-2] for i, x in enumerate(fh1) if not (i % 4))
      
      # create db
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(36))
      db1_fname = 'db1_%s' % rand
      db1 = cdb.cdbmake(db1_fname, db1_fname + '.tmp')
      for h in fh1_headers:
         db1.add(h, 'T')
      db1.finish()
      del(db1)
      
      # file 2
      fh2 = open(f2, 'r')
      fh2_headers = (x.strip()[1:-2] for i, x in enumerate(fh2) if not (i % 4))
      db2_fname = 'db2_%s' % rand
      db2 = cdb.cdbmake(db2_fname, db2_fname + '.tmp')
      for h in fh2_headers:
         db2.add(h, 'T')
      db2.finish()
      del(db2)
      
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
      
      # parse through each file and write if read is in db_common
      out1 = open(self.o[0], 'w')
      out2 = open(self.o[1], 'w')
      fh1 = open(f1, 'r')
      fh2 = open(f2, 'r')
      
      # open common db
      db_common = cdb.init(dbcommon_fname)
      # write out for read 1
      written_count = 0
      total_count = 0
      for (title, sequence, quality) in FastqGeneralIterator(fh1):
         total_count += 1
         if db_common.has_key(title[:-2]):
            out1.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
            written_count += 1
      sys.stderr.write('%s: Total %i, Written %i (%.1f%%)\n' % (self.f[0], total_count, written_count, written_count/total_count*100))
      
      # write out for read 2
      written_count = 0
      total_count = 0
      for (title, sequence, quality) in FastqGeneralIterator(fh2):
         total_count += 1
         if db_common.has_key(title[:-2]):
            out2.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
            written_count += 1
      sys.stderr.write('%s: Total %i, Written %i (%.1f%%)\n' % (self.f[1], total_count, written_count, written_count/total_count*100))
      
      out1.close()
      out2.close()
      fh1.close()
      fh2.close()
      rm_files([db1_fname, db2_fname, dbcommon_fname, f1, f2])
   
   def trim(self, paired=False):
      '''Start trimming of single end reads'''
      
      if paired:
         rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(36))
      
      for i,f in enumerate(self.f):
         written = 0
         total = 0      
         fh_in = open(f, 'r')
         
         if paired: 
            fh_out = open('tmp'+str(i)+'.'+rand, 'w')
         else:
            fh_out = open(self.o[0], 'w')
         
         for (title, sequence, quality) in FastqGeneralIterator(fh_in):
            total += 1
            (title, sequence, quality) = self.filter_adaptor(title, sequence, quality)
            (title, sequence, quality) = self.trim_qual(title, sequence, quality)
            if title != None:
               if len(sequence) != len(quality):
                  raise ValueError('sequence and quality not of the same length\n%s\n%s\n' % (sequence, quality))
               fh_out.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
               written += 1
         fh_out.close()
         fh_in.close()
               
         sys.stderr.write('%s: written %s of %s (%.1f%%)\n' % (self.f[i], written, total, written/total))
      
      if paired:
         self.write_pairs('tmp0.'+rand, 'tmp1.'+rand)

if __name__ == '__main__':
      
   parser = argparse.ArgumentParser(description='''
      Trim SE or PE files for low q bases, adaptor sequences, Ns
      ''')
   
   # add the arguments
   parser.add_argument('--i', help='input single or paired end files', nargs='+', required=True)
   parser.add_argument('--min_length', help='minimum length of a read to keep pairs [25]', type=int, default=25)
   parser.add_argument('--min_baseq', help='chomp bases with quality less than [20]', default=20, type=int)
   parser.add_argument('--min_avgq', help='minimum average quality of read [20]', default=20, type=int)
   parser.add_argument('--adaptors', help='adaptor sequence to clip [Illumina adaptors]', nargs='+', default=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'])
   parser.add_argument('--keep_n', help='do not remove sequences containing N', default=False, action='store_true')
   parser.add_argument('--min_adaptor_match', help='minimum length of match to adaptor (0=all of adaptor) [20]', default=20, type=int)
   parser.add_argument('--o', help='output files', nargs='+', required=True)
   parser.add_argument('--log', help='log level [INFO]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args(''.split())
   #args = parser.parse_args('--i kleb_test_2.fq --l 25 --q 20 --o kleb_test_2.trim.fq'.split())
   #args = parser.parse_args('--i Kleb-10-213361_2_1_sequence.txt Kleb-10-213361_2_2_sequence.txt --M 15 --o Kleb-10-213361_2_1_sequence.trim.fq Kleb-10-213361_2_2_sequence.trim.fq '.split())
   
   # create instance
   fqtrim = FastqTrim(args.i, args.o, args.min_length, args.min_baseq, args.min_avgq, args.keep_n, args.min_adaptor_match, args.adaptors)
   
   # start trimming
   if len(fqtrim.f) == 2: 
      fqtrim.trim(paired=True)
   elif len(fqtrim.f) == 1: 
      fqtrim.trim()
   else:
      raise ValueError('Only 1 (single end) or 2 (paired end) files must be given')
