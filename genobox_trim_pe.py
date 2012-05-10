#!/panvol1/simon/bin/python2.7

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import genobox_modules
import sys
import gzip
import os
import multiprocessing

def trim_reads(args):
   '''Trim fastq file(s) of reads'''
      
   def set_adaptors(alist, M):
      '''Set adaptors according to adaptor length requirements'''
      
      adaptors = []
      for a in alist:
         if M == 0:
            adaptor = a
         elif (M) > len(a):
            adaptor = a
         else:
            adaptor = a[:M]
         adaptors.append(adaptor)
      return adaptors
   
   def filter_adaptor(adaptors, min_adaptor_match, min_length, title, sequence, quality):
      '''Filters adaptor and Ns in read'''
      
      has_adaptor = False
      has_N = False
      is_short = False
      
      # search for adaptor in reads
      for adaptor in adaptors:
         hit = sequence.find(adaptor)
         if hit > -1:
            sequence = sequence[:hit]
            quality = quality[:hit]
            has_adaptor = True
      
      # check if read contains Ns and length
      if sequence.find('N') > -1: has_N = True
      if len(sequence) < min_length: is_short = True
      
      # return read
      if has_N or is_short:
         return (None, None, None)
      else:
         return (title, sequence, quality)
   
   # the actual read trimming function
   def trim_qual(min_baseq, min_avgq, min_length, fqtype, title, sequence, quality):
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
      
      # calc basesq
      if fqtype == 'Illumina':
         qual = illumina2qual(quality)
      elif fqtype == 'Sanger':
         qual = sanger2qual(quality)
      else:
         raise ValueError('Fastq must be Illumina or Sanger')
      
      # check trailing qualities, should it be trimmed?
      t = trim_index(qual, min_baseq)
      if t:
         sequence = sequence[:t]
         quality = quality[:t]
         a = average(qual[:t])
         if a < min_avgq or len(sequence) < min_length:
            return (None, None, None)            
         else:
            return (title, sequence, quality)
      else:
         # do not print out if average quality is less than m or length is less than l
         a = average(qual)
         if a < min_avgq or len(sequence) < min_length:
            return (None, None, None)
         else:
            return (title, sequence, quality)
         
   def start_trim(f, lock, conn, args):
      '''Start the actual trimming of a file'''
      
      if args.gz: fh = gzip.open(f, 'rb')
      else: fh = open(f, 'r')
      
      # set adaptors #
      adaptors = set_adaptors(args.adaptors, args.min_adaptor_match)
      
      # set fqtype #
      format = genobox_modules.set_filetype(f, args.gz)
      if format == 'fastq': fqtype = genobox_modules.set_fqtype(f, args.gz)
      else: raise ValueError('Input not fastq\n')
      
      # start trimming file #
      total = 0
      written = 0
      for (title, sequence, quality) in FastqGeneralIterator(fh):
         #lock.acquire()
         total += 1
         (title, sequence, quality) = filter_adaptor(adaptors, args.min_adaptor_match, args.min_length, title, sequence, quality)
         (title, sequence, quality) = trim_qual(args.min_baseq, args.min_avgq, args.min_length, fqtype, title, sequence, quality)
         if title != None:
            if len(sequence) != len(quality):
               raise ValueError('sequence and quality not of the same length\n%s\n%s\n' % (sequence, quality))
            written += 1
         conn.send('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
         #lock.release()
      conn.send('Stop')
      conn.close()
      return written, total
   
   # main #
   
   # open output files
   f0 = os.path.split(args.i[0])
   f1 = os.path.split(args.i[1])
   if args.gz:
      trima = gzip.open("trimmed/%s.trim.fq.gz" % f0[1], 'wb')
      trimb = gzip.open("trimmed/%s.trim.fq.gz" % f1[1], 'wb')
      trima_s = gzip.open("trimmed/%s.trim.fq.single.gz" % f0[1], 'wb')
      trimb_s = gzip.open("trimmed/%s.trim.fq.single.gz" % f1[1], 'wb')
   else:
      trima = open("trimmed/%s.trim.fq" % f0[1], 'w')
      trimb = open("trimmed/%s.trim.fq" % f1[1], 'w')
      trima_s = open("trimmed/%s.trim.single.fq" % f0[1], 'w')
      trimb_s = open("trimmed/%s.trim.single.fq" % f1[1], 'w')
   
   print >>sys.stderr, "writing %s and %s" % (trima.name, trimb.name)      
   
   # start processing #
   lock = multiprocessing.Lock()
   parent_conn1, child_conn1 = multiprocessing.Pipe()
   p1 = multiprocessing.Process(target=start_trim, args=(args.i[0], lock, child_conn1, args,), )
   p1.start()
   
   parent_conn2, child_conn2 = multiprocessing.Pipe()
   p2 = multiprocessing.Process(target=start_trim, args=(args.i[1], lock, child_conn2, args,), )
   p2.start()
      
   # parsing reads, checking
   paired, single1, single2, total, removed1, removed2 = (0, 0, 0, 0, 0, 0)
   while True:
      read1, read2 = (parent_conn1.recv(), parent_conn2.recv())
      if read1 == 'Stop':
         if read2 == 'Stop': 
            sys.stderr.write('%s, total: %i, paired: %i, single: %i, removed: %i\n' % (args.i[0], total, paired, single1, removed1))
            sys.stderr.write('%s, total: %i, paired: %i, single: %i, removed: %i\n' % (args.i[1], total, paired, single2, removed2))
            sys.exit()
         else: raise ValueError('files not equal length')
      
      total += 1
      
      # check new Illumina header
      if read1.split('\n')[0].find(' ') > -1:
         h1, h2 = (read1.split('\n')[0].split(' ')[0], read2.split('\n')[0].split(' ')[0])
      else:
         h1, h2 = (read1.split('\n')[0][:-2], read2.split('\n')[0][:-2])
      
      if read1 == '@None\nNone\n+\nNone\n' and read2 == '@None\nNone\n+\nNone\n': 
         removed1 += 1
         removed2 += 1
         continue
      else:
         if h1 == h2:
            trima.write(read1)
            trimb.write(read2)
            paired += 1
         elif read1 != '@None\nNone\n+\nNone\n': 
            trima_s.write(read1)
            single1 += 1
            removed2 += 1
         elif read2 != '@None\nNone\n+\nNone\n': 
            trimb_s.write(read2)
            single2 += 1
            removed1 += 1
         else: raise ValueError('Error in pair file parse\n')
         

if __name__ == '__main__':
      
   parser = argparse.ArgumentParser(prog='genobox_trim_pe.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50,width=140),
                                 description='''Trim PE files for 3' quals, average quals, adaptors, Ns, read length. Will output common and unique reads in seperate files.''', 
                                 usage='%(prog)s [options]')
   
   # add the arguments
   parser.add_argument('--i', help='input paired end files', nargs=2, required=True)
   parser.add_argument('--min_length', help='minimum length of a read to keep pairs [25]', type=int, default=25)
   parser.add_argument('--min_baseq', help='chomp bases with quality less than [20]', default=20, type=int)
   parser.add_argument('--min_avgq', help='minimum average quality of read [20]', default=20, type=int)
   parser.add_argument('--adaptors', help='adaptor sequence to clip [Illumina adaptors]', nargs='+', default=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'])
   parser.add_argument('--keep_n', help='do not remove sequences containing N', default=False, action='store_true')
   parser.add_argument('--min_adaptor_match', help='minimum length of match to adaptor (0=all of adaptor) [20]', default=20, type=int)
   parser.add_argument('--gz', help='input files are gzipped [False]', default=False, action='store_true')
   parser.add_argument('--log', help='log level [INFO]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--i Kleb-10-213361_2.interleaved.fastq --o Kleb-10-213361_2_1.interleaved.fastq.trim.fastq Kleb-10-213361_2_2.interleaved.fastq.trim.fastq'.split())
   #args = parser.parse_args('--i test_kleb_1.fq test_kleb_2.fq --min_length 25 --min_baseq 20 --min_avgq 20 --min_adaptor_match 20'.split())
   #args = parser.parse_args('--i tmp_1.fastq tmp_2.fastq --min_length 25 --min_baseq 20 --min_avgq 20 --min_adaptor_match 20'.split())

   trim_reads(args)
   
   
