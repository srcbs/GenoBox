#!/panvol1/simon/bin/python2.7

import argparse
import sys
from collections import defaultdict
import gzip

def create_error_matrix(i,o):
   '''Calculate transition / transversion ratio of vcf'''
   
   # open handles
   if i == 'stdin':
      fh = sys.stdin
   else:
      fh = gzip.open(i, 'r')
   
   if args.o == 'stdout':
      fh_out = sys.stdout
   else:
      fh_out = open(o, 'w')
   
   # create error matrix
   # observed A C T G and expected AA CC TT GG
   AA = {'A':0, 'C':0, 'G':0, 'T':0}
   CC = {'A':0, 'C':0, 'G':0, 'T':0}
   GG = {'A':0, 'C':0, 'G':0, 'T':0}
   TT = {'A':0, 'C':0, 'G':0, 'T':0}
   
   # parse
   total = 0
   for line in fh:
      total += 1
      fields = line.split()
      mp = fields[4]
      gt = fields[2]
      count = defaultdict(int)
      
      # skip if indels
      if mp.find('-') > -1 or mp.find('+') > -1:
         continue            
      
      # count
      count['same'] = mp.count('.') + mp.count(',')
      count['A'] = mp.count('A') + mp.count('a')
      count['C'] = mp.count('C') + mp.count('c')
      count['G'] = mp.count('G') + mp.count('g')
      count['T'] = mp.count('T') + mp.count('t')
      
      # add counts to error matrix
      if gt == 'A':
         AA['A'] = AA['A'] + count['same']
         AA['C'] = AA['C'] + count['C']
         AA['G'] = AA['G'] + count['G']
         AA['T'] = AA['T'] + count['T']
      elif gt == 'C':
         CC['C'] = CC['C'] + count['same']
         CC['A'] = CC['A'] + count['A']
         CC['G'] = CC['G'] + count['G']
         CC['T'] = CC['T'] + count['T']
      elif gt == 'G':
         GG['G'] = GG['G'] + count['same']
         GG['A'] = GG['A'] + count['A']
         GG['C'] = GG['C'] + count['C']
         GG['T'] = GG['T'] + count['T']
      elif gt == 'T':
         TT['T'] = TT['T'] + count['same']
         TT['A'] = TT['A'] + count['A']
         TT['C'] = TT['C'] + count['C']
         TT['G'] = TT['G'] + count['G']
      elif gt == 'N':
         continue
      else:
         sys.stderr.write('Non ACGTN genotype found at line %i\n' % total)
         sys.stderr.write('%s' % line)
   
   # in the end write the error matrix
   header = '\tA\tC\tG\tT\n'
   AA_line = 'AA\t%i\t%i\t%i\t%i\n' % (AA['A'], AA['C'], AA['G'], AA['T'])
   CC_line = 'CC\t%i\t%i\t%i\t%i\n' % (CC['A'], CC['C'], CC['G'], CC['T'])
   GG_line = 'GG\t%i\t%i\t%i\t%i\n' % (GG['A'], GG['C'], GG['G'], GG['T'])
   TT_line = 'TT\t%i\t%i\t%i\t%i\n' % (TT['A'], TT['C'], TT['G'], TT['T'])
   
   fh_out.write(header)
   fh_out.write(AA_line)
   fh_out.write(CC_line)
   fh_out.write(GG_line)
   fh_out.write(TT_line)
   fh_out.close()

# create the parser
parser = argparse.ArgumentParser(description='''
   Create error matrix from mpileup
   ''')

# add the arguments
parser.add_argument('--i', help='input mpileup (.gz or stdin) [stdin]', default='stdin')
parser.add_argument('--o', help='output error matrix [stdout]', default='stdout')


# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--i tmp.mpileup --o tmp.em'.split())

if __name__ == '__main__':
   create_error_matrix(args.i, args.o)
