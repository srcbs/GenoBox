#!/panvol1/simon/bin/python2.7

from __future__ import division

import argparse
import os
import logging
import re
import sys

# create the parser
parser = argparse.ArgumentParser(description='''
   Filter vcf based on read depths, quality score and allele frequency
   Reads from stdin or file and writes to stdout or file
   Filters X,Y and MT for heterozygote. Add flag for not to do that on X if female
   ''')

# add the arguments
parser.add_argument('--i', help='input vcf [stdin]', default=None)
parser.add_argument('--d', help='minimum read depth', type=int, default=10)
parser.add_argument('--D', help='maximum read depth', type=int, default=100)
parser.add_argument('--Q', help='minimum quality score', type=float, default=20.0)
parser.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser.add_argument('--refonly', help='only print reference (no SNPs) [False]', default=False, action='store_true')
parser.add_argument('--t', help='type of vcf (samtools, samtools_svn, gatk) [samtools]', default='samtools')
parser.add_argument('--o', help='output vcf [stdout]', default=None)
parser.add_argument('--log', help='log level [INFO]', default='info')

# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--i Aborigine_hg19.all.chr21.vcf --o Aborigine_hg19.all.chr21.flt.vcf --D 50 --ex /panvol1/simon/databases/hs_ref37/gi2number.build37 --refonly'.split())

logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

if args.ex:
   d = {}
   bhandle = open(args.ex, 'r')
   for line in bhandle:
      line = line.rstrip()
      fields = line.split('\t')
      d[fields[0]] = fields[1]
   bhandle.close()

if args.o:
   outhandle = open(args.o, 'w')

written = 0
linecount = 0
hetzygote = 0
homzygote = 0
hetwritten = 0
homwritten = 0
same = 0
variants = 0

# set regular expressions
reg = re.compile('DP=(\d+)')
reg_dp4 = re.compile('DP4=(\d+),(\d+),(\d+),(\d+);')
if args.t == 'samtools':
   reg_state = re.compile('\:(\d\/\d)\:')
elif args.t == 'samtools_svn':
   reg_state = re.compile('(\d\/\d):')

reg_alt = re.compile('(\w+),')
if args.i:
   line = line.rstrip()
   for line in open(args.i, 'r'):
      if line.startswith('#'):
         outhandle.write(line) if args.o else sys.stdout.write(line)
         continue
      
      linecount += 1
      fields = line.split('\t')
      chrom, pos, id, ref, alt, qual, filter, info, pl, hs = fields
      depth = int(reg.search(info).group(1))
      
      # change chromosome name
      if args.ex:
         chrom = d[chrom]
         fields[0] = chrom
         line = '\t'.join(fields)
      
      # apply filters
      if depth < args.d or depth > args.D:
         continue
      if float(qual) < args.Q:
         continue
      
      # handle no variants
      if alt == '.':
         if info.startswith('INDEL'):
            continue
         same += 1
         outhandle.write(line) if args.o else sys.stdout.write(line)
         continue
      
      if args.refonly:
         continue
      
      variants += 1
      # handle variants
      counts = reg_dp4.findall(info)
      if counts:
         ref=int(counts[0][0]) + int(counts[0][1])
         nonref= int(counts[0][2]) + int(counts[0][3])
         #if ref == 0 and nonref == 0:
         #   freq = 0
         if ref < nonref:
            freq = ref / (ref+nonref)
         else:
            freq = nonref/ (ref+nonref)
      else:
         raise ValueError('error in DP4 counting')
      
      # check if there are alternative calls that should be removed eg (A,T)
      alt_match = reg_alt.search(alt)
      if alt_match:
         alt = alt_match.group(1)
         fields[4] = alt
         line = '\t'.join(fields)
      
      # only print out if allele frequency is >0.2 or <0.2 for het and hom, respectively
      state_match = reg_state.search(hs)
      if state_match.group(1) == '0/1':
         if chrom == 'X' or chrom == 'Y' or chrom == 'MT':
            continue
         hetzygote += 1
         if freq >= 0.2:
            outhandle.write(line) if args.o else sys.stdout.write(line)
            hetwritten += 1
            written += 1
      elif state_match.group(1) == '1/1':
         homzygote += 1
         if freq < 0.2:
            outhandle.write(line) if args.o else sys.stdout.write(line)
            homwritten += 1
            written += 1
      else:
         raise ValueError('error in parsing hs state')

else:
   for line in sys.stdin:
      #line = line.rstrip()
      if line.startswith('#'):
         outhandle.write(line) if args.o else sys.stdout.write(line)
         continue
      
      linecount += 1
      fields = line.split('\t')
      chrom, pos, id, ref, alt, qual, filter, info, pl, hs = fields
      depth = int(reg.search(info).group(1))
      
      # change chromosome name
      if args.ex:
         chrom = d[chrom]
         fields[0] = chrom
         line = '\t'.join(fields)
      
      # apply filters
      if depth < args.d or depth > args.D:
         continue
      if float(qual) < args.Q:
         continue
      
      # handle no variants
      if alt == '.':
         if info.startswith('INDEL'):
            continue
         same += 1
         outhandle.write(line) if args.o else sys.stdout.write(line)
         continue
      
      if args.refonly:
         continue
      
      variants += 1
      # handle variants
      counts = reg_dp4.findall(info)
      if counts:
         ref=int(counts[0][0]) + int(counts[0][1])
         nonref= int(counts[0][2]) + int(counts[0][3])
         #if ref == 0 and nonref == 0:
         #   freq = 0
         if ref < nonref:
            freq = ref / (ref+nonref)
         else:
            freq = nonref/ (ref+nonref)
      else:
         raise ValueError('error in DP4 counting')
      
      # check if there are alternative calls that should be removed eg (A,T)
      alt_match = reg_alt.search(alt)
      if alt_match:
         alt = alt_match.group(1)
         fields[4] = alt
         line = '\t'.join(fields)
      
      # only print out if allele frequency is >0.2 or <0.2 for het and hom, respectively
      state_match = reg_state.search(hs)
      if state_match.group(1) == '0/1':
         if chrom == 'X' or chrom == 'Y' or chrom == 'MT':
            continue
         hetzygote += 1
         if freq >= 0.2:
            outhandle.write(line) if args.o else sys.stdout.write(line)
            hetwritten += 1
            written += 1
      elif state_match.group(1) == '1/1':
         homzygote += 1
         if freq < 0.2:
            outhandle.write(line) if args.o else sys.stdout.write(line)
            homwritten += 1
            written += 1
      else:
         raise ValueError('error in parsing hs state')

if args.o:
   outhandle.close()

try:
   logger.info('%s: Total=%i, Written=%i (%.1f%%), Heterozyg_written=%i (%.1f%%), Homozyg_written=%i (%.1f%%) ' % (args.i, linecount, written+same, (written+same)/linecount*100, hetwritten, hetwritten/hetzygote*100, homwritten, homwritten/homzygote*100))
except:
   pass

