#!/usr/bin/python

import argparse
import os
import logging
import re
import sys

def read_indels_vcf(indel_vcf):
   '''Read vcf of indels and return dict where keys are tuples (chr, pos)
   that are deletions'''
   
   fh = open(indel_vcf, 'r')
   deletions_dict = {}
   for line in fh:
      line = line.rstrip()
      if line.startswith('#'):
         continue
      
      fields = line.split('\t')
      
      # get deletions (alt is shorter than ref)
      diff = len(fields[3]) - len(fields[4])
      if diff > 0:
         currpos = fields[1]
         for i in range(diff):
            deletions_dict[(fields[0], str(int(currpos)+(i+1)))] = 1
      
   return deletions_dict

def parse_vcf(vcf, deletions_dict, rmpos_file, out):
   '''Remove deletion positions from vcf and print these to file'''
   
   counts = {}
   counts['total'] = 0
   counts['passed'] = 0
   counts['filtered'] = 0
   
   # open handles
   if vcf:
      fh_vcf = open(vcf, 'r')
   else:
      fh_vcf = sys.stdin
   
   if out:
      fh_out = open(out, 'w')
   else:
      fh_out = sys.stdout
   
   fh_rmpos = open(rmpos_file, 'w')
   
   # start parse
   for line in fh_vcf:
      if line.startswith('#'):
         fh_out.write(line)
         continue
      
      counts['total'] += 1
      fields = line.split('\t')
      if deletions_dict.has_key((fields[0], fields[1])):
         counts['filtered'] += 1
         fh_rmpos.write(line)
      else:
         counts['passed'] += 1
         fh_out.write(line)
   
   fh_vcf.close()
   fh_out.close()
   fh_rmpos.close()
   return counts

# create the parser
parser = argparse.ArgumentParser(description='''
   Filter indels (deletions) positions from same-as-reference files
   ''')

# add the arguments
parser.add_argument('--i', help='input vcf [stdin]', default=None)
parser.add_argument('--indels', help='indels vcf to read')
parser.add_argument('--rmpos', help='file to print removed positions to [filtered.pos]', default='filtered.pos')
parser.add_argument('--o', help='output vcf [stdout]', default=None)
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--i Aborigine_trimmed_hg19.flat.ref.chrY.flt.ann.nr.test.vcf --indels indels_for_filtering.vcf --o AAborigine_trimmed_hg19.flat.ref.chrY.flt.ann.nr.test.indelfilt.vcf --rmpos filtered_pos.chrY.vcf'.split())

logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# read indels
deletions_dict = read_indels_vcf(args.indels)

# filter positions and write if they are filtered
counts = parse_vcf(args.i, deletions_dict, args.rmpos, args.o)

logger.info('%s: total %i, passed %i, filtered %i' % (args.i, counts['total'], counts['passed'], counts['filtered']))
