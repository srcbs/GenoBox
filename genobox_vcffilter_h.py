#!/usr/bin/python

from __future__ import division

import argparse
import subprocess
import pipelinemod
import os
import logging

def get_genome(chr_file):
   '''Read genome file into list of list'''
   
   if chr_file:
      fh = open(chr_file, 'r')
      L = []
      for line in fh:
         line = line.rstrip()
         s = line.split('\t')
         L.append(s)
      return L
   else:
      raise IOError('genome file not given, input as --genome <file>')

def bcf2varfilter(bcf, genome, Q, vcf_prefix):
   '''Runs bcf through varfilter and writes to vcf'''
   
   paths = pipelinemod.setSystem()
   bcf_cmd = paths['samtools_svn_home'] + 'bcftools view'
   calls = []
   vcf_files = []
   vcfutils_cmd = paths['samtools_svn_home'] + 'vcfutils.pl'
   for chr in genome:
      d = chr[4]
      D = chr[5]
      vcf = vcf_prefix + chr[2] + '.vcf'
      vcf_files.append(vcf)
      
      bcf_arg = ' %s \"%s\"' % (bcf, chr[0])
      bcf_call = bcf_cmd + bcf_arg
      
      vcfutils_arg = ' varFilter -d%s -D%s' % (d, D)
      vcfutils_call = vcfutils_cmd + vcfutils_arg
      
      qualf_call = """ perl -ane 'if ($_ =~ m/^#/) { print $_ } else { if ($F[5] > %f) { print $_ }}' > %s""" % (Q, vcf)
      call = '%s | %s | %s' % (bcf_call, vcfutils_call, qualf_call)
      logger.info(call)
      subprocess.check_call(call, shell=True)
   return vcf_files

def cat_vcfs(vcfs, vcf_out):
   '''Concatenate vcf to a single vcf'''
   
   cmd = 'cat'
   arg = ' %s > %s' % (' '.join(vcfs), vcf_out)
   call = cmd + arg
   logger.info(call)
   subprocess.check_call(call, shell=True)
   

def vcf_filter_rmsk(vcf, rmsk, vcf_out):
   '''Removes variants called inside annotated repeat
   If no rmsk is given it simply copies the file'''
   
   if rmsk and rmsk != 'None':
      # add header
      head_call = 'head -n 1000 %s | grep "#" > %s' % (vcf, vcf_out)
      logger.info(head_call)
      subprocess.check_call(head_call, shell=True)
      
      # perform rmsk filtering
      bed_call = paths['bedtools_home'] + 'intersectBed'
      bed_arg = ' -v -a %s -b %s >> %s' % (vcf, rmsk, vcf_out)
      call = bed_call + bed_arg
      logger.info(call)
      subprocess.check_call(call, shell=True)
   else:
      call = 'cp %s %s' % (vcf, vcf_out)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def vcf_filter_haploid(vcf, genome, vcf_out):
   '''Filter heterozygote variants haploid chromosomes'''
   
   import re
   # parse ploidy to dict
   genome_dict = {}
   for c in genome:
      genome_dict[c[0]] = c[3]
   
   count = 0
   state_reg = re.compile('\t(\d)\/(\d):')
   fh = open(vcf, 'r')
   fh_out = open(vcf_out, 'w')
   for line in fh:
      if line.startswith('#'):
         fh_out.write(line)
         continue
      
      count += 1
      fields = line.split()
      if genome_dict[fields[0]] == 'diploid':
         fh_out.write(line)
         pass
      elif genome_dict[fields[0]] == 'haploid':
         state = state_reg.search(line)
         if state.groups()[0] == state.groups()[1]:
            fh_out.write(line)
         else:
            pass
      elif genome_dict[fields[0]] == 'na':
         pass
      else:
         raise ValueError('ploidy of %s is not diploid/haploid/na' % fields[0])
   fh.close()
   fh_out.close()

def vcf_filter_allelic_balance(vcf, threshold, caller, vcf_out):
   '''Filter allelic balance, heterozygote must be above, homozygote must be below'''
   
   import re
   
   # set re's
   state_re = re.compile('\t(\d)\/(\d):')
   if caller == 'samtools':
      depth_re = re.compile('DP4=(\d+),(\d+),(\d+),(\d+);')
   elif caller == 'gatk':
      depth_re = re.compile('\d\/\d:(\d+),(\d+)')
   
   if threshold != 0.5:
      # parse
      count = 0
      hetzygote = 0
      homzygote = 0
      fh = open(vcf, 'r')
      fh_out = open(vcf_out, 'w')
      for line in fh:
         if line.startswith('#'):
            fh_out.write(line)
            continue
         
         count += 1
         fields = line.split()
         
         # get allelic counts
         counts = depth_re.search(line).groups()
         if caller == 'samtools':
            ref=int(counts[0]) + int(counts[1])
            nonref= int(counts[2]) + int(counts[3])
         elif caller == 'gatk':
            ref=int(counts[0])
            nonref=int(counts[1])
         
         # calc frequency
         if ref < nonref:
            freq = ref / (ref+nonref)
         else:
            freq = nonref/ (ref+nonref)
         
         # filter based on state
         state = state_re.search(line)
         if state.groups()[0] == state.groups()[1]:
            homzygote += 1
            if freq < threshold:
               fh_out.write(line)
         elif state.groups()[0] != state.groups()[1]:
            hetzygote += 1
            if freq >= threshold:
               fh_out.write(line)
   else:
      call = 'cp %s %s' % (vcf, vcf_out)
      logger.info(call)
      subprocess.check_call(call)
      
   
def vcf_filter_prune(vcf, prune, vcf_out):
   '''Prune variants within N nt of each other'''
   
   paths = pipelinemod.setSystem()
   
   if prune != 0:
      # create header
      head_call = 'head -n 1000 %s | grep "#" > %s' % (vcf, 'genotyping/header.vcf')
      logger.info(head_call)
      subprocess.check_call(head_call, shell=True)
      
      tmp_file = vcf + '.tmp'
      prune_script = paths['pyscripts_home'] + 'saqqaq_snppruning.R'
      prune_cmd = paths['R_home'] + 'R-2.12'
      prune_arg = ' --vanilla %i %s %s < %s' % (prune, vcf, tmp_file, prune_script)
      prune_call = prune_cmd + prune_arg
      logger.info(prune_call)
      subprocess.check_call(prune_call, shell=True)
      
      # add header
      header_call = 'cat genotyping/header.vcf %s > %s' % (tmp_file, vcf_out)
      logger.info(header_call)
      subprocess.check_call(header_call, shell=True)
      
      # rm tmp_files
      rm_call = 'rm %s genotyping/header.vcf' % tmp_file
      logger.info(rm_call)
      subprocess.check_call(rm_call, shell=True)
   else:
      call = 'cp %s %s' % (vcf, vcf_out)
      logger.info(call)
      subprocess.check_call(call)

# create the parser
parser = argparse.ArgumentParser(description=
   '''Filter variants (snps, indels).
   
   Genome file must be given, format is a line for each chromosome:
   chrom\tchrom_len\tchrom_short_name\haploid/diploid\tlow_depth\thigh_depth
   
   Filtering steps:
   vcfutils.pl varFilter
   annotated repeats using rmsk
   heterozygote variants on haploid chromosomes
   allelic balance
   pruning of variants within N nt of each other
   ''')

# add the arguments
parser.add_argument('--bcf', help='input bcf var file')
parser.add_argument('--genome', help='file containing chromosomes to analyse, format: chrom\tchrom_len\tchrom_short_name\ploidy\tlow_d\thigh_d', default=None)
parser.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser.add_argument('--Q', help='minimum quality score', type=float, default=20.0)
parser.add_argument('--rmsk', help='rmsk to use', default=None)
parser.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50)
parser.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0)
parser.add_argument('--o', help='output vcffile [snpcalls.vcf]', default='snpcalls.vcf')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--bcf Aborigine_trimmed_hg19.flat.var.bcf --chr build37_rCRS.genome --Q 40 --rmsk /panvol1/simon/databases/hs_ref37_rCRS/rmsk/rmsk_build37_rCRS.sort.gi.genome --prune 5 --ab 0.2'.split())
#args = parser.parse_args('--bcf /panvol1/simon/projects/cge/test/kleb_10_213361/genotyping/kleb_10_213361.bcf --genome /panvol1/simon/projects/cge/test/kleb_10_213361/kleb_pneu.genome --caller samtools --Q 40.000000 --rmsk None --ab 20.000000 --prune 5 --o genotyping/kleb_10_213361.vcf'.split())

# set logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set queuing
paths = pipelinemod.setSystem()
home = os.getcwd()

# get genome file
genome = get_genome(args.genome)

# perform varfilter to get filtered vcf
vcf_files = bcf2varfilter(args.bcf, genome, args.Q, 'genotyping/tmp.flt.')

# combine to one file
cat_vcfs(vcf_files, 'genotyping/tmp.flt.all.vcf')

# remove in annotated repeat (rmsk)
vcf_filter_rmsk('genotyping/tmp.flt.all.vcf', args.rmsk, 'genotyping/tmp.flt.all.rmsk.vcf')

# filter haploid chromosomes for heterozygote calls
vcf_filter_haploid('genotyping/tmp.flt.all.rmsk.vcf', genome, 'genotyping/tmp.flt.all.rmsk.hetfilt.vcf')

# filter for allelic balance
vcf_filter_allelic_balance('genotyping/tmp.flt.all.rmsk.hetfilt.vcf', args.ab, args.caller, 'genotyping/tmp.flt.all.rmsk.hetfilt.abfilt.vcf')

# pruning of nearby calls
vcf_filter_prune('genotyping/tmp.flt.all.rmsk.hetfilt.abfilt.vcf', args.prune, args.o)

# remove temporary files
pipelinemod.rm_files(['genotyping/tmp.flt*'])

