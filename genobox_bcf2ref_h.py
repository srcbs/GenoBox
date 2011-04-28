#!/usr/bin/python

import argparse
import os
import logging
import re
import subprocess
import pipelinemod

def vcf_filterAll(bcf, chr_id, d, D, Q, ex, vcf_out_gz):
   '''Extracts bcf through genobox_vcffilterAll.py to vcf.gz'''
   
   paths = pipelinemod.setSystem()
   
   bcf_cmd = paths['samtools_home'] + 'bcftools view'
   if chr_id:
      bcf_arg = ''' %s \"%s\" ''' % (bcf, chr_id)
   else:
      bcf_arg = ''' %s ''' % (bcf)
   bcf_call = bcf_cmd + bcf_arg
   
   vcf_filter_cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_vcffilterAll.py'
   if ex and ex != 'None':
      vcf_filter_arg = ' --d %i --D %i --Q %f --ex %s --refonly' % (d, D, Q, ex)
   else:
      vcf_filter_arg = ' --d %i --D %i --Q %f --refonly' % (d, D, Q)
   
   vcf_filter_call = vcf_filter_cmd + vcf_filter_arg
   
   call = '%s | %s | bgzip -c > %s' % (bcf_call, vcf_filter_call, vcf_out_gz)
   logger.info(call)
   subprocess.check_call(call, shell=True)

def vcf_tabix(vcf_gz):
   '''Run tabix on vcf.gz'''
   
   paths = pipelinemod.setSystem()
   
   tabix_call = paths['bin_home'] + 'tabix -p vcf -f %s' % (vcf_gz)
   logger.info(tabix_call)
   subprocess.check_call(tabix_call, shell=True)

def vcf_annotate_dbsnp(vcfgz, dbsnp, vcf_out_gz):
   '''Annotate vcf.gz with dbsnp'''
   
   paths = pipelinemod.setSystem()
   
   if dbsnp and dbsnp != 'None':
      gunzip_call = '/usr/bin/gunzip -c %s' % vcfgz
      fill_call = paths['bin_home'] + 'fill-rsIDs -r %s | bgzip -c > %s' % (dbsnp, vcf_out_gz)
      
      dbsnp_call = '%s | %s' % (gunzip_call, fill_call)
      logger.info(dbsnp_call)
      subprocess.check_call(dbsnp_call, shell=True)
   else:
      call = 'cp %s %s' % (vcfgz, vcf_out_gz)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def vcf_filter_rmsk(vcfgz, rmsk, vcfgz_out):
   '''Removes variants called inside annotated repeat
   If no rmsk is given it simply copies the file'''
   
   import random
   import string
   import pipelinemod
   
   paths = pipelinemod.setSystem()
   if rmsk and rmsk != 'None':
      # create header
      N = 10
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))
      header = 'genotyping/tmp'+rand+'.header.vcf'
      header_call = '/usr/bin/gunzip -c %s | head -n 1000 | grep "#" > %s' % (vcfgz, header)
      logger.info(header_call)
      subprocess.check_call(header_call, shell=True)
      
      # perform rmsk filtering
      gunzip_call = '/usr/bin/gunzip -c %s' % vcfgz
      bgzip_call = paths['bin_home'] + 'bgzip -c > %s' % vcfgz_out
      
      bed_cmd = paths['bedtools_home'] + 'intersectBed'
      bed_arg = ' -v -a stdin -b %s | cat %s - | %s' % (rmsk, header, bgzip_call)
      bed_call = bed_cmd + bed_arg
      call = '%s | %s' % (gunzip_call, bed_call)
      logger.info(call)
      subprocess.check_call(call, shell=True)
      
      # rm tmp header file
      subprocess.check_call('rm %s' % header, shell=True)
   else:
      call = 'cp %s %s' % (vcfgz, vcfgz_out)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def read_rmsk(rmsk, chr):
   '''Read rmsk file and return dictionary of positions to filter
   only keep positions that matches chromosome being analyzed'''
   
   # because this is only used for MT dna only keep chromosomal positions from rmsk file
   # that matches MT (to minimize memory)
   fh = open(rmsk, 'r')
   dict = {}
   if chr.find('MT') > -1:
      for line in fh:
         line = line.rstrip()
         fields = line.split('\t')
         if fields[0] == 'MT':
            for pos in xrange(int(fields[1]), int(fields[2])+1, 1):
               dict[fields[0], str(pos)] = 1
   else:
      for line in fh:
         line = line.rstrip()
         fields = line.split('\t')
         for pos in xrange(int(fields[1]), int(fields[2])+1, 1):
            dict[fields[0], str(pos)] = 1
   
   return dict

def manual_rmsk_filter(vcfgz, chr, rmsk, vcfgz_out):
   '''Manually filter vcfgz using rmsk file'''
   
   if rmsk and rmsk != 'None':
      import gzip
      
      # get rmsk as dict
      rmsk_dict = read_rmsk(rmsk, chr)
      
      # loop over positions
      fh = gzip.open(vcfgz, 'rb')
      fh_out = gzip.open(vcfgz_out, 'wb')
      for line in fh:
         if line.startswith('#'):
            fh_out.write(line)
            continue
         
         fields = line.split('\t')
         if rmsk_dict.has_key((fields[0], fields[1])):
            pass
         else:
            fh_out.write(line)
   else:
      call = 'cp %s %s' % (vcfgz, vcfgz_out)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def vcf_filter_indels(vcfgz, chr, indels, rmpos, o):
   '''Removes same-as-ref positions covered by indels'''
   
   if indels and indels != 'None':
      paths = pipelinemod.setSystem()
      gzcat_call = '/usr/bin/gunzip -c %s' % vcfgz
      
      vcf_filter_cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_vcffilterindels.py'
      vcf_filter_arg = ' --indels %s --rmpos %s' % (indels, rmpos)
      vcf_filter_call = vcf_filter_cmd + vcf_filter_arg
      
      bgzip_call = paths['bin_home'] + 'bgzip -c > %s' % o
      
      call = '%s | %s | %s' % (gzcat_call, vcf_filter_call, bgzip_call)
      logger.info(call)
      subprocess.check_call(call, shell=True)
   else:
      call = 'cp %s %s' % (vcfgz, o)
      logger.info(call)
      subprocess.check_call(call, shell=True)

# create the parser
parser = argparse.ArgumentParser(description='''
   Filter bcf and outputs high confidence reference bases
   ''')

# add the arguments
parser.add_argument('--bcf', help='input bcf', default=None)
parser.add_argument('--chr', help='chromosome short name')
parser.add_argument('--chr_id', help='chromosome id to analyze')
parser.add_argument('--d', help='minimum read depth', type=int, default=10)
parser.add_argument('--D', help='maximum read depth', type=int, default=1000000)
parser.add_argument('--Q', help='minimum quality score', type=float, default=40.0)
parser.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser.add_argument('--dbsnp', help='dbsnp file to use (vcf format)', default=None)
parser.add_argument('--rmsk', help='rmsk to use', default=None)
parser.add_argument('--indels', help='indels vcf to read', default=None)
parser.add_argument('--o', help='output vcf.gz [all.ann.vcf.gz]', default='all.ann.vcf.gz')
parser.add_argument('--log', help='log level [INFO]', default='info')

# parse the command line
args = parser.parse_args()
#args = parser.parse_args(' --bcf genotyping/abcalls.all.bcf --chr_id "gi|251831106|ref|NC_012920.1|" --chr chrMT --d 10 --D 50 --Q 40.000000 --ex /panvol1/simon/databases/hs_ref37_rCRS/gi2number.build37_rCRS --dbsnp /panvol1/simon/databases/dbsnp/dbsnp132_hg19.vcf.gz --rmsk /panvol1/simon/databases/hs_ref37_rCRS/rmsk/rmsk_build37_rCRS.number.sort.genome --indels genotyping/indels_for_filtering.number.vcf --o abcalls.chr22.ref.ann.vcf.gz'.split())

logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

files = {}
files['filterAll'] = 'genotyping/tmp.all.bcf.%s.flt.vcf.gz' % args.chr
files['filterAll_tbi'] = 'genotyping/tmp.all.bcf.%s.flt.vcf.gz.tbi' % args.chr
files['dbsnp_ann'] = 'genotyping/tmp.all.bcf.%s.flt.ann.vcf.gz' % args.chr
files['rmsk'] = 'genotyping/tmp.all.bcf.%s.flt.ann.nr.vcf.gz' % args.chr
files['indel_filt'] = 'genotyping/tmp.indel_filtered.%s.vcf' % args.chr

# vcf_filter_All
vcf_filterAll(args.bcf, args.chr_id, args.d, args.D, args.Q, args.ex, files['filterAll'])

# tabix
vcf_tabix(files['filterAll'])

# dbsnp
vcf_annotate_dbsnp(files['filterAll'], args.dbsnp, files['dbsnp_ann'])

# rmsk filtering
if args.chr.find('MT') > -1:
   # if chromosome short name is chrMT or MT run manual filtering for MT only
   manual_rmsk_filter(files['dbsnp_ann'], args.chr, args.rmsk, files['rmsk'])
else:
   # filter for rmsk using BEDtools
   vcf_filter_rmsk(files['dbsnp_ann'], args.rmsk, files['rmsk'])

# indel filter
vcf_filter_indels(files['rmsk'], args.chr, args.indels, files['indel_filt'], args.o)

# remove tmp files
pipelinemod.rm_files(files.values())


