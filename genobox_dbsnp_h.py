#!/panvol1/simon/bin/python2.7

from __future__ import division

import argparse
import subprocess
import genobox_modules
import os
import logging

def vcf_exchange_ids(vcf, ex, vcf_out):
   '''Exchange ids in vcf'''
   print 'ex: %s' % ex
   print 'vcf: %s' % vcf
   
   if ex and ex != 'None':
      cmd = paths['genobox_home'] + 'genobox_exchangeids.py'
      if vcf.endswith('.gz'):
         call = 'gzip -dc %s | %s --x 0 --b %s --o %s' % (vcf, cmd, ex, vcf_out)
      else:
         arg = ' --a %s --x 0 --b %s --o %s' % (vcf, ex, vcf_out)
         call = cmd+arg
      logger.info(call)
      subprocess.check_call(call, shell=True)
   else:
      if vcf.endswith('.gz'):
         call = 'gzip -dc %s > %s' % (vcf, vcf_out)
      else:
         call = 'cp %s %s' % (vcf, vcf_out)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def vcf_sortbed(vcf, vcf_out):
   '''Sort vcf using bedtools'''
   
   # add header
   head_call = 'head -n 1000 %s | grep "#" > %s' % (vcf, 'header.vcf')
   logger.info(head_call)
   subprocess.check_call(head_call, shell=True)
   
   # sort
   sort_cmd = paths['bedtools_home'] + 'sortBed'
   sort_arg = ' -i %s | cat header.vcf - > %s' % (vcf, vcf_out)
   sort_call = sort_cmd + sort_arg
   logger.info(sort_call)
   subprocess.check_call(sort_call, shell=True)

def vcf_bgzip_tabix(vcf):
   '''Run bgzip and tabix on vcf'''
   
   paths = genobox_modules.setSystem()
   
   bgzip_call = paths['bin_home'] + 'bgzip -f %s' % vcf
   logger.info(bgzip_call)
   subprocess.check_call(bgzip_call, shell=True)
   
   tabix_call = paths['bin_home'] + 'tabix -p vcf -f %s.gz' % (vcf)
   logger.info(tabix_call)
   subprocess.check_call(tabix_call, shell=True)

def vcf_annotate_dbsnp(vcfgz, dbsnp, vcfgz_out):
   '''Annotate vcf.gz with dbsnp'''
   
   if dbsnp:
      gunzip_call = '/usr/bin/gunzip -c %s' % vcfgz
      fill_call = paths['bin_home'] + 'fill-rsIDs -r %s' % dbsnp
      bgzip_call = paths['bin_home'] + 'bgzip -c > %s' % vcfgz_out
      
      dbsnp_call = '%s | %s | %s' % (gunzip_call, fill_call, bgzip_call)
      logger.info(dbsnp_call)
      subprocess.check_call(dbsnp_call, shell=True)
   else:
      call = 'cp %s %s' % (vcfgz, vcfgz_out)
      logger.info(call)
      subprocess.check_call(call, shell=True)

def write_indels_for_filtering(var_vcf, ex, indel_vcf):
   '''Create indels_for_filtering file '''
   
   import genobox_modules
   import subprocess
   
   paths = genobox_modules.setSystem()
   
   grep_call = 'grep -v \"#\" %s | grep "INDEL" | cat header.vcf - > tmp_file_indels' % (var_vcf)
   logger.info(grep_call)
   subprocess.check_call(grep_call, shell=True)
   
   ex_cmd = paths['genobox_home'] + 'genobox_exchangeids.py'
   ex_arg = ' --a tmp_file_indels --x 0 --b %s --o %s' % (ex, indel_vcf)
   ex_call = cmd + arg
   logger.info(ex_call)
   subprocess.check_call(ex_call, shell=True)


# create the parser
parser = argparse.ArgumentParser(description=
   '''Annotate vcf.gz file with dbSNP,
   exchanging chromsome names to dbSNP version
   sort vcf and the input to dbSNP
   ''')

# add the arguments
parser.add_argument('--vcf', help='input vcf file')
parser.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser.add_argument('--dbsnp', help='dbsnp file to use (vcf format)', default=None)
parser.add_argument('--o', help='output vcf.gz [snpcalls.dbsnp.vcf.gz]', default='snpcalls.vcf.gz')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--vcf'.split())

# set logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

paths = genobox_modules.setSystem()

# exchange ids to numbers
vcf_exchange_ids(args.vcf, args.ex, 'genotyping/tmp.filtered.numbers.vcf')

# sort using bed
vcf_sortbed('genotyping/tmp.filtered.numbers.vcf', 'genotyping/tmp.filtered.numbers.sort.vcf')

# bgzip and index using tabix
vcf_bgzip_tabix('genotyping/tmp.filtered.numbers.sort.vcf')

# annotate with dbsnp
vcf_annotate_dbsnp('genotyping/tmp.filtered.numbers.sort.vcf.gz', args.dbsnp, args.o)

# create indels for filtering
write_indels_for_filtering(args.o, args.ex, 'indels_for_filtering.vcf')
