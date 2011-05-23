#!/usr/bin/python

import argparse
import subprocess
import genobox_modules
import os
import logging


# create the parser
parser = argparse.ArgumentParser(description=
'''Run samtools mpileup command on BAM and chromosome as input,
   BAM file MUST be indexed using samtools index or similar.
   Default is to print out bcf of all covered sites, --var sets it to variants only''')

# add the arguments
parser.add_argument('--bam', help='input file')
parser.add_argument('--chr', help='chromosome to analyze [all]', default='all')
parser.add_argument('--fa', help='reference fasta')
parser.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2)', default=None)
parser.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float)
parser.add_argument('--var', help='variants only, default all [False]', default=False, action='store_true')
parser.add_argument('--o', help='output file')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--bam stinus_d12.q25.nc.bam --chr gi|224589823|ref|NC_000024.9| --fa /panvol1/simon/databases/hs_ref36_ucsc/hg18.fa --o chrY.bcf'.split())

# set logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set paths
paths = genobox_modules.setSystem()
home = os.getcwd()

samcmd = paths['samtools_home']+'samtools'
bcfcmd = paths['samtools_home']+'bcftools'

# set mpileup command
if args.chr == 'all':
   mpileup = '%s mpileup -ugf %s %s' % (samcmd, args.fa, args.bam)
else:
   mpileup = '%s mpileup -ugf %s -r "%s" %s' % (samcmd, args.fa, args.chr, args.bam)

# set to output variants only or all
if args.var:
   bvcg = '-bvcg'
else:
   bvcg = '-bcg'

cmd = ''' %s | %s view %s -P %s -p %f - > %s''' % (mpileup, bcfcmd, bvcg, args.prior, args.pp, args.o)

logger.info(cmd)
exitcode = subprocess.check_call(cmd, shell=True)
