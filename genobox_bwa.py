#!/usr/bin/python

# python2.7 /panvol1/simon/bin/pipeline/genobox_bwa.py --i data/C6L_1_1.cleaned.truncated --bwaindex hla/A_gen.fasta --outpath hla

import argparse
import pipelinemod
import os
import logging
import subprocess
import re

parser = argparse.ArgumentParser(description='''
   Start BWA mapping for each input fqfile against bwaindex
   ''')

# add the arguments
parser.add_argument('--se', help='input fqfiles files', nargs='+', action='append')
parser.add_argument('--bwaindex', help='bwa indexes to use')
parser.add_argument('--fqtype', help='if fastq, quality type (Sanger, Illumina) [Sanger]', default='Sanger')
parser.add_argument('--outpath', help='outpath to put sam/bam [bwa_run]', default='bwa_run')
parser.add_argument('--n', help='number of threads for parallel run [16]', default='16')
parser.add_argument('--m', help='memory requirements for bwa, samse, sam2bam (if n != 16) [5gb]', default='5gb')
parser.add_argument('--q', help='queue to submit jobs to (cbs, urgent) [cbs]', default='cbs')
parser.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int)
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--i data/C6L_1_1.cleaned.truncated --bwaindex hla/A_gen.fasta --outpath hla'.split())

# set logging
logger = logging.getLogger('genobox_alignment.py')
hdlr = logging.FileHandler('genobox_alignment.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# set paths
paths = pipelinemod.setSystem()
home = os.getcwd()
cpuB = 'nodes=1:ppn=16,mem=10gb'
cpuE = 'ncpus=1,mem=%s' % args.m
if args.n != '16':
   cpuD = 'nodes=1:ppn=%s,mem=%s' % (args.n, args.m)
else:
   cpuD = cpuB

# set if args.n = 1
if args.n == 1:
   cpuD = cpuE


if not os.path.exists(args.outpath):
   os.makedirs(args.outpath)

# get filename
samfiles = []
bamfiles = []
saifiles = []
files = args.se[0]
for file in files:
   path = os.path.split(file)[0]
   filename = os.path.split(file)[1]
   fileprefix = os.path.splitext(filename)[0]
   pattern = re.compile(r'(.+?)(\.[^.]*$|$)')
   match = pattern.search(fileprefix)
   if match:
      basename = match.group(1)
   else:
      basename = fileprefix
   saifile = args.outpath + '/'+ basename+'.sai'
   samfile = args.outpath + '/'+ basename+'.sam'
   bamfile = args.outpath + '/'+ basename+'.bam'
   saifiles.append(saifile)
   samfiles.append(samfile)
   bamfiles.append(bamfile)

## create calls

runcalls = {}

# set quality type
if args.fqtype == 'Sanger':
   bwa_cmd = '%sbwa aln ' % paths['bwa_home']
elif args.fqtype == 'Illumina':
   bwa_cmd = '%sbwa aln -I ' % paths['bwa_home']

# use bwa against ref
bwa_aligncalls = []
bwa_samsecalls = []
for i in range(len(files)):
   bwa_align = '%s -t %s -q %i %s %s > %s' % (bwa_cmd, args.n, args.qtrim, args.bwaindex, files[i], saifiles[i])
   bwa_samse = '%sbwa samse %s %s %s > %s' % (paths['bwa_home'], args.bwaindex, saifiles[i], files[i], samfiles[i])
   bwa_aligncalls.append(bwa_align)
   bwa_samsecalls.append(bwa_samse)

runcalls['bwa_align'] = bwa_aligncalls
runcalls['bwa_samse'] = bwa_samsecalls

# sam2bam
sam2bamcalls = []
for i in range(len(files)):
   file = files[i]
   call = '%ssamtools view -bS -o %s %s' % (paths['samtools_home'], bamfiles[i], samfiles[i])
   sam2bamcalls.append(call)

runcalls['sam2bam'] = sam2bamcalls
## submit jobs

allids = []
if args.n != '16':
   bwa_alignids = pipelinemod.submitjob(runcalls['bwa_align'], home, paths, logger, 'run_genobox_bwaalign', 'cbs', cpuD, False)
   bwa_samseids = pipelinemod.submitjob(runcalls['bwa_samse'], home, paths, logger, 'run_genobox_bwasamse', 'cbs', cpuE, True, 'one2one', 1, True, *bwa_alignids)
else:
   bwa_alignids = pipelinemod.submitjob(runcalls['bwa_align'], home, paths, logger, 'run_genobox_bwaalign', 'cbs', cpuD, False)
   bwa_samseids = pipelinemod.submitjob(runcalls['bwa_samse'], home, paths, logger, 'run_genobox_bwasamse', 'cbs', cpuE, True, 'one2one', 1, True, *bwa_alignids)

sam2bamids = pipelinemod.submitjob(runcalls['sam2bam'], home, paths, logger, 'run_genobox_sam2bam', 'cbs', cpuE, True, 'one2one', 1, True, *bwa_samseids)
allids.extend(bwa_alignids) ; allids.extend(bwa_samseids) ; allids.extend(sam2bamids)

# add semaphore

# release jobs
releasemsg = pipelinemod.releasejobs(allids)
