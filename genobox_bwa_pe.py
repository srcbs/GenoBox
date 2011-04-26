#!/usr/bin/python

# python2.7 /panvol1/simon/bin/pipeline/genobox_bwa.py --i data/C6L_1_1.cleaned.truncated --bwaindex hla/A_gen.fasta --outpath hla

import argparse
import pipelinemod
import os
import logging
import subprocess
import re

parser = argparse.ArgumentParser(description='''
   Start BWA PE mapping for each input fqfiles against bwaindex
   ''')

# add the arguments
parser.add_argument('--pe1', help='input fqfiles files PE1', nargs='+', action='append')
parser.add_argument('--pe2', help='input fqfiles files PE2', nargs='+', action='append')
parser.add_argument('--bwaindex', help='bwa indexes to use')
parser.add_argument('--outpath', help='outpath to put sam/bam [bwa_run]', default='bwa_run')
parser.add_argument('--fqtype', help='if fastq, quality type (Sanger, Illumina) [Sanger]', default='Sanger')
parser.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int)
parser.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser.add_argument('--n', help='number of threads for parallel run [16]', default='16')
parser.add_argument('--h', help='do not release jobs [False]', default=False, action='store_true')
parser.add_argument('--q', help='queue to submit jobs to (cbs, urgent) [cbs]', default='cbs')
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --bwaindex /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --fqtype Sanger --n 4'.split())

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
cpuA = 'ncpus=1,mem=512mb'
cpuB = 'nodes=1:ppn=16,mem=10gb'
cpuC = 'ncpus=1,mem=2gb'
cpuE = 'ncpus=1,mem=5gb'
if args.n != '16':
   cpuD = 'nodes=1:ppn=%s,mem=5gb' % args.n
else:
   cpuD = cpuB

# set if args.n = 1
if args.n == '1':
   cpuD = cpuE

if not os.path.exists(args.outpath):
   os.makedirs(args.outpath)

# get filenames of pe1
saifiles1 = []
basenames1 = []
files = args.pe1[0]
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
   basenames1.append(basename)
   saifile = args.outpath + '/'+ basename+'.sai'
   saifiles1.append(saifile)

# get filenames of pe2
saifiles2 = []
basenames2 = []
files = args.pe2[0]
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
   basenames2.append(basename)
   saifile = args.outpath + '/'+ basename+'.sai'
   saifiles2.append(saifile)

# check basenames are the same for both basenames
# implies '_' in name (before 1/2)
# create samfiles and bamfiles if it is accepted
samfiles = []
bamfiles = []
for i,bn in enumerate(basenames1):
   bnshort1 = bn.split('_')[0]
   bnshort2 = basenames2[i].split('_')[0]
   if bnshort1 != bnshort2:
      raise ValueError('Paired end filenames did not match! %s != %s' % (bnshort1, bnshort2))
   else:
      samfile = args.outpath + '/'+ bnshort1+'.sam'
      bamfile = args.outpath + '/'+ bnshort1+'.bam'
      samfiles.append(samfile)
      bamfiles.append(bamfile)
      

## create calls

runcalls = {}

# set quality type
if args.fqtype == 'Sanger':
   bwa_cmd = '%sbwa aln ' % paths['bwa_home']
elif args.fqtype == 'Illumina':
   bwa_cmd = '%sbwa aln -I ' % paths['bwa_home']
else:
   raise ValueError('fqtype must be Sanger or Illumina')

# create individual calls
aligncalls1 = []
aligncalls2 = []
sampecalls = []
for i in range(len(basenames1)):
   aligncall1 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, args.n, args.qtrim, args.bwaindex, saifiles1[i], args.pe1[0][i])
   aligncall2 = '%s -t %s -q %i %s -f %s %s ' % (bwa_cmd, args.n, args.qtrim, args.bwaindex, saifiles2[i], args.pe2[0][i])
   sampecall = '%sbwa sampe -a %i %s -f %s %s %s %s %s ' % (paths['bwa_home'], args.a, args.bwaindex, samfiles[i], saifiles1[i], saifiles2[i], args.pe1[0][i], args.pe2[0][i])
   aligncalls1.append(aligncall1)
   aligncalls2.append(aligncall2)
   sampecalls.append(sampecall)

runcalls['bwa_align1'] = aligncalls1
runcalls['bwa_align2'] = aligncalls2
runcalls['bwa_sampe'] = sampecalls

# sam2bam
sam2bamcalls = []
for i in range(len(basenames1)):
   call = '%ssamtools view -bS -o %s %s' % (paths['samtools_home'], bamfiles[i], samfiles[i])
   sam2bamcalls.append(call)

runcalls['sam2bam'] = sam2bamcalls



## submit jobs
allids = []
bwa_alignids1 = pipelinemod.submitjob(runcalls['bwa_align1'], home, paths, logger, 'run_genobox_bwaalign1', args.q, cpuD, False)
bwa_alignids2 = pipelinemod.submitjob(runcalls['bwa_align2'], home, paths, logger, 'run_genobox_bwaalign2', args.q, cpuD, False)

# set jobids in the correct way
bwa_alignids = []
for i in range(len(bwa_alignids1)):
   bwa_alignids.append(bwa_alignids1[i])
   bwa_alignids.append(bwa_alignids2[i])

bwa_samseids = pipelinemod.submitjob(runcalls['bwa_sampe'], home, paths, logger, 'run_genobox_bwasampe', args.q, cpuE, True, 'conc', len(basenames1), True, *bwa_alignids)
sam2bamids = pipelinemod.submitjob(runcalls['sam2bam'], home, paths, logger, 'run_genobox_sam2bam', args.q, cpuC, True, 'one2one', 1, True, *bwa_samseids)

# release jobs
if args.h:
   pass
else:
   allids.extend(bwa_alignids) ; allids.extend(bwa_samseids) ; allids.extend(sam2bamids)
   releasemsg = pipelinemod.releasejobs(allids)

